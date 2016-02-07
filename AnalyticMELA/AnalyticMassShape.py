import re
import ROOT
import sys

#==================================
#parameters
decaymode = "Z"
spin = 0

realconstants = {
                 "M_Z": None,
                 "Ga_Z": None,
                 "M_W": None,
                 "Ga_W": None,
                }
complexconstants = {
                    "ghz1": None,
                    "ghz2": None,
                    "ghz4": None,
                    "b1": None,
                    "b5": None,
                   }
#if None, read from
#../JHUgenerator/mod_Parameters.F90

m4lmin = 60
m4lmax = 1000
#==================================

#http://arxiv.org/pdf/1208.4018v3.pdf

def spin0():
    #define RooConstVars
    mV = ROOT.RooConstVar("mV", "mV", realconstants["M_%s"%decaymode])
    gammaV = ROOT.RooConstVar("gammaV", "gammaV", realconstants["Ga_%s"%decaymode])
    Reg1 = ROOT.RooConstVar("Reg1", "Reg1", complexconstants["ghz1"].real)
    Img1 = ROOT.RooConstVar("Img1", "Img1", complexconstants["ghz1"].imag)
    Reg2 = ROOT.RooConstVar("Reg2", "Reg2", complexconstants["ghz2"].real)
    Img2 = ROOT.RooConstVar("Img2", "Img2", complexconstants["ghz2"].imag)
    Reg4 = ROOT.RooConstVar("Reg4", "Reg4", complexconstants["ghz4"].real)
    Img4 = ROOT.RooConstVar("Img4", "Img4", complexconstants["ghz4"].imag)

    #define RooRealVars
    m4l = ROOT.RooRealVar("m4l", "m4l", m4lmin, m4lmin, m4lmax)
    m1 = ROOT.RooRealVar("m1", "m1", realconstants["M_%s"%decaymode], 0, 200)
    m2 = ROOT.RooRealVar("m2", "m2", realconstants["M_%s"%decaymode], 0, 200)



    #eq. (13)
    s = ROOT.RooFormulaVar("s", "s", "(@0**2-@1**2-@2**2)/2", ROOT.RooArgList(m4l, m1, m2))

    #eq. (15)
    x = ROOT.RooFormulaVar("x", "x", "(@0/(@1*@2))**2-1", ROOT.RooArgList(s, m1, m2))

    #eq. (12)
    Rea1 = ROOT.RooFormulaVar("Rea1", "Rea1", "@0*(@1/@2)**2 + 2*@3*@4/@2**2", ROOT.RooArgList(Reg1, mV, m4l, Reg2, s, m4l))
    Ima1 = ROOT.RooFormulaVar("Ima1", "Ima1", "@0*(@1/@2)**2 + 2*@3*@4/@2**2", ROOT.RooArgList(Img1, mV, m4l, Img2, s, m4l))
    Rea2 = ROOT.RooFormulaVar("Rea2", "Rea2", "-2*@0", ROOT.RooArgList(Reg2))
    Ima2 = ROOT.RooFormulaVar("Ima2", "Ima2", "-2*@0", ROOT.RooArgList(Img2))
    Rea3 = ROOT.RooFormulaVar("Rea3", "Rea3", "-2*@0", ROOT.RooArgList(Reg4))
    Ima3 = ROOT.RooFormulaVar("Ima3", "Ima3", "-2*@0", ROOT.RooArgList(Img4))

    #eq. (14)
    ReA00 = ROOT.RooFormulaVar("ReA00", "ReA00", "-@0**2 * (@1*sqrt(1+@2) + @3*(@4*@5/@0**2)*@2)      ", ROOT.RooArgList(m4l, Rea1, x, Rea2, m1, m2))
    ImA00 = ROOT.RooFormulaVar("ImA00", "ImA00", "-@0**2 * (@1*sqrt(1+@2) + @3*(@4*@5/@0**2)*@2)      ", ROOT.RooArgList(m4l, Ima1, x, Ima2, m1, m2))
    ReApp = ROOT.RooFormulaVar("ReApp", "ReApp", " @0**2 * (@1            - @3*(@4*@5/@0**2)*sqrt(@2))", ROOT.RooArgList(m4l, Rea1, x, Ima3, m1, m2))
    ImApp = ROOT.RooFormulaVar("ImApp", "ImApp", " @0**2 * (@1            + @3*(@4*@5/@0**2)*sqrt(@2))", ROOT.RooArgList(m4l, Ima1, x, Rea3, m1, m2))
    ReAmm = ROOT.RooFormulaVar("ReAmm", "ReAmm", " @0**2 * (@1            + @3*(@4*@5/@0**2)*sqrt(@2))", ROOT.RooArgList(m4l, Rea1, x, Ima3, m1, m2))
    ImAmm = ROOT.RooFormulaVar("ImAmm", "ImAmm", " @0**2 * (@1            - @3*(@4*@5/@0**2)*sqrt(@2))", ROOT.RooArgList(m4l, Ima1, x, Rea3, m1, m2))

    #Phase space, eq. (6)
    P = ROOT.RooFormulaVar("P", "P",
                             "  sqrt(1 - (@0+@1)**2/@2**2)"
                             "* sqrt(1 - (@0-@1)**2/@2**2)"
                             "* @0**3/((@0**2-@3**2)**2 + @3**2 * @4**2)"
                             "* @1**3/((@1**2-@3**2)**2 + @3**2 * @4**2)",
                           ROOT.RooArgList(m1, m2, m4l, mV, gammaV),
                          )

    #final distribution, eq. (7)
    pdf = ROOT.RooGenericPdf("pdf", "pdf", "(@0**2 + @1**2 + @2**2 + @3**2 + @4**2 + @5**2) * @6", ROOT.RooArgList(ReA00, ImA00, ReApp, ImApp, ReAmm, ImAmm, P))

    frame = m4l.frame()
    c1 = ROOT.TCanvas.MakeDefCanvas()
    pdf.createProjection(ROOT.RooArgSet(m1, m2)).plotOn(frame)
    frame.Draw()
    c1.SaveAs("AnalyticMassShape.png")
    c1.SaveAs("AnalyticMassShape.eps")
    c1.SaveAs("AnalyticMassShape.pdf")
    c1.SaveAs("AnalyticMassShape.root")


def spin2():
    #define RooConstVars
    mV = ROOT.RooConstVar("mV", "mV", realconstants["M_%s"%decaymode])
    gammaV = ROOT.RooConstVar("gammaV", "gammaV", realconstants["Ga_%s"%decaymode])
    Reb1 = ROOT.RooConstVar("Reb1", "Reb1", complexconstants["b1"].real)
    Imb1 = ROOT.RooConstVar("Imb1", "Imb1", complexconstants["b1"].imag)
    Reb5 = ROOT.RooConstVar("Reb5", "Reb5", complexconstants["b5"].real)
    Imb5 = ROOT.RooConstVar("Imb5", "Imb5", complexconstants["b5"].imag)

    #define RooRealVars
    m4l = ROOT.RooRealVar("m4l", "m4l", m4lmin, m4lmin, m4lmax)
    m1 = ROOT.RooRealVar("m1", "m1", realconstants["M_%s"%decaymode], 4, 200)
    m2 = ROOT.RooRealVar("m2", "m2", realconstants["M_%s"%decaymode], 4, 200)



    #eq. (13)
    s = ROOT.RooFormulaVar("s", "s", "(@0**2-@1**2-@2**2)/2", ROOT.RooArgList(m4l, m1, m2))

    #eq. (15)
    x = ROOT.RooFormulaVar("x", "x", "(@0/(@1*@2))**2-1", ROOT.RooArgList(s, m1, m2))

    #eq. (20)
    Rec1 = ROOT.RooFormulaVar("Rec1", "Rec1", "2*@0 + 2*@1*@2**2/@3", ROOT.RooArgList(Reb1, Reb5, mV, s))
    Imc1 = ROOT.RooFormulaVar("Imc1", "Imc1", "2*@0 + 2*@1*@2**2/@3", ROOT.RooArgList(Imb1, Imb5, mV, s))
    Rec2 = ROOT.RooFormulaVar("Rec2", "Rec2", "-@0/2", ROOT.RooArgList(Reb1))
    Imc2 = ROOT.RooFormulaVar("Imc2", "Imc2", "-@0/2", ROOT.RooArgList(Imb1))
    Rec41 = ROOT.RooFormulaVar("Rec41", "Rec41", "-@0", ROOT.RooArgList(Reb1))
    Imc41 = ROOT.RooFormulaVar("Imc41", "Imc41", "-@0", ROOT.RooArgList(Imb1))
    Rec42 = ROOT.RooFormulaVar("Rec42", "Rec42", "-@0", ROOT.RooArgList(Reb1))
    Imc42 = ROOT.RooFormulaVar("Imc42", "Imc42", "-@0", ROOT.RooArgList(Imb1))

    #eq. (14)
    A00string = """

      m4l**4/(m1*m2*sqrt(6)) * c1/8
    + m1*m2/sqrt(6) * (
        c1 * 0.5 * (1+x) - c2*2*x + c41*2*x + c42*2*x
      )
    - (m1**4+m2**4)/(m1*m2*sqrt(6)) * c1/4
    + m1*m2*(m1**2-m2**2)/(m4l**2*sqrt(6)) * (c41-c42)*2*x
    + (m1**8+m2**8)/(m4l**4*m1*m2*sqrt(6)) * c1/8
    + m1**3*m2**3/(m4l**4*sqrt(6)) * (
        c1*(3/4+x) - c2*(4*x+8*x**2)
      )
    + m1*m2*(m1**4+m2**4)/(m4l**4*sqrt(6)) * (
        -c1*0.5*(1+x)+c2*2*x
      )
    """

    Appstring = """
      m4l**2/sqrt(6) * c1/4
    - (m1**4+m2**4)/(m4l**2*sqrt(6)) * c1/4
    + m1**2*m2**2/(m4l**2*sqrt(6)) * (
        c1*(.5+x) + c2*8*x
      )
    """
    Ammstring = Appstring.replace("i", "-i")

    Ap0string = """
      m4l**3/(m2*sqrt(2)) * c1/8
    + m4l*m2/sqrt(2) * (1-m1**2/m2**2) * c1/8
    + m1**2*m2/(m4l*sqrt(2)) * (
        c1 * (
            .25 + .5*x - m2**2/(8*m1**2) - m1**2/(8*m2**2)
        )
      + c41*2*x
      )
    + m1**2*m2**3/(m4l**3*sqrt(2))*c1 * (
        .125 * (m1**4/m2**4 - m2**2/m1**2)
      + (1-m1**2/m2**2) * (.375+.5*x)
      )
    """

    A0pstring = (Ap0string.replace("c41", "c42")
                          .replace("m1", "TMPTMP")
                          .replace("m2", "m1")
                          .replace("TMPTMP", "m2")
                )

    Am0string = Ap0string.replace("i", "-i")
    A0mstring = A0pstring.replace("i", "-i")

    Apmstring = Ampstring = """
      m4l**2*c1/4
    + m1**2*m2**2/m4l**2 * c1*x
    - (m1**2-m2**2)**2/m4l**2 * c1/4
    """

    realvars = [m4l, m1, m2, s, x]
    complexvars = [(Rec1, Imc1), (Rec2, Imc2), (Rec41, Imc41), (Rec42, Imc42)]
    ReA = {}
    ImA = {}
    ReA[0,0], ImA[0,0] = string2RooFormulaVar("A00", "A00", A00string, realvars, complexvars)
    ReA[1,1], ImA[1,1] = string2RooFormulaVar("App", "App", Appstring, realvars, complexvars)
    ReA[-1,-1], ImA[-1,-1] = string2RooFormulaVar("Amm", "Amm", Ammstring, realvars, complexvars)
    ReA[1,0], ImA[1,0] = string2RooFormulaVar("Ap0", "Ap0", Ap0string, realvars, complexvars)
    ReA[0,1], ImA[0,1] = string2RooFormulaVar("A0p", "A0p", A0pstring, realvars, complexvars)
    ReA[-1,0], ImA[-1,0] = string2RooFormulaVar("Am0", "Am0", Am0string, realvars, complexvars)
    ReA[0,-1], ImA[0,-1] = string2RooFormulaVar("A0m", "A0m", A0mstring, realvars, complexvars)
    ReA[1,-1], ImA[1,-1] = string2RooFormulaVar("Apm", "Apm", Apmstring, realvars, complexvars)
    ReA[-1,1], ImA[-1,1] = string2RooFormulaVar("Amp", "Amp", Ampstring, realvars, complexvars)

    #Phase space, eq. (6)
    P = ROOT.RooFormulaVar("P", "P",
                             "  sqrt(1 - (@0+@1)**2/@2**2)"
                             "* sqrt(1 - (@0-@1)**2/@2**2)"
                             "* @0**3/((@0**2-@3**2)**2 + @3**2 * @4**2)"
                             "* @1**3/((@1**2-@3**2)**2 + @3**2 * @4**2)",
                           ROOT.RooArgList(m1, m2, m4l, mV, gammaV),
                          )

    #final distribution, eq. (7)
    tlist = ROOT.TList()
    pdfformula = "("
    i = 0
    for amplitude in ReA.values() + ImA.values():
        tlist.Add(amplitude)
        if i > 0:
            pdfformula += " + "
        pdfformula += "@%i**2" % i
        i += 1
    tlist.Add(P)
    pdfformula += ") * %i" % i
    pdf = ROOT.RooGenericPdf("pdf", "pdf", pdfformula, ROOT.RooArgList(tlist))

    frame = m4l.frame()
    c1 = ROOT.TCanvas.MakeDefCanvas()
    pdf.createProjection(ROOT.RooArgSet(m1, m2)).plotOn(frame)
    frame.Draw()
    c1.SaveAs("AnalyticMassShape_spin2.png")
    c1.SaveAs("AnalyticMassShape_spin2.eps")
    c1.SaveAs("AnalyticMassShape_spin2.pdf")
    c1.SaveAs("AnalyticMassShape_spin2.root")


def string2RooFormulaVar(name, title, formula, realvars, complexvars):
    formula = formula.replace("\n", "")
    args = []
    index = 0
    for var in realvars:
        regex = r"\b" + var.GetName() + r"\b"
        if re.search(regex, formula):
            formula = re.sub(regex, "@%i"%len(args), formula)
            args.append(var)

    reformula = formula
    imformula = formula
    reargs = args[:]
    imargs = args[:]
    for revar, imvar in complexvars:
        rename, imname = revar.GetName(), imvar.GetName()
        varname = rename.replace("Re", "")
        assert varname == imname.replace("Im", "")
        regex = r"\b" + varname + r"\b"
        regexi = r"\bi[*]" + varname + r"\b"
        if re.search(regexi, reformula):
            reformula = re.sub(regexi, "@%i"%len(reargs), reformula)
            reargs.append(imvar)
            imformula = re.sub(regexi, "(-@%i)"%len(imargs), imformula)
            imargs.append(revar)
        if re.search(regex, reformula):
            reformula = re.sub(regex, "@%i"%len(reargs), reformula)
            reargs.append(revar)
            imformula = re.sub(regex, "@%i"%len(imargs), imformula)
            imargs.append(imvar)

    rearglist = ROOT.RooArgList(*reargs)
    imarglist = ROOT.RooArgList(*imargs)

    return (ROOT.RooFormulaVar("Re"+name, "Re"+title, reformula, rearglist),
            ROOT.RooFormulaVar("Im"+name, "Im"+title, imformula, imarglist))

def getparameters():
    global decaymode, spin, realconstants, complexconstants

    with open("../JHUGenerator/mod_Parameters.F90") as f:
        for line in f:
            line = line.split("!")[0]
            for constantname in realconstants:
                if realconstants[constantname] is None:
                    regex = r"\W" + constantname + r"\s*=([^!*]*)"
                    result = re.search(regex, line)
                    if result:
                        value = result.group(1)
                        value = value.replace("d", "e").replace("D", "E")  #fortran format --> normal format
                        realconstants[constantname] = float(value)
            for constantname in complexconstants:
                if complexconstants[constantname] is None:
                    regex = r"\W" + constantname + (r"\s*=\s*[(]"
                                                     "([^,]*)"     #real part
                                                     ","
                                                     "([^)]*)"   ) #imaginary part
                    result = re.search(regex, line)
                    if result:
                        real = result.group(1).replace("d","e").replace("D","E")
                        imag = result.group(2).replace("d","e").replace("D","E")
                        complexconstants[constantname] = float(real)+1j*float(imag)

    for arg in sys.argv[1:]:
        success = False
        for constantname in realconstants:
            if arg.split("=")[0] == constantname:
                try:
                    realconstants[constantname] = float(arg.split("=")[1])
                except (ValueError, IndexError):
                    raise ValueError("Syntax: " + constantname + "=123.456")
                success = True
        for constantname in complexconstants:
            if arg.split("=")[0] == constantname:
                try:
                    Re = float(arg.split("=")[1].split(",")[0])
                    Im = float(arg.split("=")[1].split(",")[1])
                    complexconstants[constantname] = Re + Im*1j
                except (ValueError, IndexError):
                    raise ValueError("Syntax: " + constantname + "=1.23,4.56")
                success = True
        if arg.split("=")[0] == "decaymode":
            decaymode = arg.split("=")[1]
            success = True
        if arg.split("=")[0] == "spin":
            spin = int(arg.split("=")[1])
            success = True
        if not success:
            raise ValueError("Unknown argument: " + arg)

    for constants in realconstants, complexconstants:
        for constantname in constants:
            if constants[constantname] is None:
                raise IOError("Value of %s not found in mod_Parameters and not specified!" % constantname)


    if decaymode not in ("Z", "W"):
        raise ValueError("Invalid decay mode!  Needs to be Z or W.")

if __name__ == '__main__':
    getparameters()
    if spin == 0:
        spin0()
    elif spin == 2:
        spin2()
    else:
        raise ValueError("Spin must be 0 or 2!")
