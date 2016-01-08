import re
import ROOT

#==================================
#parameters
decaymode = "Z"

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
                   }
#if None, read from
#../JHUgenerator/mod_Parameters.F90

m4lmin = 0
m4lmax = 1000
#==================================

#http://arxiv.org/pdf/1208.4018v3.pdf

if decaymode not in ("Z", "W"):
    raise ValueError("Invalid decay mode!  Needs to be Z or W.")


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

for constants in realconstants, complexconstants:
    for constantname in constants:
        if constants[constantname] is None:
            raise IOError("Value of %s not found in mod_Parameters!" % constantname)

#define RooConstVars
mV = ROOT.RooConstVar("mV", "mV", realconstants["M_%s"%decaymode])
gammaV = ROOT.RooConstVar("gammaV", "gammaV", realconstants["Ga_%s"%decaymode])
Reg1 = ROOT.RooConstVar("Reg1", "Reg1", complexconstants["ghz1"].real)
Img1 = ROOT.RooConstVar("Img1", "Img1", complexconstants["ghz1"].imag)
Reg2 = ROOT.RooConstVar("Reg2", "Reg2", complexconstants["ghz2"].real)
Img2 = ROOT.RooConstVar("Img2", "Img2", complexconstants["ghz2"].imag)
Reg4 = ROOT.RooConstVar("Reg3", "Reg3", complexconstants["ghz4"].real)
Img4 = ROOT.RooConstVar("Img3", "Img3", complexconstants["ghz4"].imag)

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
