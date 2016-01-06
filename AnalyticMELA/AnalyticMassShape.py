import ROOT

#============================
#parameters
decaymode = "Z"
constants = {
             "M_Reso": None,
             "Ga_Reso": None,
             "M_Z": None,
             "Ga_Z": None,
             "M_W": None,
             "Ga_W": None,
             "vev": 246,
             "ghz1": 2,
             "ghz2": 0,
             "ghz4": 0,
            }
#if None, read from
#../JHUgen/mod_Parameters.F90
#============================

#http://arxiv.org/pdf/1208.4018v3.pdf

if decaymode not in ("Z", "W"):
    raise ValueError("Invalid decay mode!  Needs to be Z or W.")
if constants["ghz1"] != 2 or constants["ghz2"] != 0 or constants["ghz4"] != 0:
    raise NotImplementedError("anomalous couplings not yet implemented!")


with open("../JHUGenerator/mod_Parameters.F90") as f:
    for line in f:
        line = line.split("!")[0]
        for constantname in constants:
            if constants[constantname] is None and constantname in line:
                line = line.split("=")[1].split("*")[0]          #isolate the number
                line = line.replace("d", "e").replace("D", "E")  #fortran format --> normal format
                constants[constantname] = float(line)

for constantname in constants:
    if constants[constantname] is None:
        raise IOError("Value of %s not found in mod_Parameters!" % constantname)


vev = ROOT.RooConstVar("vev", "vev", constants["vev"])
mV = ROOT.RooConstVar("mV", "mV", constants["M_%s"%decaymode])
mH = ROOT.RooConstVar("mH", "mH", constants["M_Reso"])
gammaV = ROOT.RooConstVar("gammaV", "gammaV", constants["Ga_%s"%decaymode])
gammaH = ROOT.RooConstVar("gammaH", "gammaH", constants["Ga_Reso"])
a1 = ROOT.RooConstVar("a1", "a1", constants["ghz1"])
a2 = ROOT.RooConstVar("a1", "a1", constants["ghz2"])
a3 = ROOT.RooConstVar("a1", "a1", constants["ghz4"])

m4l = ROOT.RooRealVar("m4l", "m4l", constants["M_Reso"], 0, 2000)
m1 = ROOT.RooRealVar("m1", "m1", constants["M_%s"%decaymode], 0, 200)
m2 = ROOT.RooRealVar("m2", "m2", constants["M_%s"%decaymode], 0, 200)

x = ROOT.RooFormulaVar("x", "x", "((@0**2-@1**2-@2**2)/(2*@1*@2))**2-1", ROOT.RooArgList(mH, m1, m2))  #eq. (15)

#Anomalous couplings will go here, eq. (14)
A00 = ROOT.RooFormulaVar("A00", "A00", "-@0**2/@1 * (@2*sqrt(1+@3))", ROOT.RooArgList(m4l, vev, a1, x))
App = ROOT.RooFormulaVar("App", "App", "@0**2/@1 * @2", ROOT.RooArgList(m4l, vev, a1))
Amm = ROOT.RooFormulaVar("Amm", "Amm", "@0**2/@1 * @2", ROOT.RooArgList(m4l, vev, a1))

#Phase space, eq. (6)
P = ROOT.RooFormulaVar("P", "P",
                         "  sqrt(1 - (@0+@1)**2/@2**2)"
                         "* sqrt(1 - (@0-@1)**2/@2**2)"
                         "* @0**3/((@0**2-@3**2)**2 + @3**2 * @4**2)"
                         "* @1**3/((@1**2-@3**2)**2 + @3**2 * @4**2)",
                       ROOT.RooArgList(m1, m2, m4l, mV, gammaV),
                      )

#final distribution, eq. (7)
pdf = ROOT.RooGenericPdf("pdf", "pdf", "(@1**2 + @2**2 + @3**2) * @4", ROOT.RooArgList(A00, App, Amm, P))

frame = m4l.frame()
c1 = ROOT.TCanvas.MakeDefCanvas()
pdf.createProjection(ROOT.RooArgSet(m1, m2)).plotOn(frame)
frame.Draw()
c1.SaveAs("AnalyticMassShape.png")
