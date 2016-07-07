import ROOT

fspin0 = ROOT.TFile.Open("AnalyticMassShape_workspace.root")
wspin0 = fspin0.Get("workspace")
pdfspin0 = wspin0.pdf("spin0pdf_Proj[m1,m2]")
m4lspin0 = wspin0.var("m4l")
fspin2 = ROOT.TFile.Open("AnalyticMassShape_spin2_workspace.root")
wspin2 = fspin2.Get("workspace")
pdfspin2 = wspin2.pdf("spin2pdf_Proj[m1,m2]")
m4lspin2 = wspin2.var("m4l")

print pdfspin0
print pdfspin2

argsetspin0 = ROOT.RooArgSet(m4lspin0)
argsetspin2 = ROOT.RooArgSet(m4lspin2)

for i in range(4, 100):
    pdfspin0.getVal(argsetspin0)
    pdfspin2.getVal(argsetspin2)
    bincenter = (i+.5)*15
    m4lspin0.setVal(bincenter)
    m4lspin2.setVal(bincenter)
    print pdfspin2.getVal(argsetspin2)/pdfspin0.getVal(argsetspin0), ","
