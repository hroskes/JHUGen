import ROOT
from array import array
def function():
        for a in "VBF", "HJJ", "ZH", "WH":
                    for b in "SM", "g2", "g4", "L1":
                        if a == "HJJ" and (b == "g2" or b == "L1"):
                                    continue
                        dictarea, dictwidth =  mass(a,b)
                        x = array("d", dictarea.keys())
                        y = array("d", dictarea.values())
                        g = ROOT.TGraph(len( x ), x, y)
                        c = ROOT.TCanvas()
                        g.Draw("AP")
                        d = a + b +".png"
                        c.SaveAs(d)
def mass(a,b):
        dictarea = {}
        dictwidth = {}
        for m in [70]:
            area = 0
            width = 0
            filename = "{}/{}/m{}.out".format(a, b, m)
            for line in open(filename):
                 if "accum. integral" in line:
                     line = line.replace("*","")
                     area = float(line.split()[7])
                 if "accum. std. dev" in line:
                     line =  line.replace("*","")
                     width  = float(line.split()[9])
            dictarea[m] = area
            dictwidth[m] = width
        return dictarea, dictwidth

function()    
