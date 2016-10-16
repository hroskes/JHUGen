import os
import ROOT
import style
from array import array
from collections import OrderedDict
def function():
        for a in "VBF", "HJJ", "ZH", "WH", "spin2_gg":
                    glist = []
                    mg = ROOT.TMultiGraph("mg", "mg")
                    l = ROOT.TLegend(.6, .7, .9, .9)
                    l.SetFillStyle(0)
                    l.SetBorderSize(0)
                    color = 1
                    for b in "SM", "g2", "g4", "L1", "2bplus", "SM_withcut", "g2_withcut", "g4_withcut", "L1_withcut":
                        if a == "HJJ" and (b == "g2" or b == "L1") or ((a == "spin2_gg") != (b == "2bplus")):
                                    continue
                        dictarea, dicterror =  mass(a,b)
                        if not dictarea: print a, b; continue
                        x = array("d", dictarea.keys())
                        y = array("d", dictarea.values())
                        errx = array ("d", [0] * len(dicterror))
                        erry = array("d", dicterror.values())
                        g = ROOT.TGraphErrors(len( x ), x, y,errx,erry)
                        g.SetMarkerColor(color)
                        g.SetLineColor(color)
                        color += 1
                        l.AddEntry(g,b,"lp")
                        c = ROOT.TCanvas()
                        g.Draw("APEZ")
                        mg.Add(g)
                        glist.append(g)
                        d = a + b + ".png"
                        e = a + b + ".root"
                        f = a + b + ".pdf"
                        c.SaveAs(d)
                        tfile = ROOT.TFile(e, "recreate")
                        g.Write()
                        c.SaveAs(f)
                    c = ROOT.TCanvas()
                    mg.Draw("APEZ")
                    l.Draw()
                    d = a + ".png"
                    e = a + ".root"
                    f = a + ".pdf"
                    c.SaveAs(d)
                    c.SaveAs(e)
                    c.SaveAs(f)

def hmass():
        masses = [125]
        HMASS = 70
        while True:
           if HMASS > 3000:
              break
           masses.append(HMASS) 
           if HMASS < 500:
              HMASS += 2
           elif HMASS < 1000:
              HMASS += 5
           elif HMASS <  1500:
              HMASS += 10
           else:
              HMASS += 50
        return sorted(masses)

def mass(a,b):
        
        dictarea = OrderedDict()
        dicterror = OrderedDict()
        Hmasses = hmass()
        for m in Hmasses:
	    area = 0
            width = 0
            num  = 0
            denom = 0
            filename = "{}/{}/m{}.out".format(a, b, m)
            try:
              f = open(filename) 
            except IOError:
              continue
            for line in f:
                 if "integral" in line:
                     area  = float(line.split()[3])
                 if "std. dev" in line:
                     width  = float(line.split()[4])
                     num += area/(width*width)
                     denom += 1/(width*width)
            dictarea[m] = num/denom
            dicterror[m] = 1/(denom**0.5)
        #area125 = dictarea[125]
        #for m in Hmasses:
        #   dictarea[m] /= area125
        #   dicterror[m] /= area125
        return dictarea, dicterror

function()    
