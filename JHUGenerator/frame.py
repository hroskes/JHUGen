import ROOT
import style
from array import array
def function():
        for a in "VBF", "HJJ", "ZH", "WH":
                    for b in "SM", "g2", "g4", "L1":
                        if a == "HJJ" and (b == "g2" or b == "L1"):
                                    continue
                        dictarea, dicterror =  mass(a,b)
                        x = array("d", dictarea.keys())
                        y = array("d", dictarea.values())
                        errx = array ("d", [0] * len(dicterror))
                        erry = array("d", dicterror.values())
                        g = ROOT.TGraphErrors(len( x ), x, y,errx,erry)
                        c = ROOT.TCanvas()
                        g.Draw("APEZ")
                        d = a + b + ".png"
                        e = a + b + ".root"
                        f = a + b + ".pdf"
                        c.SaveAs(d)
                        c.SaveAs(e)
                        c.SaveAs(f)

def hmass():
        masses = []
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
        return masses

def mass(a,b):
        
        dictarea = {}
        dicterror = {}
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
        return dictarea, dicterror

function()    
