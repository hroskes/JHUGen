from abc import ABCMeta, abstractmethod, abstractproperty
import errno
import os
import pipes
import subprocess

def call(function):
    return function()

def mkdir_p(path):
    """http://stackoverflow.com/a/600612/5228524"""
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

@call
def masses():
    result = []
    for m in range(70, 500, 2): result.append(m)
    for m in range(500, 1000, 5): result.append(m)
    for m in range(1000, 1500, 10): result.append(m)
    for m in range(1500, 3000, 50): result.append(m)
    result.append(125)
    result.append(3000)
    result.sort()
    return result

@call
def allarguments():
    result = set()
    with open("main.F90") as f:
        for line in f:
            line = line.split("!")[0]
            if "ReadCommandLineArgument" in line:
                split = line.split('"')
                assert len(split) == 3
                result.add(split[1])
    return result

@call
def couplingnames():
    result = set()
    for _ in allarguments:
        if ("ghg" in _ or "ghz" in _ or "ghw" in _ or "zprime_" in _ or "graviton_" in _ or "kappa" in _
         or any(_ == "a{}".format(i) for i in range(1,6)) or any(_ == "b{}".format(i) for i in range(1, 11))):
            result.add(_)
    return result

def submit_job(outdir, jobbasename, MReso, **kwargs):
    MReso = float(MReso)
    if MReso.is_integer():
        MReso = int(MReso)
    MReso = str(MReso)

    for k, v in list(kwargs.iteritems()):
        if k not in allarguments:
            raise TypeError("Unknown kwarg {}={}!".format(k, v))
        if k in couplingnames:
            if isinstance(v, basestring):
                v = complex(v.replace("i", "j"))
            kwargs[k] = "{},{}".format(v.real, v.imag)

    otherargs = " ".join(pipes.quote("{}={}".format(k, v)) for k, v in sorted(kwargs.iteritems()))

    mkdir_p(outdir)
    outfile = os.path.join(os.getcwd(), outdir, "m{}.out".format(MReso))
    if os.path.exists(outfile):
        return
    jobname = "{}_{}".format(jobbasename, MReso)
    bjobs = subprocess.check_output(["bjobs"])
    if jobname in bjobs.split():
        return

    subprocess.check_call(["sbatch", "-o", outfile, "-e", outfile, "--job-name", jobname, "template.slurm.sh", str(MReso), otherargs])

def submit_masses(outdir, jobbasename, masses, **kwargs):
    for mass in masses:
        submit_job(outdir, jobbasename, mass, **kwargs)

class XsecScanBase(object):
    __metaclass__ = ABCMeta
    @abstractproperty
    def outdir(self):
        pass
    @abstractproperty
    def jobbasename(self):
        pass
    @abstractproperty
    def kwargs(self):
        pass
    def submit_masses(self):
        submit_masses(self.outdir, self.jobbasename, masses, **self.kwargs)

class XsecScanProcessCouplings(XsecScanBase):
    def __init__(self, process, coupling, usecut=False):
        self.process = process
        self.coupling = coupling
        self.usecut = usecut
        if coupling not in self.allowedcouplings(process):
            raise ValueError("Bad coupling {} for process {}".format(coupling, process))
    @property
    def kwargs(self):
        result = {"Process": self.processid}
        result.update(self.couplings)
        result.update(self.decaymodes)
        result.update(self.cut)
        return result
    @property
    def couplings(self):
        if self.processid == 61:
            if self.coupling == "SM": return {"ghg2": 1}
            elif self.coupling == "g4": return {"ghg2": 0, "ghg4": 1}
        else:
            if self.coupling == "SM": return {"ghz1": 1}
            elif self.coupling == "g2": return {"ghz1": 0, "ghz2": 1}
            elif self.coupling == "g4": return {"ghz1": 0, "ghz4": 1}
            elif self.coupling == "L1": return {"ghz1": 0, "ghz1_prime2": 1}
        assert False
    @property
    def processid(self):
        if self.process == "VBF": return 60
        elif self.process == "HJJ": return 61
        elif self.process == "ZH": return 50
        elif self.process == "WH": return 50
        elif self.process == "HZZ2e2mu": return 0
        assert False
    @property
    def decaymodes(self):
        if self.process in ["VBF", "HJJ"]: return {}
        if self.process == "ZH": return {"DecayMode1": 9}
        if self.process == "WH": return {"DecayMode1": 11}
        if self.process == "HZZ2e2mu": return {"DecayMode1": 0, "DecayMode2": 0, "Interf": 0}
        assert False
    @property
    def cut(self):
        if self.process == "HZZ" and self.usecut: raise ValueError("Can't run {} with cut!".format(self.process))
        result = {}
        if self.process == "VBF": result.update(pTjetcut=0, deltaRcut=0)
        if self.usecut: result.update(pTjetcut=25, etajetcut=4.7, pTlepcut=5, etalepcut=2.5)
        return result
    @property
    def jobbasename(self):
        return self.outdir.replace("/", "_")
    @property
    def outdir(self):
        return os.path.join(self.process, self.coupling + ("_withcut" if self.usecut else ""))

    @staticmethod
    def allowedcouplings(process):
        if process in ["HZZ2e2mu", "VBF", "ZH", "WH"]:
            return ["SM" ,"g2", "g4", "L1"]
        if process == "HJJ":
            return ["SM", "g4"]
        raise ValueError("Bad process {}!".format(process))

    @staticmethod
    def cancut(process):
        if process == "HZZ2e2mu": return False
        if process in ["HJJ", "VBF", "ZH", "WH"]: return True
        assert False

if __name__ == "__main__":
    for process in "HZZ2e2mu", "VBF", "HJJ", "ZH", "WH":
        for couplings in XsecScanProcessCouplings.allowedcouplings(process):
            XsecScanProcessCouplings(process, couplings, False).submit_masses()
            if XsecScanProcessCouplings.cancut(process):
                XsecScanProcessCouplings(process, couplings, True).submit_masses()

