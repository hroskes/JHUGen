#!/usr/bin/env python

import argparse, array, collections, itertools, os

p = argparse.ArgumentParser()
p.add_argument("infile")
p.add_argument("outfile")
p.add_argument("--vbf", action="store_true")
g = p.add_mutually_exclusive_group()
g.add_argument("--zh")
g.add_argument("--wh")
p.add_argument("--decay", action="store_true")
args = p.parse_args()

if os.path.exists(args.outfile):
  p.error(args.outfile+" already exists")
if not os.path.exists(args.infile):
  p.error(args.infile+" doesn't exist")

if not any((args.vbf, args.zh, args.wh, args.decay)):
  p.error("Have to specify at least one of --vbf --zh --wh --decay")

import ROOT

import lhefile

newf = ROOT.TFile(args.outfile, "CREATE")
t = ROOT.TTree("SelectedTree", "SelectedTree")

branches = collections.OrderedDict()

if args.decay:
  branches.update(collections.OrderedDict(
    (name, array.array("i", [-99999])) for name in (
      "GenLep1Id", "GenLep2Id", "GenLep3Id", "GenLep4Id"
    )
  ))

  branches.update(collections.OrderedDict(
    (name, array.array("f", [-99999])) for name in (
      "GenHMass", "Gencosthetastar", "GenhelcosthetaZ1", "GenhelcosthetaZ2", "Genhelphi", "GenphistarZ1", "GenZ1Mass", "GenZ2Mass",
    )
  ))

if args.vbf:
  branches.update(collections.OrderedDict(
    (name, array.array("f", [-99999])) for name in (
      "Gencosthetastar_VBF", "GenhelcosthetaV1_VBF", "GenhelcosthetaV2_VBF", "Genhelphi_VBF", "GenphistarV1_VBF", "GenQ_V1", "GenQ_V2", "GenDiJetMass", "GenDRjet"
    )
  ))

if args.zh or args.wh:
  branches.update(collections.OrderedDict(
    (name, array.array("f", [-99999])) for name in (
      "Gencosthetastar_VH", "GenhelcosthetaV1_VH", "GenhelcosthetaV2_VH", "Genhelphi_VH", "GenphistarV1_VH",
    )
  ))

for name, branch in branches.iteritems():
  t.Branch(name, branch, name+"/"+{"i": "I", "f": "F"}[branch.typecode])

if args.decay:
  lhefileclass = lhefile.LHEFile_Hwithdecay
else:
  lhefileclass = lhefile.LHEFile_JHUGenVBFVH

with lhefileclass(args.infile, isgen=True) as f:
  for i, event in enumerate(f, start=1):
    if args.decay:
      event.computeP()
    else:
      event.computeProdP()
    if args.vbf:
      Q2V1, Q2V2, branches["GenhelcosthetaV1_VBF"][0], branches["GenhelcosthetaV2_VBF"][0], branches["Genhelphi_VBF"][0], branches["Gencosthetastar_VBF"][0], branches["GenphistarV1_VBF"][0] = event.computeVBFAngles()
      branches["GenQ_V1"][0] = Q2V1 ** 0.5
      branches["GenQ_V2"][0] = Q2V2 ** 0.5
      branches["GenDiJetMass"][0] = (event.associated[0].second + event.associated[1].second).M()
      branches["GenDRjet"][0] = event.associated[0].second.DeltaR(event.associated[1].second)
    if args.zh or args.wh:
      if args.zh:
        if any(11 <= a.first <= 16 for a in event.associated): prod = TVar.Lep_ZH
        else: prod = TVar.Had_ZH
      else:
        if any(11 <= a.first <= 16 for a in event.associated): prod = TVar.Lep_WH
        else: prod = TVar.Had_WH
      branches["GenhelcosthetaV1_VH"][0], branches["GenhelcosthetaV2_VH"][0], branches["Genhelphi_VH"][0], branches["Gencosthetastar_VH"][0], branches["phistarV1_VH"] = event.computeVHAngles(prod)

    if args.decay:
      branches["GenHMass"][0], branches["GenZ1Mass"][0], branches["GenZ2Mass"][0], branches["GenhelcosthetaZ1"][0], branches["GenhelcosthetaZ2"][0], branches["Genhelphi"][0], branches["Gencosthetastar"][0], branches["GenphistarZ1"][0] = event.computeDecayAngles()

      daughters = event.daughters

      for daughter1, daughter2 in itertools.combinations(daughters, 2):
        daughter1, daughter2 = sorted((daughter1, daughter2), key=lambda x: -x.first)
        otherdaughters = sorted((d for d in daughters if d != daughter1 and d != daughter2), key=lambda x: -x.first)
        if len(otherdaughters) != 2:
          raise ValueError("Wrong number of daughters:\n" + "\n".join("  {} {} {} {} {}".format(d.first, d.second.Px(), d.second.Py(), d.second.Pz(), d.second.E()) for d in [daughter1, daughter2]+otherdaughters))
        daughter3, daughter4 = otherdaughters
        if (
          abs((daughter1.second + daughter2.second).M() - branches["GenZ1Mass"][0]) < 1e-5
        ) and (
          abs((daughter3.second + daughter4.second).M() - branches["GenZ2Mass"][0]) < 1e-5
        ):
          break
      else:
        raise ValueError(
            "Couldn't match daughters to masses.  Daughters:\n"
          + "\n".join("  {} {} {} {} {}".format(d.first, d.second.Px(), d.second.Py(), d.second.Pz(), d.second.E()) for d in daughters)
          + "Masses:\n  {}\n  {}".format(branches["GenZ1Mass"][0], branches["GenZ2Mass"][0])
        )

      branches["GenLep1Id"][0] = daughter1.first
      branches["GenLep2Id"][0] = daughter2.first
      branches["GenLep3Id"][0] = daughter3.first
      branches["GenLep4Id"][0] = daughter4.first

    for name, branch in branches.iteritems():
      if not isinstance(branch, array.array):
        raise TypeError(name+' is no longer an array, maybe you assigned to branches["'+name+'"] instead of branches["'+name+'"][0]?')
      if branch[0] == -99999:
        raise ValueError('Forgot to assign to branches["'+name+'"]')
    t.Fill()

    if i % 1000 == 0:
      print "processed", i, "events"

newf.Write()
newf.Close()
