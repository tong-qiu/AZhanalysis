import ROOT
import numpy as np
import os
import json
import sys
lib_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.append(lib_path)
from package.events import *
import package
from package.curveplot import curveplot

def ptll_cut(mvh):
    if mvh < 320.:
        return 0.
    return 20. + 9. * pow(mvh - 320., 0.5)

def resolvedselection(entry):
    if entry.nsigjet < 2:
        return False
    if entry.nbjets == 0:
        return False
    ptl1 = entry.ptl1
    ptl2 = entry.ptl2
    if ptl2 > ptl1:
        ptl1, ptl2 = ptl2, ptl1
    if ptl2/1000. < 20:
        return False
    if ptl1/1000. < 27:
        return False
    if entry.ptb1/1000. < 45:
        return False
    lowmll = max(40, 87 - 0.03 * entry.mVHres/1000.)
    highmll = 97 + 0.013 * entry.mVHres/1000.
    if entry.mll/1000. < lowmll:
        return False
    if entry.mll/1000. > highmll:
        return False
    if entry.mbbres/1000. < 100:
        return False
    if entry.mbbres/1000. > 145:
        return False
    if entry.ptll/1000. < ptll_cut(entry.mVHres/1000.):
        return False
    # print(entry.met/1000. / (entry.ht/1000.)**0.5, 1.15 + 8 * (10**-3) * entry.mVHres/1000.)
    # if entry.met/1000. / (entry.ht/1000.)**0.5 > 1.15 + 8 * (10**-3) * entry.mVHres/1000.:
    #     return False
    # print("5")
    return True

def mergedselection(entry):
    if entry.ptfj1/1000. < 250:
        return False

    ptl1 = entry.ptl1
    ptl2 = entry.ptl2
    if ptl2/1000. < 25:
        return False
    if ptl1/1000. < 27:
        return False
    if entry.mfj1/1000. < 75:
        return False
    if abs(entry.etafj1) > 2:
        return False
    if entry.mfj1/1000. > 145:
        return False

    lowmll = max(40, 87 - 0.03 * entry.mVHmerg/1000.)
    highmll = 97 + 0.013 * entry.mVHmerg/1000.
    if entry.mll/1000. < lowmll:
        return False
    if entry.mll/1000. > highmll:
        return False
    if entry.ptll/1000. < ptll_cut(entry.mVHmerg/1000.):
        return False
    return True

def loadnumber(path):
    resolved = 0
    merged = 0
    f = ROOT.TFile(path)
    t1 = f.Get("data")
    dumped = 0
    for entry in t1:
        if entry.nlepton != 2:
            continue
        if resolvedselection(entry):
            resolved += 1
        elif mergedselection(entry):
            merged += 1
        else:
            dumped += 1
    print(resolved, merged, dumped)
    return (resolved, merged)


def main():
    loadfile = False
    acceptance = {}
    if loadfile:
        massset = set()
        for each in os.listdir("ntuples"):
            if ".root" in each:
                massset.add(int(each.split("w")[0]))
        for each_mass in sorted(massset):
            print("loading" + str(each_mass))
            total_resolved, total_merged = loadnumber(
                "ntuples/" + str(each_mass) + "w" + str(0) + ".root")
            acceptance[each_mass] = {}
            for each_width in [1, 2, 5, 10, 20]:
                path = "ntuples/" + str(each_mass) + \
                    "w" + str(each_width) + ".root"
                each_resolved, each_merged = loadnumber(path)
                acceptance[each_mass][each_width] = [
                    each_resolved/total_resolved, each_merged/total_merged]
        with open('acceptance.json', 'w') as f:
            json.dump(acceptance, f)
        exit(1)
    else:
        with open('acceptance.json', 'r') as f:
            for eachline in f:
                acceptance = json.loads(eachline)
                break
    masses = []
    w1 = []
    w2 = []
    w5 = []
    w10 = []
    w20 = []
    fakeerror = []
    for each_mass in acceptance.keys():
        masses.append(int(each_mass))
        w1.append(acceptance[each_mass]['1'][0])
        w2.append(acceptance[each_mass]['2'][0])
        w5.append(acceptance[each_mass]['5'][0])
        w10.append(acceptance[each_mass]['10'][0])
        w20.append(acceptance[each_mass]['20'][0])
        fakeerror.append(0)
    xlist = [w1, w2, w5, w10, w20]
    label = ["1%", "2%", "5%", "10%", "20%"]
    curveplot([masses] * len((xlist)), xlist, labels=label, ylabel="acceptance ratio", filename="acceptance_resolved", title3="resolved",
    ylimit=[0, 3] ,yshift=0.05)


    masses = []
    w1 = []
    w2 = []
    w5 = []
    w10 = []
    w20 = []
    fakeerror = []
    for each_mass in acceptance.keys():
        masses.append(int(each_mass))
        w1.append(acceptance[each_mass]['1'][1])
        w2.append(acceptance[each_mass]['2'][1])
        w5.append(acceptance[each_mass]['5'][1])
        w10.append(acceptance[each_mass]['10'][1])
        w20.append(acceptance[each_mass]['20'][1])
        fakeerror.append(0)
    xlist = [w1, w2, w5, w10, w20]
    label = ["1%", "2%", "5%", "10%", "20%"]
    curveplot([masses] * len((xlist)), xlist, labels=label, ylabel="acceptance ratio", filename="acceptance_merged", title3="merged",
    ylimit=[0, 3] ,yshift=0.05)

if __name__ == "__main__":
    main()
