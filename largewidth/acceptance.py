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

def loadnumber(path):
    resolved = 0
    merged = 0
    f = ROOT.TFile(path)
    t1 = f.Get("data")
    for entry in t1:

        if entry.havetau == 1:
            continue
        ptl1 = entry.ptl1
        ptl2 = entry.ptl2
        if ptl2 > ptl1:
            ptl1, ptl2 = ptl2, ptl1
        if ptl2/1000. < 20:
            continue
        if ptl1/1000. < 27:
            continue
        if entry.ptb1/1000. < 45:
            continue


        # resolved
        if ptl2/1000. > 25 and entry.pth/1000. > 250:
            merged += 1
        # merged
        if not(ptl2/1000. > 25 and entry.pth/1000. > 250):
            resolved += 1

        # if entry.ptl2/1000. < 20:
        #     continue
        # if entry.ptl1/1000. < 27:
        #     continue
        # if entry.ptb1/1000. < 45:
        #     continue
        # if entry.ptl2/1000. > 25 and entry.ptb1/1000. > 250:
        #     merged += 1
        # if not (entry.ptl2/1000. > 25 and entry.ptb1/1000. > 250):
        #     resolved += 1
    print(resolved, merged)
    return (resolved, merged)


def main():
    loadfile = False
    acceptance = {}
    if loadfile:
        massset = set()
        for each in os.listdir("ntuples"):
            if ".root" in each:
                massset.add(int(each.split("w")[0]))
        massset = [2000]
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
    ylimit=[0, 10] ,yshift=0.05)


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
