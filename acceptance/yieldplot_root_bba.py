import uproot
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.stats
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from overloads import ROOT, SinglePadCanvas, niceTLegend
from helpers import ATLASInternal, EnergyAndLumi, DrawText, ATLASPre
import re
legend_title = {
    "r1": "1 #font[52]{b}-tag, resolved",
    "r2": "2 #font[52]{b}-tags, resolved",
    "r3": "3+ #font[52]{b}-tags, resolved",
    "m1": "1 #font[52]{b}-tag, 0 add., merged",
    "m2": "2 #font[52]{b}-tags, 0 add., merged",
    "m1a": "1 #font[52]{b}-tag, 1+ add., merged",
    "m2a": "2 #font[52]{b}-tags, 1+ add., merged",
    "m12a": "1+2 #font[52]{b}-tags, 1+ add., merged",
    "sum": "All signal regions",
}
graph_style = {
    "r1": {"color": ROOT.kTeal - 8, "markerstyle": 8, "markersize": 0.75, "linestyle": 1, "linewidth": 2},
    "r2": {"color": ROOT.kTeal - 8, "markerstyle": 8, "markersize": 0.75, "linestyle": 2, "linewidth": 2},
    "r3": {"color": ROOT.kTeal - 8, "markerstyle": 8, "markersize": 0.75, "linestyle": 3, "linewidth": 2},
    "m1": {
        "color": ROOT.kAzure - 8,
        "markerstyle": 8,
        "markersize": 0.75,
        "linestyle": 1,
        "linewidth": 2,
    },
    "m2": {
        "color": ROOT.kAzure - 8,
        "markerstyle": 8,
        "markersize": 0.75,
        "linestyle": 2,
        "linewidth": 2,
    },
    "m12a": {
        "color": ROOT.kAzure - 8,
        "markerstyle": 8,
        "markersize": 0.75,
        "linestyle": 3,
        "linewidth": 2,
    },
    "m1a": {
        "color": ROOT.kAzure - 8,
        "markerstyle": 8,
        "markersize": 0.75,
        "linestyle": 3,
        "linewidth": 2,
    },
    "m2a": {
        "color": ROOT.kAzure - 8,
        "markerstyle": 8,
        "markersize": 0.75,
        "linestyle": 4,
        "linewidth": 2,
    },
    "sum": {"color": ROOT.kBlack, "markerstyle": 8, "markersize": 1.2, "linestyle": 1, "linewidth": 3},
}
text = {
    "HVTZHllqq": {"title": "2L channel, HVT Z' #rightarrow Zh #rightarrow llbb, llcc", "xtitle": r"#font[52]{m_{Z'}} [GeV]"},
    "AZhllbb": {
        "title": "2L channel, gg #rightarrow A #rightarrow Zh #rightarrow llbb",
        "xtitle": r"#font[52]{A} resonance mass [GeV]",
    },
    "bbAZhllbb": {
        "title": "2L channel, gg #rightarrow b#bar{b}A, A #rightarrow Zh",
        "xtitle": r"#font[52]{A} resonance mass [GeV]",
    },
}
def apply_style(graph, region):
    graph.SetLineStyle(graph_style[region]["linestyle"])
    graph.SetLineWidth(graph_style[region]["linewidth"])
    graph.SetLineColor(graph_style[region]["color"])

    graph.SetMarkerStyle(graph_style[region]["markerstyle"])
    graph.SetMarkerSize(graph_style[region]["markersize"])
    graph.SetMarkerColor(graph_style[region]["color"])
    
def main():
    # samplenameinhist = "AZhllbb"
    # samplenameininfo = "ggA"
    # signalname = "AZhllbb"

    # samplenameinhist = "HVTZHllqq"
    # samplenameininfo = "HVT"
    # signalname = "HVTZHllqq"
    
    samplenameinhist = "bbAZhllbb"
    samplenameininfo = "bbA"
    signalname = "bbAZhllbb"
    ids = []
    masses = []
    yields = {}
    selected = {}
    selectedr1 = {}
    selectedr2 = {}
    selectedr3 = {}
    selectedm1 = {}
    selectedm2 = {}
    selectedm1a = {}
    selectedm2a = {}
    selectedm12a = {}

    with open("sample_info.txt") as f:
        for eachline in f:
            sample_tem = eachline.split(" ")
            sample = []
            for each in sample_tem:
                if each != "" and each != "\n":
                    sample.append(each)
            if samplenameininfo in sample[1]:
                ids.append(int(sample[0]))
            else:
                continue
            if samplenameininfo == "ggA" or samplenameininfo == "bbA":
                for each in sample[2].split("_"):
                    if samplenameininfo in each:
                        masses.append(int(each.replace(samplenameininfo, "")))
                        break
            elif samplenameininfo == "HVT":
                masses.append(int(re.findall("\d+", sample[2].split("_")[-1])[0]))

    
    yieldfiles = ["a.txt", "d.txt", "e.txt"]
    for eachfile in yieldfiles:
        with open(eachfile) as f:
            for eachline in f:
                sample_tem = eachline.split(" ")
                sample = []
                for each in sample_tem:
                    if each != "" and each != "\n":
                        sample.append(each)
                if int(sample[0]) in ids:
                    i = ids.index(int(sample[0]))
                    if masses[i] not in yields:
                        yields[masses[i]] = int(sample[1])
                    else:
                        yields[masses[i]] += int(sample[1])
    
    f = uproot.open("run2dbl_bbA.root")
    # for each_histname in f.keys():
    #     if "bbA" in each_histname.decode("utf-8") and "300_" in each_histname.decode("utf-8") and "_SR_" in each_histname.decode("utf-8"):
    #         print(each_histname)
    # exit(1)
    print(masses)
    for each_histname in f.keys():
        each_histname_b = each_histname
        each_histname = each_histname.decode("utf-8")
        for each_mass in masses:
            if samplenameinhist  + str(each_mass) + "_" in each_histname and "_SR_" in each_histname:
                # if each_mass == 500:
                #     print(each_histname)
                if "0tag" in each_histname:
                    continue
                if "bbA" not in each_histname:
                    continue
                if samplenameininfo != "bbA":
                    if "topaddbjetcr" in each_histname or "4ptag2pjet" in each_histname or "3tag2pjet" in each_histname or "3ptag2pjet" in each_histname:
                        continue
                if each_mass not in selected:
                    selected[each_mass] = f[each_histname_b]._fEntries
                else:
                    selected[each_mass] += f[each_histname_b]._fEntries
                if "1tag2pjet" in each_histname:
                    if each_mass not in selectedr1:
                        selectedr1[each_mass] = f[each_histname_b]._fEntries
                    else:
                        print("error1")
                        exit(1)
                if "2tag2pjet" in each_histname:
                    if each_mass not in selectedr2:
                        selectedr2[each_mass] = f[each_histname_b]._fEntries
                    else:
                        print("error2")
                        exit(1)
                if "3ptag2pjet" in each_histname:
                    if each_mass not in selectedr3:
                        selectedr3[each_mass] = f[each_histname_b]._fEntries
                    else:
                        print("error3")
                        exit(1)
                if "1tag1pfat0pjet_0ptv_SR_noaddbjetsr" in each_histname:
                    if each_mass not in selectedm1:
                        selectedm1[each_mass] = f[each_histname_b]._fEntries
                    else:
                        print("error4")
                        exit(1)
                if "2tag1pfat0pjet_0ptv_SR_noaddbjetsr_mVH" in each_histname:
                    if each_mass not in selectedm2:
                        selectedm2[each_mass] = f[each_histname_b]._fEntries
                    else:
                        print("error5")
                        exit(1)
                if "1tag1pfat0pjet_0ptv_SR_topaddbjetcr" in each_histname:
                    if each_mass not in selectedm1a:
                        selectedm1a[each_mass] = f[each_histname_b]._fEntries
                    else:
                        print("error6")
                        exit(1)
                if "2tag1pfat0pjet_0ptv_SR_topaddbjetcr_mVH" in each_histname:
                    if each_mass not in selectedm2a:
                        selectedm2a[each_mass] = f[each_histname_b]._fEntries
                    else:
                        print("error4")
                        exit(1)
                if "2tag1pfat0pjet_0ptv_SR_topaddbjetcr_mVH" in each_histname or "1tag1pfat0pjet_0ptv_SR_topaddbjetcr" in each_histname:
                    if each_mass not in selectedm12a:
                        selectedm12a[each_mass] = f[each_histname_b]._fEntries
                    else:
                        selectedm12a[each_mass] += f[each_histname_b]._fEntries
    print(selected)
    masses_tem = sorted(masses)
    masses = []
    for each in masses_tem:
        if each < 5010:
            masses.append(each)

    acceptance = []
    for each in masses:
        if each in selected:
            acceptance.append(selected[each]/yields[each])
        print(each, selected[each], yields[each], selected[each]/yields[each])

    acceptancer1 = []
    for each in masses:
        if each in selectedr1:
            acceptancer1.append(selectedr1[each]/yields[each])
        else:
            acceptancer1.append(0)
    
    acceptancer2 = []
    for each in masses:
        if each in selectedr2:
            acceptancer2.append(selectedr2[each]/yields[each])
        else:
            acceptancer2.append(0)

    acceptancer3 = []
    for each in masses:
        if each in selectedr3:
            acceptancer3.append(selectedr3[each]/yields[each])
        else:
            acceptancer3.append(0)

    acceptancem1 = []
    for each in masses:
        if each in selectedm1:
            acceptancem1.append(selectedm1[each]/yields[each])
        else:
            acceptancem1.append(0)

    acceptancem2 = []
    for each in masses:
        if each in selectedm2:
            acceptancem2.append(selectedm2[each]/yields[each])
        else:
            acceptancem2.append(0)


    acceptancem1a = []
    for each in masses:
        if each in selectedm1a:
            acceptancem1a.append(selectedm1a[each]/yields[each])
        else:
            acceptancem1a.append(0)

    acceptancem2a = []
    for each in masses:
        if each in selectedm2a:
            acceptancem2a.append(selectedm2a[each]/yields[each])
        else:
            acceptancem2a.append(0)

    acceptancem12a = []
    for each in masses:
        if each in selectedm12a:
            acceptancem12a.append(selectedm12a[each]/yields[each])
        else:
            acceptancem12a.append(0)

    efficiencies = {}
    efficiencies["r1"] = acceptancer1
    efficiencies["r2"] = acceptancer2
    efficiencies["r3"] = acceptancer3
    efficiencies["m1"] = acceptancem1
    efficiencies["m2"] = acceptancem2
    efficiencies["m1a"] = acceptancem1a
    efficiencies["m2a"] = acceptancem2a
    efficiencies["m12a"] = acceptancem12a
    efficiencies["sum"] = acceptance
    
    canv = SinglePadCanvas()
    canv.cd()
    xax, yax = None, None
    legend = niceTLegend()
    __root_objects = []
    __root_objects.append(legend)
    sumlist = ["sum"] + ["r1", "r2", "m1", "m2"]
    if samplenameininfo == "bbA":
        sumlist = ["sum"] + ["r1", "r2", "r3", "m1", "m2", "m12a"]
    for region in sumlist:
        graph = ROOT.TGraph()
        __root_objects.append(graph)
        graph.SetTitle(region)
        for idx, (mass, eff) in enumerate(zip(masses, efficiencies[region])):
            graph.SetPoint(idx, mass, eff)
        if xax is None:
            graph.Draw("APL")
            xax = graph.GetXaxis()
            yax = graph.GetYaxis()
            graph.SetMinimum(0)
            graph.SetMaximum(1.8 * np.max(efficiencies["sum"]))
        else:
            graph.Draw("PLsame")
        apply_style(graph, region)
        legend.AddEntry(graph, legend_title[region], "PL")

    canv.paper_style()
    canv.xtitle = text[signalname]["xtitle"]
    canv.ytitle = "Acceptance #times efficiency"
    canv.title = ""
    canv.xrange = (100, 5100)
    if "HVT" not in signalname:
        canv.xrange = (250, 2100)

    yax.SetNdivisions(505)

    legend.UpdateCoords(0.52, 0.65, 0.9, 0.93)
    legend.Draw()

    ATLASInternal(x=0.15, y=0.88)
    #ATLASPre(x=0.15, y=0.88)
    EnergyAndLumi(x=0.15, y=0.82, lumi=139, size=0.035)
    DrawText(text[signalname]["title"], x=0.15, y=0.78, size=0.035)

    canv.save(f"signal_efficiency_{signalname}")

if __name__ == "__main__":
    main()
