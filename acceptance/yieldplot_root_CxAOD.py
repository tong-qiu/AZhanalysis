import uproot
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.stats
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from overloads import ROOT, SinglePadCanvas, niceTLegend
from helpers import ATLASInternal, EnergyAndLumi, DrawText
import re
legend_title = {
    "r1": "1 #font[52]{b}-tag, resolved",
    "r2": "2 #font[52]{b}-tags, resolved",
    "m1": "1 #font[52]{b}-tag, merged",
    "m2": "2 #font[52]{b}-tags, merged",
    "sum": "All signal regions",
}
graph_style = {
    "r1": {"color": ROOT.kTeal - 8, "markerstyle": 8, "markersize": 0.75, "linestyle": 1, "linewidth": 2},
    "r2": {"color": ROOT.kTeal - 8, "markerstyle": 8, "markersize": 0.75, "linestyle": 2, "linewidth": 2},
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
    "sum": {"color": ROOT.kBlack, "markerstyle": 8, "markersize": 1.2, "linestyle": 1, "linewidth": 3},
}
text = {
    "HVTZHllqq": {"title": "2L channel, HVT Z'#rightarrow Zh", "xtitle": "#font[52]{Z'} resonance mass [GeV]"},
    "AZhllbb": {
        "title": "2L channel, gg #rightarrow A #rightarrow Zh",
        "xtitle": "#font[52]{A} resonance mass [GeV]",
    },
    "bbAZhllbb": {
        "title": "2L channel, gg #rightarrow b#bar{b}A, A #rightarrow Zh",
        "xtitle": "#font[52]{A} resonance mass [GeV]",
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
    samplenameinhist = "AZhllbb"
    samplenameininfo = "ggA"
    signalname = "AZhllbb"
    # samplenameinhist = "HVTZHllqq"
    # samplenameininfo = "HVT"
    # signalname = "HVTZHllqq"
    ids = []
    masses = []
    yields = {}
    selected = {}
    selectedr1 = {}
    selectedr2 = {}
    selectedm1 = {}
    selectedm2 = {}

    yields = dict([(300, 2594.636392509799), (1000, 3754.4275185899446), (2000, 4022.7003195593543), (400, 2980.8809280766345), (420, 3034.4605561353465), 
    (440, 3033.2871666791457), (460, 3123.5699001664343), (500, 3211.256911045622), (600, 3399.3625576627287), (700, 3520.623322146133), (800, 3638.826112486981), 
    (900, 3724.7439040558133), (1200, 3871.747200890897), (1400, 3896.416927709315), (1600, 3965.91154826843)])
    masses = sorted(yields.keys())

    
    f = uproot.open("run2dbl.root")
    for each_histname in f.keys():
        each_histname_b = each_histname
        each_histname = each_histname.decode("utf-8")
        for each_mass in masses:
            if samplenameinhist  + str(each_mass) + "_" in each_histname and "_SR_" in each_histname:
                # if each_mass == 500:
                #     print(each_histname)
                if "topaddbjetcr" in each_histname or "4ptag2pjet" in each_histname or "3tag2pjet" in each_histname or "0tag" in each_histname or "bbA" in each_histname:
                    continue
                if each_mass not in selected:
                    selected[each_mass] = sum(f[each_histname_b].values)
                else:
                    selected[each_mass] += sum(f[each_histname_b].values)
                if "1tag2pjet" in each_histname:
                    if each_mass not in selectedr1:
                        selectedr1[each_mass] = sum(f[each_histname_b].values)
                    else:
                        selectedr1[each_mass] += sum(f[each_histname_b].values)
                if "2tag2pjet" in each_histname:
                    if each_mass not in selectedr2:
                        selectedr2[each_mass] = sum(f[each_histname_b].values)
                    else:
                        selectedr2[each_mass] += sum(f[each_histname_b].values)
                if "1tag1pfat0pjet_0ptv" in each_histname:
                    if each_mass not in selectedm1:
                        selectedm1[each_mass] = sum(f[each_histname_b].values)
                    else:
                        selectedm1[each_mass] += sum(f[each_histname_b].values)
                if "2tag1pfat0pjet_0ptv" in each_histname:
                    if each_mass not in selectedm2:
                        selectedm2[each_mass] = sum(f[each_histname_b].values)
                    else:
                        selectedm2[each_mass] += sum(f[each_histname_b].values)
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
        #print(each, selected[each]/yields[each])
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



    efficiencies = {}
    efficiencies["r1"] = acceptancer1
    efficiencies["r2"] = acceptancer2
    efficiencies["m1"] = acceptancem1
    efficiencies["m2"] = acceptancem2
    efficiencies["sum"] = acceptance
    
    canv = SinglePadCanvas()
    canv.cd()
    xax, yax = None, None
    legend = niceTLegend()
    __root_objects = []
    __root_objects.append(legend)
    for region in ["sum"] + ["r1", "r2", "m1", "m2"]:
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

    legend.UpdateCoords(0.6, 0.65, 0.9, 0.93)
    legend.Draw()

    ATLASInternal(x=0.15, y=0.88)
    EnergyAndLumi(x=0.15, y=0.82, lumi=139, size=0.035)
    DrawText(text[signalname]["title"], x=0.15, y=0.78, size=0.035)

    canv.save(f"signal_efficiency_{signalname}_CxAOD")

if __name__ == "__main__":
    main()