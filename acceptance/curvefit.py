import uproot
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import math
import scipy.stats
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import ROOT
import re
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
matplotlib.use('Agg')
def gaussian(x, norm, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))) * norm

def rebin(inputs, factor, order=1):
    out = []
    i = 0
    sums = 0
    for each in inputs:
        sums += each**order
        i += 1
        if i == factor:
            out.append(sums**(1./order))
            sums = 0
            i = 0
    return out


def main():
    samplenameinhist = "AZhllbb"
    samplenameininfo = "ggA"
    samplenameinhist = "HVTZHllqq"
    samplenameininfo = "HVT"
    ids = []
    masses = []
    yields = {}
    selectedr1 = {}
    selectedr2 = {}
    selectedm1 = {}
    selectedm2 = {}
    
    with open("sample_info.txt") as f:
        for eachline in f:
            sample_tem = eachline.split(" ")
            sample = []
            for each in sample_tem:
                if each != "" and each != "\n":
                    sample.append(each)
            #print(sample[1])
            if samplenameininfo in sample[1]:
                ids.append(int(sample[0]))
            else:
                continue
            if samplenameininfo == "ggA":
                for each in sample[2].split("_"):
                        if samplenameininfo in each:
                            masses.append(int(each.replace(samplenameininfo, "")))
                            break
            elif samplenameininfo == "HVT":
                masses.append(int(re.findall("\d+", sample[2].split("_")[-1])[0]))

    masses = sorted(masses)
    masses_tem = masses
    masses = []
    for each in masses_tem:
        if each < 5001:
            masses.append(each)
    
    f = uproot.open("run2dblnominal.root")
    for each_histname in f.keys():
        each_histname_b = each_histname
        each_histname = each_histname.decode("utf-8")
        for each_mass in masses:
            if samplenameinhist  + str(each_mass) + "_" in each_histname and "_SR_" in each_histname and "_mVHresolution" in each_histname:
                if each_mass == 500:
                    print(each_histname)
                if "topaddbjetcr" in each_histname or "4ptag2pjet" in each_histname or "3tag2pjet" in each_histname or "0tag" in each_histname:
                    continue
                if "1tag2pjet" in each_histname or "2tag2pjet" in each_histname:
                    if each_mass not in selectedr1:
                        selectedr1[each_mass] = [(np.array(f[each_histname_b].edges[0:-1]) + np.array(f[each_histname_b].edges[1:]))/2, np.array(f[each_histname_b].values), np.array(f[each_histname_b].variances)**0.5]
                    else:
                        selectedr1[each_mass][1] += np.array(f[each_histname_b].values)
                        selectedr1[each_mass][2] = (selectedr1[each_mass][2]**2 + np.array(f[each_histname_b].variances))**0.5

                if "1tag1pfat0pjet_0ptv_SR_noaddbjetsr" in each_histname or "2tag1pfat0pjet_0ptv_SR_noaddbjetsr_mVH" in each_histname:
                    if each_mass not in selectedm1:
                        selectedm1[each_mass] = [(np.array(f[each_histname_b].edges[0:-1]) + np.array(f[each_histname_b].edges[1:]))/2, f[each_histname_b].values, np.array(f[each_histname_b].variances)**0.5]
                    else:
                        selectedm1[each_mass][1] += np.array(f[each_histname_b].values)
                        selectedm1[each_mass][2] = (selectedm1[each_mass][2]**2 + np.array(f[each_histname_b].variances))**0.5
    for eachkey in selectedm1.keys():
        selectedm1[eachkey][0] = rebin(selectedm1[eachkey][0], 3)
        selectedm1[eachkey][1] = rebin(selectedm1[eachkey][1], 3)
        selectedm1[eachkey][2] = rebin(selectedm1[eachkey][2], 3, 2)
    for eachkey in selectedr1.keys():
        selectedr1[eachkey][0] = rebin(selectedr1[eachkey][0], 3)
        selectedr1[eachkey][1] = rebin(selectedr1[eachkey][1], 3)
        selectedr1[eachkey][2] = rebin(selectedr1[eachkey][2], 3, 2)
    effs = []
    errors = []
    xlow = -0.5
    xhigh = 0.5
    masses_resolved = []
    for each_mass in masses:
        if each_mass > 3000:
            continue
        masses_resolved.append(each_mass)
        graph = ROOT.TGraphErrors()
        for i in range(len(selectedr1[each_mass][0])):
            graph.SetPoint(i, selectedr1[each_mass][0][i], selectedr1[each_mass][1][i])
            graph.SetPointError(i, 0, selectedr1[each_mass][2][i])
        fit1 = ROOT.TF1("gaussfit","gaus" ,0, 6000)
        graph.Fit(fit1)
        result = fit1.GetParameters()
        error = fit1.GetParErrors()
        effs.append(result[2])
        errors.append(error[2])

        plotx = np.linspace(-1, 1, 1000)
        plt.plot(plotx, gaussian(plotx, result[0], result[1], result[2]), zorder=2)
        plt.errorbar(selectedr1[each_mass][0], selectedr1[each_mass][1], yerr=selectedr1[each_mass][2], fmt='k.', zorder=1)
        plt.ylim([0, max(selectedr1[each_mass][1])*1.5 ])
        # xlow = each_mass - 6 * result[2]
        # if xlow < 0:
        #     xlow = 0
        # xhigh = each_mass + 6 * result[2]
        # if xhigh > 9000:
        #     xhigh = 9000
        plt.xlim([xlow, xhigh])
        title1 = r"ATLAS"
        title1_1 = r"Internal"
        title3 = r"$\sqrt{s}=13\:TeV,139\:fb^{-1}$"
        ax = plt.gca()
        plt.xlabel(r"resolution", fontsize=17)
        plt.ylabel("entries", fontsize=17)
        plt.text(0.05, 0.9, title1, fontsize=18, transform=ax.transAxes, style='italic', fontweight='bold')
        plt.text(0.25, 0.9, title1_1, fontsize=18, transform=ax.transAxes)
        plt.text(0.05, 0.82, title3, fontsize=14, weight='bold', style='italic', transform=ax.transAxes)
        plt.text(0.05, 0.75, "$m_{Z'}$ = " + str(each_mass) + " GeV, 2 lep. resolved", fontsize=14, transform=ax.transAxes)
        plt.savefig("fitplot/hvt" + str(each_mass) + "resolved"+".pdf" ,bbox_inches='tight', pad_inches = 0.02)
        plt.clf()
        plt.cla()
        plt.close()

    effsm = []
    errorsm = []
    xlow = -1
    xhigh = 1
    masses_merged = []
    for each_mass in masses:
        if each_mass < 700:
            continue
        masses_merged.append(each_mass)
        graph = ROOT.TGraphErrors()
        for i in range(len(selectedm1[each_mass][0])):
            graph.SetPoint(i, selectedm1[each_mass][0][i], selectedm1[each_mass][1][i])
            graph.SetPointError(i, 0, selectedm1[each_mass][2][i])
        fit1 = ROOT.TF1("gaussfit","gaus" ,0, 6000)
        graph.Fit(fit1)
        result = fit1.GetParameters()
        error = fit1.GetParErrors()
        effsm.append(result[2])
        errorsm.append(error[2])

        plotx = np.linspace(-1, 1, 1000)
        plt.plot(plotx, gaussian(plotx, result[0], result[1], result[2]), zorder=2)
        plt.errorbar(selectedm1[each_mass][0], selectedm1[each_mass][1], yerr=selectedm1[each_mass][2], fmt='k.', zorder=1)
        plt.ylim([0, max(selectedm1[each_mass][1])*1.5 ])
        # xlow = each_mass - 6 * result[2]
        # if xlow < 0:
        #     xlow = 0
        # xhigh = each_mass + 6 * result[2]
        # if xhigh > 9000:
        #     xhigh = 9000
        plt.xlim([xlow, xhigh])
        title1 = r"ATLAS"
        title1_1 = r"Internal"
        title3 = r"$\sqrt{s}=13\:TeV,139\:fb^{-1}$"
        ax = plt.gca()
        plt.xlabel(r"resolution", fontsize=17)
        plt.ylabel("entries", fontsize=17)
        plt.text(0.05, 0.9, title1, fontsize=18, transform=ax.transAxes, style='italic', fontweight='bold')
        plt.text(0.25, 0.9, title1_1, fontsize=18, transform=ax.transAxes)
        plt.text(0.05, 0.82, title3, fontsize=14, weight='bold', style='italic', transform=ax.transAxes)
        plt.text(0.05, 0.75, "$m_{Z'}$ = " + str(each_mass) + " GeV, 2 lep. merged", fontsize=14, transform=ax.transAxes)
        plt.savefig("fitplot/hvt" + str(each_mass) + "merged"+".pdf" ,bbox_inches='tight', pad_inches = 0.02)
        plt.clf()
        plt.cla()
        plt.close()

    plt.errorbar(masses_resolved, effs, yerr=errors, fmt='b.-', label = "resolved")
    plt.errorbar(masses_merged, effsm, yerr=errorsm, fmt='r.-', label = "merged")
    plt.legend(frameon=False, prop={'size': 13})
    plt.ylim([0, max(effs + effsm)*1.5 ])
    title1 = r"ATLAS"
    title1_1 = r"Internal"
    title3 = r"$\sqrt{s}=13\:TeV,139\:fb^{-1}$"
    ax = plt.gca()
    plt.ylabel(r"resolution", fontsize=17)
    plt.xlabel("$m_{Zh}$ [GeV]", fontsize=17)
    plt.text(0.05, 0.9, title1, fontsize=18, transform=ax.transAxes, style='italic', fontweight='bold')
    plt.text(0.25, 0.9, title1_1, fontsize=18, transform=ax.transAxes)
    plt.text(0.05, 0.82, title3, fontsize=14, weight='bold', style='italic', transform=ax.transAxes)
    plt.text(0.05, 0.75, "HVT Z', 2 lep", fontsize=14, transform=ax.transAxes)
    plt.savefig("fitplot/overall.pdf" ,bbox_inches='tight', pad_inches = 0.02)
    plt.clf()
    plt.cla()
    plt.close()


if __name__ == "__main__":
    main()