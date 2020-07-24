import uproot
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.stats
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import re
def main():
    samplenameinhist = "AZhllbb"
    samplenameininfo = "ggA"
    samplenameinhist = "HVTZHllqq"
    samplenameininfo = "HVT"
    ids = []
    masses = []
    yields = {}
    selected = {}
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
    
    f = uproot.open("run2dbl.root")
    for each_histname in f.keys():
        each_histname_b = each_histname
        each_histname = each_histname.decode("utf-8")
        for each_mass in masses:
            if samplenameinhist  + str(each_mass) + "_" in each_histname and "_SR_" in each_histname:
                if each_mass == 500:
                    print(each_histname)
                if "topaddbjetcr" in each_histname or "4ptag2pjet" in each_histname or "3tag2pjet" in each_histname or "0tag" in each_histname:
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
                if "1tag1pfat0pjet_0ptv_SR_noaddbjetsr" in each_histname:
                    if each_mass not in selectedm1:
                        selectedm1[each_mass] = f[each_histname_b]._fEntries
                    else:
                        print("error3")
                        exit(1)
                if "2tag1pfat0pjet_0ptv_SR_noaddbjetsr_mVH" in each_histname:
                    if each_mass not in selectedm2:
                        selectedm2[each_mass] = f[each_histname_b]._fEntries
                    else:
                        print("error4")
                        exit(1)
    print(selected)
    masses = sorted(masses)
    acceptance = []
    for each in masses:
        if each in selected:
            acceptance.append(selected[each]/yields[each])
        print(each, selected[each]/yields[each])

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


    # zipped_pairs = zip(masses, acceptance) 
    # acceptance = [x for _, x in sorted(zipped_pairs)] 
    # masses = sorted(masses)


    plt.plot(masses, acceptance, ".-k", label="All signal regions")
    plt.plot(masses, acceptancer1, ".-k", label="1 $b$-tag, resolved")
    plt.plot(masses, acceptancer2, ".-k", label="2 $b$-tags, resolved")
    plt.plot(masses, acceptancem1, ".-k", label="1 $b$-tag, merged")
    plt.plot(masses, acceptancem2, ".-k", label="2 $b$-tags, merged")

    plt.xlabel("mass [GeV]", fontsize=17)
    plt.ylabel("efficency x acceptance", fontsize=17)
    plt.ylim([0, 1])
    plt.legend(frameon=False, prop={'size': 15})
    #plt.xlim([masses[0], masses[-1]])
    #plt.yscale("log")
    ax = plt.gca()
    #ax.xaxis.set_minor_locator(MultipleLocator(200))
    title1 = r"ATLAS"
    title1_1 = r"Internal"
    title3 = r"$\sqrt{s}=13\:TeV,139\:fb^{-1}$"
    plt.text(0.05, 0.9, title1, fontsize=18, transform=ax.transAxes, style='italic', fontweight='bold')
    plt.text(0.25, 0.9, title1_1, fontsize=18, transform=ax.transAxes)
    plt.text(0.05, 0.82, title3, fontsize=14, weight='bold', style='italic', transform=ax.transAxes)
    plt.text(0.05, 0.75, samplenameininfo, fontsize=14, weight='bold', style='italic', transform=ax.transAxes)
    plt.savefig(samplenameininfo+".pdf" ,bbox_inches='tight', pad_inches = 0.02)
    plt.show()

if __name__ == "__main__":
    main()