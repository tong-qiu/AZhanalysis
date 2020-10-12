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

def poplist(inputs, index):
    out = []
    for i in range(len(inputs)):
        if i not in index:
            out.append(inputs[i])
    return out

def formatnumber(value, error):
    valuestring = "{:.2e}".format(value)
    valuedigit = int(valuestring.split("e")[1])
    errorstring = "{:.2e}".format(error)
    errordigit = int(errorstring.split("e")[1])
    error = float(errorstring.split("e")[0]) * 10.**(-valuedigit + errordigit)

    # error1 = str(error).split(".")[0]
    # error2 = str(error).split(".")[1]
    # if len(error2) > 15:
    #     error2 = error2[0:14]
    #     error2 = str(int(error2)).rstrip("0")
    # error = float(error1 + "." + error2)
    error = round(error,2)
    return("(" + valuestring.split("e")[0] + "+/-" + str(error) + ")e" + valuestring.split("e")[1])

def texformatnumber(value, error):
    origin = formatnumber(value, error).replace("+/-", r"\pm")
    p1 = origin.split("e")[0]
    p2 = origin.split("e")[1]
    return "$" + p1 + r" \times 10^{" + p2 + "}$"

# 0 xp, 1 sigma, 2 ksai, 3 p1, 4 p2
def fitfunction_root(x, par):
    # source: https://github.com/root-project/root/blob/master/roofit/roofit/src/RooBukinPdf.cxx
    Xp = par[0] # x_p
    sigp = par[1] # \sigma_p
    xi = par[2] # \xi
    rho1 = par[3] # \rho_1
    rho2 = par[4] # \rho_2
    ap = par[5] # \ \A_p
    consts = 2*math.sqrt(2*math.log(2.0))
    r1=0
    r2=0
    r3=0
    r4=0
    r5=0
    hp=0
    x1 = 0
    x2 = 0
    fit_result = 0

    hp=sigp*consts
    r3=ROOT.TMath.Log(2.)
    r4=math.sqrt(xi*xi+1)
    r1=xi/r4

    if abs(xi) > math.exp(-6.):
        r5=xi/ROOT.TMath.Log(r4+xi)
    else:
        r5=1

    x1 = Xp + (hp / 2) * (r1-1)
    x2 = Xp + (hp / 2) * (r1+1)

    #--- Left Side
    if x[0] < x1:
        r2=rho1*(x[0] - x1) * (x[0] - x1)/(Xp - x1)/(Xp - x1)-r3 + 4 * r3 * (x[0] - x1)/hp * r5 * r4/(r4 - xi)/(r4 - xi)
    #--- Center
    elif x[0] < x2:
        if abs(xi) > math.exp(-6.):
            r2=ROOT.TMath.Log(1 + 4 * xi * r4 * (x[0] - Xp)/hp)/ROOT.TMath.Log(1 + 2*xi*(xi - r4))
            r2=-r3*r2*r2
        else:
            r2=-4*r3*(x[0] - Xp)*(x[0] - Xp)/hp/hp
    #--- Right Side
    else:
        r2=rho2*(x[0] - x2)*(x[0] - x2)/(Xp - x2)/(Xp - x2)-r3 - 4 * r3 * (x[0] - x2)/hp * r5 * r4/(r4 + xi)/(r4 + xi)

    if abs(r2) > 100:
        fit_result = 0
    else:
        #---- Normalize the result
        fit_result = math.exp(r2)
    return fit_result *ap


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
    doplot = True
    samplenameinhist = "AZhllbb"
    samplenameininfo = "ggA"
    signalSymbol = "A"
    # samplenameinhist = "HVTZHllqq"
    # samplenameininfo = "HVT"
    # signalSymbol = "Z'"
    # samplenameinhist = "bbAZhllbb"
    # samplenameininfo = "bbA"
    # signalSymbol = "A"
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
            elif samplenameininfo == "bbA":
                for each in sample[2].split("_"):
                        if samplenameininfo in each:
                            masses.append(int(each.replace(samplenameininfo, "")))
                            break
    masses = sorted(masses)
    masses_tem = masses
    masses = []
    for each in masses_tem:
        if each < 5001:
            masses.append(each)
    
    f = uproot.open("dblnominal.root")
    for each_histname in f.keys():
        each_histname_b = each_histname
        each_histname = each_histname.decode("utf-8")
        for each_mass in masses:
            if samplenameinhist  + str(each_mass) + "_" in each_histname and "_SR_" in each_histname and "_mVHresolution" in each_histname:
                if samplenameininfo != "bbA":
                    if "bbA" in each_histname:
                        continue
                    if "topaddbjetcr" in each_histname or "4ptag2pjet" in each_histname or "3tag2pjet" in each_histname or "3ptag2pjet" in each_histname or "0tag" in each_histname:
                        continue
                else:
                    if "4ptag2pjet" in each_histname or "3tag2pjet" in each_histname or "0tag" in each_histname:
                        continue
                bbAtag = samplenameininfo == "bbA" and "3ptag2pjet" in each_histname
                if "1tag2pjet" in each_histname or "2tag2pjet" in each_histname or bbAtag:
                    if each_mass not in selectedr1:
                        selectedr1[each_mass] = [(np.array(f[each_histname_b].edges[0:-1]) + np.array(f[each_histname_b].edges[1:]))/2, np.array(f[each_histname_b].values), np.array(f[each_histname_b].variances)**0.5]
                    else:
                        selectedr1[each_mass][1] += np.array(f[each_histname_b].values)
                        selectedr1[each_mass][2] = (selectedr1[each_mass][2]**2 + np.array(f[each_histname_b].variances))**0.5
                bbAtag = samplenameininfo == "bbA" and "SR_topaddbjetcr" in each_histname
                if "1tag1pfat0pjet_0ptv_SR_noaddbjetsr" in each_histname or "2tag1pfat0pjet_0ptv_SR_noaddbjetsr_mVH" in each_histname or bbAtag:
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
    #ROOT.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(90000000)
    #ROOT.Math.MinimizerOptions.SetDefaultTolerance(10000); 
    for each_mass in masses:
        if each_mass > 3000:
            continue
        masses_resolved.append(each_mass)
        graph = ROOT.TGraphErrors()
        for i in range(len(selectedr1[each_mass][0])):
            graph.SetPoint(i, selectedr1[each_mass][0][i], selectedr1[each_mass][1][i])
            graph.SetPointError(i, 0, selectedr1[each_mass][2][i])
        fit1 = ROOT.TF1("fit1", fitfunction_root,-1, 1, 6)
        fit1.SetParameters(0,0.1,0,0,0,100)
        print("Fit for resolved:", each_mass)
        Ntry = 0
        while 1:
            fitStatus = graph.Fit(fit1)
            result = fit1.GetParameters()
            error = fit1.GetParErrors()
            if error[0] > 0.3:
                print("not converge, try again...")
                fit1.SetParameters(result[0], result[1], result[2], result[3], result[4], result[5])
                Ntry += 1
                if Ntry > 4:
                    print("fit failed")
                    break
                continue
            break
        
        effs.append(result[1])
        errors.append(error[1])

        if doplot:

            plotx = np.linspace(-1, 1, 1000)
            ploty = []
            for each in plotx:
                ploty.append(fitfunction_root([each], [result[0], result[1], result[2], result[3], result[4], result[5]]))
            plt.plot(plotx, ploty, zorder=2)
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
            plt.text(0.05, 0.75, "$m_{" + signalSymbol + "}$ = " + str(each_mass) + " GeV, 2 lep. resolved", fontsize=11, transform=ax.transAxes)

            plt.text(0.57, 0.95, r"$x_p$ = " + texformatnumber(result[0], error[0]), fontsize=10, transform=ax.transAxes, style='italic')
            plt.text(0.57, 0.90, r"$\sigma_p$ = " + texformatnumber(result[1], error[1]), fontsize=10, transform=ax.transAxes, style='italic')
            plt.text(0.57, 0.85, r"$\xi$ = " + texformatnumber(result[2], error[2]), fontsize=10, transform=ax.transAxes, style='italic')
            plt.text(0.57, 0.80, r"$\rho_1$ = " + texformatnumber(result[3], error[3]), fontsize=10, transform=ax.transAxes, style='italic')
            plt.text(0.57, 0.75, r"$\rho_2$ = " + texformatnumber(result[4], error[4]), fontsize=10, transform=ax.transAxes, style='italic')
            plt.text(0.57, 0.70, r"$A_p$ = " + texformatnumber(result[5], error[5]), fontsize=10, transform=ax.transAxes, style='italic')

            plt.savefig("fitplot/" + str(each_mass) + "resolved"+".pdf" ,bbox_inches='tight', pad_inches = 0.02)
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
        fit1 = ROOT.TF1("fit1", fitfunction_root,-1, 1, 6)
        fit1.SetParameters(0,0.1,0,0,0,100)
        print("Fit for merged:", each_mass)
        Ntry = 0
        while 1:
            graph.Fit(fit1)
            result = fit1.GetParameters()
            error = fit1.GetParErrors()
            if error[0] > 0.3:
                print("not converge, try again...")
                fit1.SetParameters(result[0], result[1], result[2], result[3], result[4], result[5])
                Ntry += 1
                if Ntry > 4:
                    print("fit failed")
                    break
                continue
            break
        effsm.append(result[1])
        errorsm.append(error[1])
        if doplot:
            plotx = np.linspace(-1, 1, 1000)
            ploty = []
            for each in plotx:
                ploty.append(fitfunction_root([each], [result[0], result[1], result[2], result[3], result[4], result[5]]))
            plt.plot(plotx, ploty, zorder=2)
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
            plt.text(0.05, 0.75, "$m_{" + signalSymbol + "}$ = " + str(each_mass) + " GeV, 2 lep. merged", fontsize=11, transform=ax.transAxes)

            plt.text(0.57, 0.95, r"$x_p$ = " + texformatnumber(result[0], error[0]), fontsize=10, transform=ax.transAxes, style='italic')
            plt.text(0.57, 0.90, r"$\sigma_p$ = " + texformatnumber(result[1], error[1]), fontsize=10, transform=ax.transAxes, style='italic')
            plt.text(0.57, 0.85, r"$\xi$ = " + texformatnumber(result[2], error[2]), fontsize=10, transform=ax.transAxes, style='italic')
            plt.text(0.57, 0.80, r"$\rho_1$ = " + texformatnumber(result[3], error[3]), fontsize=10, transform=ax.transAxes, style='italic')
            plt.text(0.57, 0.75, r"$\rho_2$ = " + texformatnumber(result[4], error[4]), fontsize=10, transform=ax.transAxes, style='italic')
            plt.text(0.57, 0.70, r"$A_p$ = " + texformatnumber(result[5], error[5]), fontsize=10, transform=ax.transAxes, style='italic')
            plt.savefig("fitplot/" + str(each_mass) + "merged"+".pdf" ,bbox_inches='tight', pad_inches = 0.02)
            plt.clf()
            plt.cla()
            plt.close()

    plt.errorbar(masses_resolved, effs, yerr=errors, fmt='b.-', label = "resolved")
    notconverge = []
    for i in range(len(masses_merged)):
        if errorsm[i] > 0.1:
            notconverge.append(i)
    masses_merged = poplist(masses_merged, notconverge)
    effsm = poplist(effsm, notconverge)
    errorsm = poplist(errorsm, notconverge)
    plt.errorbar(masses_merged, effsm, yerr=errorsm, fmt='r.-', label = "merged")
    plt.legend(frameon=False, prop={'size': 13})
    plt.ylim([0, max(effs + effsm)*1.5 ])
    title1 = r"ATLAS"
    title1_1 = r"Internal"
    title3 = r"$\sqrt{s}=13\:TeV,139\:fb^{-1}$"
    ax = plt.gca()
    plt.ylabel(r"resolution", fontsize=17)
    plt.xlabel("$m_{" + signalSymbol + "}$ [GeV]", fontsize=17)
    plt.text(0.05, 0.9, title1, fontsize=18, transform=ax.transAxes, style='italic', fontweight='bold')
    plt.text(0.25, 0.9, title1_1, fontsize=18, transform=ax.transAxes)
    plt.text(0.05, 0.82, title3, fontsize=14, weight='bold', style='italic', transform=ax.transAxes)
    plt.text(0.05, 0.75, samplenameininfo + ", 2 lep", fontsize=14, transform=ax.transAxes)
    plt.savefig("fitplot/overall.pdf" ,bbox_inches='tight', pad_inches = 0.02)
    plt.clf()
    plt.cla()
    plt.close()


if __name__ == "__main__":
    main()