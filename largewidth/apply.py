import ROOT
import numpy as np
import os
import pickle
ROOT.gROOT.ProcessLine(".L My_modified_bw.cxx+")
ROOT.gInterpreter.ProcessLine('#include "My_modified_bw.h"')
ROOT.gSystem.Load('My_modified_bw_cxx.so')
from ROOT import My_modified_bw
import json
import multiprocessing
import tqdm

# This is a Python3 script

x = ROOT.RooRealVar("x1","mA" , -2000, 4000.)

def loadhist(inFile):
    out = {}
    for key in inFile.GetListOfKeys():
        tem_hist = key.ReadObj()
        if isinstance(tem_hist, ROOT.TH1F):
            temName = key.GetName()
            # if temName[0:len(filter)] == filter:
            tem_hist.SetDirectory(0)
            out[temName] = tem_hist
    return out

def getupdownhsit(nominal, name):
    hup = nominal.Clone(name + "temup")
    hdown = nominal.Clone(name + "temdown")
    nbins = nominal.GetNbinsX()
    for i in range(nbins):
        staterror = nominal.GetBinError(i)
        value = nominal.GetBinContent(i)
        up = value + staterror
        down = value - staterror
        if down < 0:
            down = 0
        hup.SetBinContent(i, up)
        hdown.SetBinContent(i, down)
    return (hup, hdown)

def runMC(pdf, x, hist, newname, totalweight):
    outhist = hist.Clone(newname)
    outhist.Reset()
    data = pdf.generate(ROOT.RooArgSet(x), 100000)
    for i in range(0, 100000-1):
        outhist.Fill(data.get(i).getRealValue("x1"))
        scale = totalweight/outhist.Integral()
        outhist.Scale(scale)
    return outhist

def smearing(ps, nominalhist, smearingtype, name, newname, SF, outdic):
    hup, hdown = getupdownhsit(nominalhist, newname + "before")
    # x = ROOT.RooRealVar("x1","mA" , -2000, 4000.)
    datahist = ROOT.RooDataHist("datahist" + newname, "datahist", ROOT.RooArgList(x), nominalhist)
    datahistup = ROOT.RooDataHist("datahistup" + newname, "datahistup", ROOT.RooArgList(x), hup)
    datahistdown = ROOT.RooDataHist("datahistdown" + newname, "datahistdown", ROOT.RooArgList(x), hdown)

    sum_nominal = nominalhist.Integral() * SF
    sum_up = hup.Integral() * SF
    sum_down = hdown.Integral() * SF

    pdf = ROOT.RooHistPdf("pdf" + newname, "your pdf", ROOT.RooArgSet(x), datahist)
    pdfup = ROOT.RooHistPdf("pdfup" + newname, "your pdfup", ROOT.RooArgSet(x), datahistup)
    pdfdown = ROOT.RooHistPdf("pdfdown" + newname, "your pdfdown", ROOT.RooArgSet(x), datahistdown)
    smearingpdf = None
    if smearingtype == "BW":
        bwm0 = ROOT.RooRealVar("bwm0" + newname, "bwm0", 0)
        bww = ROOT.RooRealVar("bwwidth" + newname, "bwwidth", ps[2])
        smearingpdf = ROOT.RooBreitWigner("bworiginal" + newname, "bworiginal", x, bwm0, bww)
    else:
        m0 = ROOT.RooRealVar("m0" + newname, "m0", ps[3])
        w = ROOT.RooRealVar("width" + newname, "width", ps[2])
        lnm0 = ROOT.RooRealVar("lnm0" + newname, "lnm0", ps[0])
        lnk = ROOT.RooRealVar("lnmk" + newname, "lnmk", ps[1])
        smearingpdf = My_modified_bw("bwmodified" + newname, "bwmodified", x, w, lnm0, lnk, m0)
    outpdf = ROOT.RooFFTConvPdf("outpdf" + newname, "outpdf", x, pdf, smearingpdf)
    outpdfup = ROOT.RooFFTConvPdf("outpdfup" + newname, "outpdfup", x, pdfup, smearingpdf)
    outpdfdown = ROOT.RooFFTConvPdf("outpdfdown" + newname, "outpdfdown", x, pdfdown, smearingpdf)

    outhist = runMC(outpdf, x, nominalhist, name, sum_nominal)
    outhistup = runMC(outpdfup, x, hup, newname + "up", sum_up)
    outhistdown = runMC(outpdfdown, x, hdown, newname + "down", sum_down)

    nbins = nominalhist.GetNbinsX()
    for i in range(nbins):
        upvalue = outhistup.GetBinContent(i)
        downvalue = outhistdown.GetBinContent(i)
        error = abs(upvalue - downvalue) / 2.
        outhist.SetBinError(i, error)
    outdic[name] = outhist
    return outhist

def multiprocess(hists, pvalues_resolved, pvalues_merged, SFs, width):
    manager = multiprocessing.Manager()
    all_sample = manager.dict()
    processes = []
    missinghist = 0
    with tqdm.tqdm(total=len(hists)) as pbar:
        for i, each_key in enumerate(hists.keys()):
            mass = int(each_key.split("_")[0].split("AZhllbb")[1])
            if mass not in pvalues_merged or mass not in pvalues_resolved or str(mass) not in SFs:
                missinghist += 1
                # print("Error: cannot find parameters for " + each_key)

            new_name = each_key.replace("AZhllbb", "AZhllbb" + str(width))
            if "1pfat0pjet" in new_name:
                sf_tem = SFs[str(mass)][str(width)][1]
                smearingtype = pvalues_merged[mass][width][4]
                p_tem = [pvalues_merged[mass][width][0][0], pvalues_merged[mass][width][1][0], pvalues_merged[mass][width][2][0], pvalues_merged[mass][width][3][0]]
            else:
                sf_tem = SFs[str(mass)][str(width)][0]
                smearingtype = pvalues_resolved[mass][width][4]
                p_tem = [pvalues_resolved[mass][width][0][0], pvalues_resolved[mass][width][1][0], pvalues_resolved[mass][width][2][0], pvalues_resolved[mass][width][3][0]]
            t = multiprocessing.Process(target=smearing, args=(p_tem, hists[each_key], smearingtype, each_key, new_name, sf_tem, all_sample))
            processes.append(t)
            t.start()
            if (i+1)%6 == 0:
                for each in processes:
                    each.join()
                    pbar.update(1)
                    processes = []
        for each in processes:
            each.join()
            pbar.update(1)
    if len(all_sample) != len(hists):
        print("Error. Smearing failed. Following histograms missing:")
        for each in hists.keys():
            if each not in all_sample:
                print(each)
        print(str(missinghist) + " histograms failed due to missing parameters.")
    return all_sample

def getsignalhist(hists, ifbba, debug):
    signalhist = {}
    keys = [*hists.keys()]
    i = 0
    for each in keys:
        if not ifbba:
            if "AZhllbb" == each[0:len("AZhllbb")]:
                signalhist[each] = hists[each]
                hists.pop(each)
                i += 1
        if ifbba:
            if "bbAZhllbb" == each[0:len("bbAZhllbb")]:
                signalhist[each] = hists[each]
                hists.pop(each)
                i += 1
        if debug:
            if i == 14:
                break
    return signalhist

def savehist(nominalHist, sysHist, filename):
    newhistofile = ROOT.TFile(filename + ".root","recreate")
    newhistofile.cd()
    for each_key in nominalHist:
        nominalHist[each_key].Write(each_key)
    cwd = ROOT.gDirectory.mkdir("Systematics")
    cwd.cd()
    for each_key in sysHist:
        sysHist[each_key].Write(each_key)
    newhistofile.Close()


def main():
    debug = False
    bbA = False

    print("Loading parameters...")
    with open("fitvalues_resolved.pickle",'rb') as f:
        pvalues_resolved = pickle.load(f)
    with open("fitvalues_merged.pickle",'rb') as f:
        pvalues_merged = pickle.load(f)
    acceptance = []
    with open('acceptance.json', 'r') as f:
        for eachline in f:
            acceptance = json.loads(eachline)
            break

    print("Loading histograms...")
    nominalFile = ROOT.TFile("run2dbl_nosys.root","read")
    nominalHist = loadhist(nominalFile)
    sysFile = nominalFile.Get("Systematics")
    sysHist = loadhist(sysFile)
    singalhistnominal = getsignalhist(nominalHist, bbA, debug)
    singalhissys = getsignalhist(sysHist, bbA, debug)

    # print(r"Doing 1% smearing..")
    # print(" Smearning nominal...")
    # outdicnominal = multiprocess(singalhistnominal, pvalues_resolved, pvalues_merged, acceptance, 1)
    # print(" Smearning systematics...")
    # outdicsys = multiprocess(singalhissys, pvalues_resolved, pvalues_merged, acceptance, 1)
    # print(r"Saving 1% smearing output..")
    # savehist({**outdicnominal, **nominalHist}, {**outdicsys, **sysHist}, "run2dblw1")

    # print(r"Doing 2% smearing..")
    # print(" Smearning nominal...")
    # outdicnominal = multiprocess(singalhistnominal, pvalues_resolved, pvalues_merged, acceptance, 2)
    # print(" Smearning systematics...")
    # outdicsys = multiprocess(singalhissys, pvalues_resolved, pvalues_merged, acceptance, 2)
    # print(r"Saving 2% smearing output..")
    # savehist({**outdicnominal, **nominalHist}, {**outdicsys, **sysHist}, "run2dblw2")

    # print(r"Doing 5% smearing..")
    # print(" Smearning nominal...")
    # outdicnominal = multiprocess(singalhistnominal, pvalues_resolved, pvalues_merged, acceptance, 5)
    # print(" Smearning systematics...")
    # outdicsys = multiprocess(singalhissys, pvalues_resolved, pvalues_merged, acceptance, 5)
    # print(r"Saving 5% smearing output..")
    # savehist({**outdicnominal, **nominalHist}, {**outdicsys, **sysHist}, "run2dblw5")

    print(r"Doing 10% smearing..")
    print(" Smearning nominal...")
    outdicnominal = multiprocess(singalhistnominal, pvalues_resolved, pvalues_merged, acceptance, 10)
    print(" Smearning systematics...")
    outdicsys = multiprocess(singalhissys, pvalues_resolved, pvalues_merged, acceptance, 10)
    print(r"Saving 10% smearing output..")
    savehist({**outdicnominal, **nominalHist}, {**outdicsys, **sysHist}, "run2dblw10")



    # sysFile = nominalFile.Get("Systematics")
    # sysHist = loadhist(sysFile)
    # p_tem = [pvalues_resolved[400][1][0][0], pvalues_resolved[400][1][1][0], pvalues_resolved[400][1][2][0], pvalues_resolved[400][1][3][0]]
    # testout = smearing(p_tem, singalhistnominal['AZhllbb400_1tag2pjet_0ptv_SR_mVH'], "BW", 'AZhllbb400_1tag2pjet_0ptv_SR_mVH', 'AZhllbbsmear400_1tag2pjet_0ptv_SR_mVH', 1)
    # c = ROOT.TCanvas()
    # testout.Draw()
    # c.Print("testhist.pdf")
    # # c2 = ROOT.TCanvas()
    # # nominalHist['AZhllbb400_1tag2pjet_0ptv_SR_mVH'].Draw()
    # c.Print("testhistnominal.pdf")
    # # print(nominalHist.keys())

if __name__ == "__main__":
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    # ROOT.RooTrace.active(1)
    main()

