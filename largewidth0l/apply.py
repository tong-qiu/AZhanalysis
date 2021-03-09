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

# x = ROOT.RooRealVar("x1","mA" , -2000, 4000.)

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
    total = outhist.Integral()
    if total > 0:
        scale = totalweight/total
    else:
        scale = 0
    outhist.Scale(scale)
    return outhist

def runMC2(pdf, x, hist, newname, totalweight):
    outhist = hist.Clone(newname)
    outhist.Reset()
    nbins = outhist.GetNbinsX()
    for i in range(nbins):
        x.setVal(outhist.GetBinCenter(i))
        outhist.SetBinContent(i, pdf.getVal())
    # print(newname, outhist.Integral())
    total = outhist.Integral()
    if total > 0:
        scale = totalweight/total
    else:
        scale = 0
    outhist.Scale(scale)
    # for i in range(nbins):
    #     if outhist.GetBinContent(i) < 0.5:
    #         # print("here", outhist.GetBinContent(i), outhist.GetBinCenter(i))
    #         outhist.SetBinContent(i, 0)
    # scale = totalweight/outhist.Integral()
    # outhist.Scale(scale)
    return outhist

def smearing(ps, nominalhist, smearingtype, name, newname, SF, outdic=None):
    hup, hdown = getupdownhsit(nominalhist, newname + "before")
    # x = ROOT.RooRealVar("x1" + newname,"mA" , -2000, 4000.)
    x = ROOT.RooRealVar("x1","mA" , -2000, 4000.)
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

    outhist = runMC2(outpdf, x, nominalhist, name, sum_nominal)
    outhistup = runMC2(outpdfup, x, hup, newname + "up", sum_up)
    outhistdown = runMC2(outpdfdown, x, hdown, newname + "down", sum_down)

    nbins = nominalhist.GetNbinsX()
    for i in range(nbins):
        upvalue = outhistup.GetBinContent(i)
        downvalue = outhistdown.GetBinContent(i)
        error = abs(upvalue - downvalue) / 2.
        if outhist.GetBinContent(i) - error < 0:
            error = outhist.GetBinContent(i)
        outhist.SetBinError(i, error)
    if outdic is not None:
        outdic[name] = outhist
    return outhist

def smearing_batch(ps_b, nominalhist_b, smearingtype_b, name_b, newname_b, SF_b, outdic=None):
    for ps, nominalhist, smearingtype, name, newname, SF in zip(ps_b, nominalhist_b, smearingtype_b, name_b, newname_b, SF_b):
        hup, hdown = getupdownhsit(nominalhist, newname + "before")
        # x = ROOT.RooRealVar("x1" + newname,"mA" , -2000, 4000.)
        x = ROOT.RooRealVar("x1","mA" , -2000, 4000.)
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

        outhist = runMC2(outpdf, x, nominalhist, name, sum_nominal)
        outhistup = runMC2(outpdfup, x, hup, newname + "up", sum_up)
        outhistdown = runMC2(outpdfdown, x, hdown, newname + "down", sum_down)

        nbins = nominalhist.GetNbinsX()
        for i in range(nbins):
            upvalue = outhistup.GetBinContent(i)
            downvalue = outhistdown.GetBinContent(i)
            error = abs(upvalue - downvalue) / 2.
            if outhist.GetBinContent(i) - error < 0:
                error = outhist.GetBinContent(i)
            outhist.SetBinError(i, error)
        if outdic is not None:
            outdic[name] = outhist
        else:
            print("error 1")
            exit(1)
    # return outhist

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
                print("Error: cannot find parameters for " + each_key)

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
            if (i+1)%10 == 0:
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

def multiprocess_batch(hists, pvalues_resolved, pvalues_merged, SFs, width):
    ncore = 10
    nbatch = 200
    manager = multiprocessing.Manager()
    all_sample = manager.dict()
    processes = []
    missinghist = 0
    with tqdm.tqdm(total=len(hists)) as pbar:
        p_tem_b = []
        hists_b = []
        smearingtype_b = []
        each_key_b = []
        new_name_b = []
        sf_tem_b = []
        for i, each_key in enumerate(hists.keys()):
            mass = int(each_key.split("_")[0].split("AZhllbb")[1])
            if mass not in pvalues_merged or mass not in pvalues_resolved or str(mass) not in SFs:
                missinghist += 1
                print("Error: cannot find parameters for " + each_key)

            new_name = each_key.replace("AZhllbb", "AZhllbb" + str(width))
            if "1pfat0pjet" in new_name:
                sf_tem = SFs[str(mass)][str(width)][1]
                smearingtype = pvalues_merged[mass][width][4]
                p_tem = [pvalues_merged[mass][width][0][0], pvalues_merged[mass][width][1][0], pvalues_merged[mass][width][2][0], pvalues_merged[mass][width][3][0]]
            else:
                sf_tem = SFs[str(mass)][str(width)][0]
                smearingtype = pvalues_resolved[mass][width][4]
                p_tem = [pvalues_resolved[mass][width][0][0], pvalues_resolved[mass][width][1][0], pvalues_resolved[mass][width][2][0], pvalues_resolved[mass][width][3][0]]
            
            p_tem_b.append(p_tem)
            hists_b.append(hists[each_key])
            smearingtype_b.append(smearingtype)
            each_key_b.append(each_key)
            new_name_b.append(new_name)
            sf_tem_b.append(sf_tem)
            if len(p_tem_b) < ncore * nbatch:
                continue
            for i_core in range(ncore):
                lowi = i_core * nbatch
                highi = (i_core + 1) * nbatch + 1
                t = multiprocessing.Process(target=smearing_batch, args=(p_tem_b[lowi:highi], hists_b[lowi:highi], smearingtype_b[lowi:highi], each_key_b[lowi:highi], new_name_b[lowi:highi], sf_tem_b[lowi:highi], all_sample))
                processes.append(t)
                t.start()

            for each in processes:
                each.join()
                pbar.update(nbatch)
            processes = []
            p_tem_b = []
            hists_b = []
            smearingtype_b = []
            each_key_b = []
            new_name_b = []
            sf_tem_b = []
        if len(sf_tem_b) != 0:
            processes = []
            if int(len(sf_tem_b)/ncore):
                t = multiprocessing.Process(target=smearing_batch, args=(p_tem_b, hists_b, smearingtype_b, each_key_b, new_name_b, sf_tem_b, all_sample))
                t.start()
                t.join()
                pbar.update(len(p_tem_b))
            else:
                for i_core in range(ncore):
                    lowi = i_core * int(len(sf_tem_b)/ncore)
                    highi = (i_core + 1) * int(len(sf_tem_b)/ncore) + 1
                    if i_core == ncore - 1:
                        t = multiprocessing.Process(target=smearing_batch, args=(p_tem_b[lowi:], hists_b[lowi:], smearingtype_b[lowi:], each_key_b[lowi:], new_name_b[lowi:], sf_tem_b[lowi:], all_sample))
                    else:
                        t = multiprocessing.Process(target=smearing_batch, args=(p_tem_b[lowi:highi], hists_b[lowi:highi], smearingtype_b[lowi:highi], each_key_b[lowi:highi], new_name_b[lowi:highi], sf_tem_b[lowi:highi], all_sample))
                    processes.append(t)
                    t.start()
                for each in processes:
                    each.join()
                    pbar.update(int(len(sf_tem_b)/ncore))
    
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
    nominalFile = ROOT.TFile("run2dbl_hvt.root","read")
    nominalHist = loadhist(nominalFile)
    sysFile = nominalFile.Get("Systematics")
    sysHist = loadhist(sysFile)
    singalhistnominal = getsignalhist(nominalHist, bbA, debug)
    singalhissys = getsignalhist(sysHist, bbA, debug)

    # for each in singalhistnominal.keys():
    #     print(each)

    print(r"Doing 1% smearing..")
    print(" Smearning nominal...")
    outdicnominal = multiprocess_batch(singalhistnominal, pvalues_resolved, pvalues_merged, acceptance, 1)
    print(" Smearning systematics...")
    outdicsys = multiprocess_batch(singalhissys, pvalues_resolved, pvalues_merged, acceptance, 1)
    print(r"Saving 1% smearing output..")
    savehist({**outdicnominal, **nominalHist}, {**outdicsys, **sysHist}, "run2dblggAw1")

    print(r"Doing 2% smearing..")
    print(" Smearning nominal...")
    outdicnominal = multiprocess_batch(singalhistnominal, pvalues_resolved, pvalues_merged, acceptance, 2)
    print(" Smearning systematics...")
    outdicsys = multiprocess_batch(singalhissys, pvalues_resolved, pvalues_merged, acceptance, 2)
    print(r"Saving 2% smearing output..")
    savehist({**outdicnominal, **nominalHist}, {**outdicsys, **sysHist}, "run2dblggAw2")

    print(r"Doing 5% smearing..")
    print(" Smearning nominal...")
    outdicnominal = multiprocess_batch(singalhistnominal, pvalues_resolved, pvalues_merged, acceptance, 5)
    print(" Smearning systematics...")
    outdicsys = multiprocess_batch(singalhissys, pvalues_resolved, pvalues_merged, acceptance, 5)
    print(r"Saving 5% smearing output..")
    savehist({**outdicnominal, **nominalHist}, {**outdicsys, **sysHist}, "run2dblggAw5")

    print(r"Doing 10% smearing..")
    print(" Smearning nominal...")
    outdicnominal = multiprocess_batch(singalhistnominal, pvalues_resolved, pvalues_merged, acceptance, 10)
    print(" Smearning systematics...")
    outdicsys = multiprocess_batch(singalhissys, pvalues_resolved, pvalues_merged, acceptance, 10)
    print(r"Saving 10% smearing output..")
    savehist({**outdicnominal, **nominalHist}, {**outdicsys, **sysHist}, "run2dblggAw10")

    print(r"Doing 20% smearing..")
    print(" Smearning nominal...")
    outdicnominal = multiprocess_batch(singalhistnominal, pvalues_resolved, pvalues_merged, acceptance, 20)
    print(" Smearning systematics...")
    outdicsys = multiprocess_batch(singalhissys, pvalues_resolved, pvalues_merged, acceptance, 20)
    print(r"Saving 20% smearing output..")
    savehist({**outdicnominal, **nominalHist}, {**outdicsys, **sysHist}, "run2dblggAw20")


    # print("testing")
    # p_tem = [pvalues_resolved[400][1][0][0], pvalues_resolved[400][1][1][0], pvalues_resolved[400][1][2][0], pvalues_resolved[400][1][3][0]]
    # testout = smearing(p_tem, singalhistnominal['AZhllbb400_1tag2pjet_0ptv_SR_mVH'], "BW", 'AZhllbb400_1tag2pjet_0ptv_SR_mVH', 'AZhllbbsmear400_1tag2pjet_0ptv_SR_mVH', 1)
    # c = ROOT.TCanvas()
    # testout.Draw()
    # c.Print("testhist.pdf")
    # c2 = ROOT.TCanvas()
    # singalhistnominal['AZhllbb400_1tag2pjet_0ptv_SR_mVH'].Draw()
    # c2.Print("testhistnominal.pdf")
    # print(nominalHist.keys())

if __name__ == "__main__":
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    # ROOT.RooTrace.active(1)
    main()

