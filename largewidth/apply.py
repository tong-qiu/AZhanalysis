import ROOT
import numpy as np
import os
import pickle
ROOT.gROOT.ProcessLine(".L My_modified_bw.cxx+")
ROOT.gInterpreter.ProcessLine('#include "My_modified_bw.h"')
ROOT.gSystem.Load('My_modified_bw_cxx.so')
from ROOT import My_modified_bw
import math

x = ROOT.RooRealVar("x1","mA" , -2000, 4000.)

def loadhist(inFile, filter):
    out = {}
    for key in inFile.GetListOfKeys():
        tem_hist = key.ReadObj()
        if isinstance(tem_hist, ROOT.TH1F):
            temName = key.GetName()
            if temName[0:len(filter)] == filter:
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
    data = pdf.generate(ROOT.RooArgSet(x), 1000000)
    for i in range(0, 1000000-1):
        outhist.Fill(data.get(i).getRealValue("x1"))
        scale = totalweight/outhist.Integral()
        outhist.Scale(scale)
    return outhist

def smearing(ps, nominalhist, smearingtype, name, newname, SF):
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

    outhist = runMC(outpdf, x, nominalhist, newname, sum_nominal)
    outhistup = runMC(outpdfup, x, hup, newname + "up", sum_up)
    outhistdown = runMC(outpdfdown, x, hdown, newname + "down", sum_down)

    nbins = nominalhist.GetNbinsX()
    for i in range(nbins):
        upvalue = outhistup.GetBinContent(i)
        downvalue = outhistdown.GetBinContent(i)
        error = abs(upvalue - downvalue) / 2.
        outhist.SetBinError(i, error)
    return outhist


def main():
    with open("fitvalues_resolved.pickle",'rb') as f:
        pvalues = pickle.load(f)
        # print(pvalues)

    nominalFile = ROOT.TFile("run2dbl_nosys.root","read")
    nominalHist = loadhist(nominalFile, "AZhllbb")
    # sysFile = nominalFile.Get("Systematics")
    # sysHist = loadhist(sysFile)
    p_tem = [pvalues[400][1][0][0], pvalues[400][1][1][0], pvalues[400][1][2][0], pvalues[400][1][3][0]]
    testout = smearing(p_tem, nominalHist['AZhllbb400_1tag2pjet_0ptv_SR_mVH'], "BW", 'AZhllbb400_1tag2pjet_0ptv_SR_mVH', 'AZhllbbsmear400_1tag2pjet_0ptv_SR_mVH', 1)
    c = ROOT.TCanvas()
    testout.Draw()
    c.Print("testhist.pdf")
    c2 = ROOT.TCanvas()
    nominalHist['AZhllbb400_1tag2pjet_0ptv_SR_mVH'].Draw()
    c2.Print("testhistnominal.pdf")
    # print(nominalHist.keys())

if __name__ == "__main__":
    main()
