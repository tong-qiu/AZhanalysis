from numpy.core.function_base import linspace
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
import os
import pickle
import multiprocessing

ROOT.gROOT.ProcessLine(".L My_modified_bw.cxx+")
ROOT.gInterpreter.ProcessLine('#include "My_modified_bw.h"')
ROOT.gSystem.Load('My_modified_bw_cxx.so')
# ROOT.gSystem.Load("My_modified_bw.cxx")
from ROOT import My_modified_bw

def ptll_cut(mvh):
    if mvh < 320.:
        return 0.
    return 20. + 9. * pow(mvh - 320., 0.5)

def resolvedselection(entry):
    if entry.nsigjet < 2:
        return False
    if entry.nbjets == 0:
        return False
    if entry.ptb1/1000. < 45:
        return False
    if entry.mbbres/1000. < 110:
        return False
    if entry.mbbres/1000. > 140:
        return False
    if entry.met/1000. < 150:
        return False
    if abs(entry.phiBB) > 7. * math.pi / 9.:
        return False
    if entry.phiHMETres < 2. * math.pi / 3.:
        return False
    if entry.minJetMET < math.pi / 9 and (entry.nsigjet == 2 or entry.nsigjet == 3):
        return False
    if entry.minJetMET < math.pi / 6 and entry.nsigjet >= 4:
        return False
    return True

def mergedselection(entry):
    if entry.ptfj1/1000. < 250:
        return False
    if entry.mfj1/1000. < 75:
        return False
    if abs(entry.etafj1) > 2:
        return False
    if entry.mfj1/1000. > 145:
        return False

    if entry.met/1000. < 200:
        return False
    if entry.phiHMETmerg < 2. * math.pi / 3.:
        return False
    # if entry.minJetMET < math.pi / 9 and (entry.nsigjet == 2 or entry.nsigjet == 3):
    #     return False
    # if entry.minJetMET < math.pi / 6 and entry.nsigjet >= 4:
    #     return False
    return True


def histplot_withsub_raw(datas, bins, weights=None, usererror = None, labels = None, scale=1., removenorm = None, **kwargs,):
    settings = {
        "xlabel" : r"$m_{Vh}[GeV]$",
        "ylabel": 'Number of Events',
        "title1": r"$\mathbf{ATLAS}$",# \newline Ptl next-leading, full cuts, 2 b-tags $",
        "title1_1": r"$\mathit{Internal}$",
        "title2": r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$",# Ptl next-leading, full cuts, 2 b-tags $",
        #"title3": r"$\mathbf{2\;lep.,2\;b-tag}$",
        "title3": "0 lep., 2 b-tag",
        "filename": "deltatest2",
        "log_y":False,
        "norm":False,
        "central":"none",
        "upper_y": 1.5, 
        "do_errorbar": False
        }
    for each_key in kwargs.items():
        settings[each_key[0]] = kwargs[each_key[0]]
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios':[3, 1]}, figsize=(10, 10))
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.06)

    if weights is None:
        weights = []
        for each in datas:
            weights.append(np.ones(len(each)))

    if removenorm:
        for i in range(len(weights)):
            totalweightbefore = np.sum(weights[i])
            weights[i] = np.array(weights[i]) / np.sum(weights[i])
            if usererror is not None:
                if usererror[i] is not None:
                    usererror[i] = usererror[i]/totalweightbefore

    sigmas = []
    weight_in_binses = []
    for i in range(len(datas)):
        datas[i] = np.array(datas[i])
        event_location = np.digitize(datas[i]/scale, bins)
        sigma2 = []
        weight_in_bins = []
        for j in range(np.size(bins) - 1):
            bin_weight = weights[i][np.where(event_location == j+1)[0]]
            if not (usererror is not None and usererror[i] is not None):
                sigma2.append(np.sum(bin_weight**2.))
            weight_in_bins.append(np.sum(bin_weight))
            if usererror is not None:
                if usererror[i] is not None:
                    binederror = usererror[i][np.where(event_location == j+1)[0]]
                    sigma2.append(np.sum(binederror**2))

        sigmas.append(np.array(sigma2)**0.5)
        weight_in_binses.append(np.array(weight_in_bins))

    colors = ['b', 'g', 'r', 'c', 'm', 'y']
    # print(np.array(datas).shape, np.array(weights).shape)
    # print(np.array(datas[1]).shape, np.array(weights[1]).shape)
    # print(np.array(datas[0]).shape, np.array(weights[0]).shape)
    for i in range(len(datas)):
        datas[i] = np.array(datas[i])/scale
    ax1.hist(datas, bins, histtype='step', fill=False, color=colors[0:len(datas)], weights=weights)
    bins = np.array(bins)
    bin_centre = []
    for i in range(len(datas)):
        bin_centre = (bins[0:-1] + bins[1:])/2
        ax1.errorbar(bin_centre, weight_in_binses[i], xerr=0.0001, yerr=sigmas[i], fmt='.', color=colors[i], label=str(labels[i]))
    ax1.legend(loc='upper right',prop={'size': 20}, frameon=False)
    
    ymin, ymax = ax1.get_ylim()
    ax1.set_ylim([0,ymax* settings["upper_y"]])
    ax1.text(0.05, 1.55 / 1.7, settings['title1'], fontsize=25, transform=ax1.transAxes, style='italic', fontweight='bold')
    ax1.text(0.227, 1.55/ 1.7, settings['title1_1'], fontsize=25, transform=ax1.transAxes)
    ax1.text(0.05, 1.40 / 1.7, settings['title2'], fontsize=20, transform=ax1.transAxes, style='italic', fontweight='bold')
    ax1.text(0.05, 1.26 / 1.7, settings['title3'], fontsize=18, weight='bold', style='italic', transform=ax1.transAxes)
    ax1.set_ylabel(settings['ylabel'], fontsize=20)
    if settings['log_y']:
        ax1.set_yscale('log')
        ax1.set_ylim([0.1, 10**(math.log10(ymax) * settings["upper_y"])])
        ax1.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10,numticks=100))
        ax1.minorticks_on()
    ax1.get_xaxis().set_ticks([])

    i = 0
    centerheight = []
    for each_height, each_label in zip(weight_in_binses, labels):
        if i == 0:
            centerheight = np.array(each_height)
        for j in range(len(centerheight)):
            if centerheight[j] == 0:
                centerheight[j] = 0.00001
        new_each_height = np.array(each_height)/centerheight
        if i != 0:
            new_each_height[np.isnan(new_each_height)] = -10
            new_each_height[np.isposinf(new_each_height)] = 2
            new_each_height[np.isneginf(new_each_height)] = -2
        else:
            new_each_height[np.isnan(new_each_height)] = 1
            new_each_height[np.isinf(new_each_height)] = 1
        ax2.hist(bin_centre, bins, weights=new_each_height-1, label=each_label, histtype=u'step', color=colors[i])
        i += 1
    ax2.set_ylim([-1, 1])
    ax2.plot([bins[0], bins[np.size(bins)-1]], [0, 0], linestyle='--', color='k')
    ax2.set_ylabel(settings["central"] + "/MC" + "-1", fontsize=20)
    ax2.set_xlabel(settings['xlabel'], fontsize=20)

    fig.savefig(settings['filename'] + '.pdf', bbox_inches='tight', pad_inches = 0.25)
    plt.close(fig)

xg = ROOT.RooRealVar("x1","mA", -2000, 4000.)

class loadroot():
    def __init__(self, path):
        self.path = path
        self.ptl1 = []
        self.ptl2 = []
        self.ptb1 = []
        self.ma = []
        self.path = path
        self.name = path.split("/")[-1].replace(".root", "")
        self.mass = float(self.name.split("w")[0])
        self.loadfile()
    
    def loadfile(self):
        self.hist = ROOT.TH1F("hist" + self.name, "", 1000, -2000, 4000.)
        self.hist.Sumw2()
        print("loading " + self.path + " ...")
        f = ROOT.TFile(self.path)
        t1 = f.Get("data")
        for entry in t1:
            # if resolvedselection(entry):
            #     self.hist.Fill(entry.mVHres/1000.)

            if resolvedselection(entry):
                continue
            elif mergedselection(entry):
                self.hist.Fill(entry.mVHmerg/1000.)   
        f.Close()
        # self.x = ROOT.RooRealVar("x1" + self.name,"mA", -2000, 4000.)
        self.x = xg
        self.datahist = ROOT.RooDataHist("datahist" + self.name, "datahist", ROOT.RooArgList(self.x), self.hist)
        self.pdf = ROOT.RooHistPdf("pdf" + self.name, "your pdf", ROOT.RooArgSet(self.x), self.datahist)
    
    def getDatahist(self):
        return self.datahist

    def gethist(self):
        # hist = self.pdf.createHistogram("x1", 1000)
        hist = self.hist
        nbin = hist.GetNbinsX()
        center = []
        height = []
        error = []
        for i in range(nbin):
            center.append(hist.GetBinCenter(i))
            height.append(hist.GetBinContent(i))
            error.append(hist.GetBinError(i))
        return(np.array(center), np.array(height), np.array(error))

    def convbw(self, width):
        widthp = int(width * 100)
        width = self.mass * (width)
        self.x.setBins(1000, "cache")
        self.bwm0 = ROOT.RooRealVar("bwm0" + self.name + str(widthp), "bwm0", 0)
        self.bww = ROOT.RooRealVar("bwwidth" + self.name + str(widthp), "bwwidth", width)
        self.bworiginal = ROOT.RooBreitWigner("bworiginal" + self.name + str(widthp), "bworiginal", self.x, self.bwm0, self.bww)
        out = ROOT.RooFFTConvPdf("outbw" + self.name + str(widthp), "outbw", self.x, self.bworiginal, self.pdf)
        # data= out.generate(ROOT.RooArgSet(self.x), 1000000)
        # outlist = []
        # for i in range(0, 1000000-1):
        #     outlist.append(data.get(i).getRealValue("x1"))
        return out
    

    def convbwmodified(self, width):
        widthp = int(width * 100)
        width = self.mass * (width)
        self.x.setBins(1000, "cache")
        # self.m0 = ROOT.RooRealVar("m0" + self.name + str(widthp), "m0", -self.mass)
        # self.m0.setConstant()
        self.w = ROOT.RooRealVar("width" + self.name + str(widthp), "width", width)
        self.w.setConstant()
        # self.massv = ROOT.RooRealVar("massv" + self.name + str(widthp), "massv", self.mass)
        # self.massv.setConstant()

        # if widthp >= 5: # merged
        #     # self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(widthp), "lnm0", 0, 1050)
        #     self.lnk = ROOT.RooRealVar("lnmk" + self.name + str(widthp), "lnmk", 1.001, 1.7)
        #     if widthp > 10:
        #         self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(widthp), "lnm0", 0, 500)
        #     else:
        #         self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(widthp), "lnm0", 0, 2000)
        # elif widthp == 1:
        #     self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(widthp), "lnm0", 0, 2000)
        #     self.lnk = ROOT.RooRealVar("lnmk" + self.name + str(widthp), "lnmk", 1.001, 1.2)
        # else:

        # resolved
        # if widthp == 5:
        #     self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(widthp), "lnm0", 0, 2000)
        #     self.lnk = ROOT.RooRealVar("lnmk" + self.name + str(widthp), "lnmk", 1.001, 2.2)
        # elif widthp == 20:
        #     if self.mass > 375:
        #         self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(widthp), "lnm0", 0, 600)
        #     else:
        #         self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(widthp), "lnm0", 0, 475)
        #     self.lnk = ROOT.RooRealVar("lnmk" + self.name + str(widthp), "lnmk", 1.001, 2.2)
        # else:
        #     self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(widthp), "lnm0", 0, 2000)
        #     self.lnk = ROOT.RooRealVar("lnmk" + self.name + str(widthp), "lnmk", 1.001, 30)


        # merged
        # if widthp == 1:
        #     print(widthp, "here", self.mass)
        #     self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(widthp), "lnm0", 0, self.mass*0.3)
        # elif widthp == 2:
        #     self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(widthp), "lnm0", 0, self.mass*0.5)
        # else:
        #     self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(widthp), "lnm0", 0, self.mass)
        self.lnk = ROOT.RooRealVar("lnmk" + self.name + str(widthp), "lnmk", 1.001, 30)
        if widthp < 6:
            self.m0 = ROOT.RooRealVar("m0" + self.name + str(widthp), "m0", -self.mass * 0.8, self.mass * widthp * 2 / 100)
        else: 
            self.m0 = ROOT.RooRealVar("m0" + self.name + str(widthp), "m0", -self.mass, -self.mass * 0.5)
        self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(widthp), "lnm0", 0, 3000)
        self.bw = My_modified_bw("bwmodified" + self.name + str(widthp), "bwmodified", self.x, self.w, self.lnm0, self.lnk, self.m0)

        # self.m0test = ROOT.RooRealVar("m0test" + self.name + str(width), "m0test", 0)
        # self.testbw = ROOT.RooBreitWigner("bwtest" + self.name + str(width), "bwtest", self.x, self.m0test, self.w)

        out = ROOT.RooFFTConvPdf("out" + self.name + str(widthp), "out", self.x, self.pdf, self.bw)

        # c = ROOT.TCanvas("rf208_convolution", "rf208_convolution", 600, 600)
        # ROOT.gPad.SetLeftMargin(0.15)
        # # xframe = glabalx.frame()
        # xframe = self.x.frame(ROOT.RooFit.Title("landau (x) gauss convolution"))
        # xframe.GetYaxis().SetTitleOffset(1.4)
        # #datahist.plotOn(xframe)
        # out.plotOn(xframe)
        # xframe.Draw()
        # # input()
        # c.SaveAs("test.png")

        return out

    def convbw(self, width):
        widthp = int(width * 100)
        width = self.mass * (width)
        self.x.setBins(1000, "cache")
        self.bwm0 = ROOT.RooRealVar("bwm0" + self.name + str(widthp), "bwm0", 0)
        self.bww = ROOT.RooRealVar("bwwidth" + self.name + str(widthp), "bwwidth", width)
        self.bworiginal = ROOT.RooBreitWigner("bworiginal" + self.name + str(widthp), "bworiginal", self.x, self.bwm0, self.bww)
        out = ROOT.RooFFTConvPdf("outbw" + self.name + str(widthp), "outbw", self.x, self.bworiginal, self.pdf)
        # data= out.generate(ROOT.RooArgSet(self.x), 1000000)
        # outlist = []
        # for i in range(0, 1000000-1):
        #     outlist.append(data.get(i).getRealValue("x1"))
        return out
    

    def convbwmodifiedp(self, width, pm0, plnm0, plnk, i):
        widthp = int(width * 100)
        width = self.mass * (width)
        self.x.setBins(1000, "cache")
        self.w = ROOT.RooRealVar("width" + self.name + str(widthp) + str(i), "width", width)
        self.w.setConstant()
        self.lnk = ROOT.RooRealVar("lnmk" + self.name + str(widthp) + str(i), "lnmk", plnk, 1.001, 30)
        # if widthp < 6:
        #     self.m0 = ROOT.RooRealVar("m0" + self.name + str(widthp) + str(i), "m0", pm0, -self.mass * 0.8, self.mass * widthp * 2 / 100)
        # else: 
        self.m0 = ROOT.RooRealVar("m0" + self.name + str(widthp) + str(i), "m0", pm0, -self.mass, self.mass * 0.5)
        self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(widthp) + str(i), "lnm0", plnm0, 0, 3000)
        self.bw = My_modified_bw("bwmodified" + self.name + str(widthp) + str(i), "bwmodified", self.x, self.w, self.lnm0, self.lnk, self.m0)

        # self.m0test = ROOT.RooRealVar("m0test" + self.name + str(width), "m0test", 0)
        # self.testbw = ROOT.RooBreitWigner("bwtest" + self.name + str(width), "bwtest", self.x, self.m0test, self.w)

        out = ROOT.RooFFTConvPdf("out" + self.name + str(widthp) + str(i), "out", self.x, self.pdf, self.bw)


        return out

def bw(x, width):
    return 1/(x * x + 0.25 * width * width)

def lognorm(x, m0, k, mass):
    ln_k =  ROOT.TMath.Abs(ROOT.TMath.Log(k))
    ln_m0 = ROOT.TMath.Log(m0)
    return ROOT.Math.lognormal_pdf(x, ln_m0, ln_k, mass)

def getloss(data1, w1, error, data2, w2, binning):
    bin_heights1, bin_edges1 = np.histogram(data1, bins=binning, weights=w1)
    bin_heights2, bin_edges2 = np.histogram(data2, bins=binning, weights=w2)
    bin_heights1 = bin_heights1 / sum(bin_heights1)
    bin_heights2 = bin_heights2 / sum(bin_heights2)

    # event_location = np.digitize(data1, binning)
    # error_in_bins = []
    # for j in range(np.size(binning) - 1):
    #     bin_error = error[np.where(event_location == j+1)[0]]
    #     errorsum = np.sum(bin_error**2)
    #     if errorsum == 0:
    #         errorsum = 1
    #     error_in_bins.append(errorsum)
    # return sum((bin_heights1 - bin_heights2) ** 2/error_in_bins) ** 0.5
    return np.sum((bin_heights1 - bin_heights2) ** 2) ** 0.5

def crazyfit(NWidthobj, datahist, width, center, height, error, bins, mass):
    mass = int(mass)
    # if width < 6:
    #     lowm0 = -mass * 0.8
    #     highm0 = mass * width * 2 / 1000.
    # else:
    lowm0 = -mass
    highm0 = mass * 0.5
    # outpdf1 = NWidthobj.convbwmodified(width/100.)
    minloss = 99999999999
    optlnk = 0
    optlnm0 = 0
    optm0 = 0
    optdatalist = []
    i = 0
    for each_lnk in np.linspace(1.1, 20, 6):
        for each_lnm0 in np.linspace(100, 2000-1, 8):
            for each_m0 in np.linspace(lowm0+10, highm0-10, 8):
                i += 1
                # NWidthobj.lnk.setVal(each_lnk, 1.001, 30)
                # NWidthobj.lnm0.setVal(each_lnm0, 0, 3000)
                # NWidthobj.m0.setVal(each_m0, lowm0, highm0)
                outpdf1 = NWidthobj.convbwmodifiedp(width, each_m0, each_lnm0, each_lnk, i)
                outpdf1.fitTo(datahist)
                datalist1 = []
                dataerror = []
                for eachc in center:
                    dataerror.append(0)
                    if eachc > 0:
                        NWidthobj.x.setVal(eachc)
                        datalist1.append(outpdf1.getVal())
                    else:
                        datalist1.append(0)
                datalist1 = np.array(datalist1) * ( np.sum(height) / np.sum(datalist1))
                loss1 = getloss(center, height, error, center, datalist1, bins)
                # exit(1)
                if loss1 < minloss:
                    # print(datalist1)
                    # print(height)
                    minloss = loss1
                    
                    optdatalist = datalist1
                    optlnk = NWidthobj.lnk.getValV()
                    optlnm0 = NWidthobj.lnm0.getValV()
                    optm0 = NWidthobj.m0.getValV() 
                    # print(minloss, optlnk, optlnm0, optm0)   
    NWidthobj.lnk.setVal(optlnk)
    NWidthobj.lnm0.setVal(optlnm0)
    NWidthobj.m0.setVal(optm0)
    return optdatalist, minloss, dataerror

def plot_mass_modified(mass, bins):
    mass = str(mass)
    NWidthobj = loadroot("ntuples/" + mass + "w0.root")
    outdic = {}
    for each_width in [1, 2, 5, 10, 20]:
        LWidthobj = loadroot("ntuples/" + mass + "w" + str(each_width) + ".root")
        datahist = LWidthobj.getDatahist()
        center, height, error = LWidthobj.gethist()

        # modified bw
        # outpdf1 = NWidthobj.convbwmodified(each_width/100.)
        # outpdf1.fitTo(datahist)
        # datalist1 = []
        # dataerror = []
        # for eachc in center:
        #     dataerror.append(0)
        #     if eachc > 0:
        #         NWidthobj.x.setVal(eachc)
        #         datalist1.append(outpdf1.getVal())
        #     else:
        #         datalist1.append(0)
        # datalist1 = np.array(datalist1) * ( np.sum(height) / np.sum(datalist1))
        # loss1 = getloss(center, height, error, center, datalist1, bins)
        datalist1, loss1, dataerror = crazyfit(NWidthobj, datahist, each_width/100., center, height, error, bins, mass)

        # original bw
        outpdf2 = NWidthobj.convbw(each_width/100.)
        datalist2 = []
        for eachc in center:
            if eachc > 0:
                NWidthobj.x.setVal(eachc)
                datalist2.append(outpdf2.getVal())
            else:
                datalist2.append(0)
        datalist2 = np.array(datalist2) * ( np.sum(height) / np.sum(datalist2))
        loss2 = getloss(center, height, error, center, datalist2, bins)
        # print("BW loss ", loss2)



        # # modified bw
        # outpdf1 = NWidthobj.convbwmodified(each_width/100.)
        # outpdf1.fitTo(datahist)
        # data1 = outpdf1.generate(ROOT.RooArgSet(NWidthobj.x), 1000000)
        # datalist1 = []
        # for i in range(0, 1000000-1):
        #     datalist1.append(data1.get(i).getRealValue("x1"))
        # loss1 = getloss(center, height, error, datalist1, [1]*len(datalist1), bins)
        # # ROOT.RooTrace.printObjectCounts()
        # # original bw
        # outpdf2 = NWidthobj.convbw(each_width/100.)
        # # outpdf2.fitTo(datahist)
        # data2 = outpdf2.generate(ROOT.RooArgSet(NWidthobj.x), 1000000)
        # datalist2 = []
        # for i in range(0, 1000000-1):
        #     datalist2.append(data2.get(i).getRealValue("x1"))
        # loss2 = getloss(center, height, error, datalist2, [1]*len(datalist2), bins)

        if loss2 >= loss1:
            method = "modified BW"
            altmethod = "BW"
            datalist = datalist1
            datalistalt = datalist2
        if loss2 < loss1:
            method = "BW"
            altmethod = "modified BW"
            datalist = datalist2
            datalistalt = datalist1

        # histplot_withsub_raw([center, datalist], bins, weights=[np.array(height), np.array([1]*len(datalist))], usererror=[error, None], 
        #                         labels = ["MC", "smearing"], removenorm=True, filename="output/com"+mass+"w"+str(each_width),
        #                         title2=r"$\mathit{\sqrt{s}=13\:TeV}$ ", title3="mass = "+mass+" GeV, W = " + str(each_width) + "%, " + method + " smearing", central="smearing")
        # histplot_withsub_raw([center, datalistalt], bins, weights=[np.array(height), np.array([1]*len(datalistalt))], usererror=[error, None], 
        #                         labels = ["MC", "smearing"], removenorm=True, filename="output/altcom"+mass+"w"+str(each_width),
        #                         title2=r"$\mathit{\sqrt{s}=13\:TeV}$ ", title3="mass = "+mass+" GeV, W = " + str(each_width) + "%, " + altmethod + " smearing", central="smearing")
        histplot_withsub_raw([center, center], bins, weights=[np.array(height), datalist], usererror=[error, dataerror], 
                                labels = ["MC", "smearing"], removenorm=True, filename="output/com"+mass+"w"+str(each_width),
                                title2=r"$\mathit{\sqrt{s}=13\:TeV}$ ", title3="mass = "+mass+" GeV, W = " + str(each_width) + "%, " + method + " smearing", central="smearing")
        histplot_withsub_raw([center, center], bins, weights=[np.array(height), datalistalt], usererror=[error, dataerror], 
                                labels = ["MC", "smearing"], removenorm=True, filename="output/altcom"+mass+"w"+str(each_width),
                                title2=r"$\mathit{\sqrt{s}=13\:TeV}$ ", title3="mass = "+mass+" GeV, W = " + str(each_width) + "%, " + altmethod + " smearing", central="smearing")
        # pdfplot
        bwlist = []
        lognormlist = []
        # getError () 
        v_width = NWidthobj.w.getValV()
        v_m0 = NWidthobj.m0.getValV()
        v_lnm0 = NWidthobj.lnm0.getValV()
        v_lnk = NWidthobj.lnk.getValV()

        e_width = NWidthobj.w.getError()
        e_m0 = NWidthobj.m0.getError()
        e_lnm0 = NWidthobj.lnm0.getError()
        e_lnk = NWidthobj.lnk.getError()
        for each_x in range(-2000, 1000):
            bwlist.append(bw(each_x, v_width))
            lognormlist.append(lognorm(each_x, v_lnm0, v_lnk, v_m0))
        fig, (ax1) = plt.subplots(figsize=(10, 10))
        ax1.plot(range(-2000, 1000), bwlist, label="Breit–Wigner")
        ax1.plot(range(-2000, 1000), lognormlist, label="Log-normal")
        ymin, ymax = ax1.get_ylim()
        factor = int(ymax / max(np.array(lognormlist) * np.array(bwlist)) * 0.8)
        ax1.plot(range(-2000, 1000), np.array(lognormlist) * np.array(bwlist) * factor, label=r"product $\times$ " + str(factor))
        ax1.set_xlabel(r"m_{A} [GeV]", fontsize=20)
        ax1.set_ylabel(r"arbitrary unit", fontsize=20)
        ax1.legend(loc='upper right',prop={'size': 20}, frameon=False)
        ymin, ymax = ax1.get_ylim()
        ax1.set_ylim([0,ymax * 1.4])
        ax1.text(0.05, 1.55 / 1.7, "ATLAS", fontsize=25, transform=ax1.transAxes, style='italic', fontweight='bold')
        ax1.text(0.227, 1.55/ 1.7, "Internal", fontsize=25, transform=ax1.transAxes)
        ax1.text(0.05, 1.40 / 1.7, "m = " + mass + " GeV, width = " + str(each_width) + "%", fontsize=17, transform=ax1.transAxes)
        fig.savefig("output/pdfm" + mass + "w" + str(each_width) + ".pdf", bbox_inches='tight', pad_inches = 0.25)
        plt.close(fig)

        outdic[each_width] = ((v_lnm0, e_lnm0), (v_lnk, e_lnk), (v_width, e_width), (v_m0, e_m0), method)
    with open('pickle/fitvalues' + str(mass) + '.pickle', 'wb') as f:
        pickle.dump(outdic, f)
    return outdic

def main():
    massset = set()
    for each in os.listdir("ntuples"):
        if ".root" in each:
            massset.add(int(each.split("w")[0]))
    i = 0
    processes = []
    for each_mass in sorted(massset):
        if each_mass < 300:
            continue
        maxv = int(each_mass * 2.5)
        if maxv > 3000:
            maxv = 3000
        i += 1
        # plot_mass_modified(each_mass, linspace(0, int(maxv), 25))
        t = multiprocessing.Process(target=plot_mass_modified, args=(each_mass, linspace(0, int(maxv), 25)))
        processes.append(t)
        t.start()
        if (i+1)%8 == 0:
            for each in processes:
                each.join()
            processes = []
    for each in processes:
        each.join()
        

    paras = {}
    allfiles = os.listdir("pickle")
    for each in allfiles:
        mass = int(each.replace("fitvalues", "").replace(".pickle", ""))
        with open("pickle/"  + each,'rb') as f:
            pvalues = pickle.load(f)
            paras[mass] = pvalues
    with open('fitvalues.pickle', 'wb') as f:
        pickle.dump(paras, f)
    # with open('pickle/fitvalues.pickle', 'wb') as f:
    #     pickle.dump(paras, f)

if __name__ == '__main__':
    # ROOT.RooMsgService.instance().setSilentMode(True)
    # ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.FATAL)
    # ROOT.RooTrace.active(1)
    ROOT.RooMsgService.instance().setSilentMode(True)
    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)
    ROOT.gErrorIgnoreLevel = ROOT.kError
    main()
