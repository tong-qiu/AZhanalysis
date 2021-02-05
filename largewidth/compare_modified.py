from numpy.core.function_base import linspace
import ROOT
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
# ROOT.gROOT.ProcessLine(".L My_modified_bw.cxx+")
# ROOT.gInterpreter.ProcessLine('#include "My_modified_bw.h"')
# ROOT.gSystem.Load('My_modified_bw_cxx.so')
# # ROOT.gSystem.Load("My_modified_bw.cxx")
# from ROOT import My_modified_bw

def histplot_withsub_raw(datas, bins, weights=None, usererror = None, labels = None, scale=1., removenorm = None, **kwargs,):
    settings = {
        "xlabel" : r"$m_{Vh}[GeV]$",
        "ylabel": 'Number of Events',
        "title1": r"$\mathbf{ATLAS}$",# \newline Ptl next-leading, full cuts, 2 b-tags $",
        "title1_1": r"$\mathit{Internal}$",
        "title2": r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$",# Ptl next-leading, full cuts, 2 b-tags $",
        #"title3": r"$\mathbf{2\;lep.,2\;b-tag}$",
        "title3": "2 lep., 2 b-tag",
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
    ax1.hist(np.array(datas)/scale, bins, histtype='step', fill=False, color=colors[0:len(datas)], weights=weights)
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
    ax1.text(0.05, 1.40 / 1.7, settings['title2'], fontsize=23, transform=ax1.transAxes, style='italic', fontweight='bold')
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

glabalx = None
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
        f = ROOT.TFile(self.path)
        t1 = f.Get("data")
        for entry in t1:
            if entry.ptl2/1000. < 20:
                continue
            if entry.ptl2/1000. < 27:
                continue
            if entry.ptb1/1000. < 45:
                continue
            self.hist.Fill(entry.mA/1000.)
        self.x = glabalx
        self.datahist = ROOT.RooDataHist("datahist" + self.name, "datahist", ROOT.RooArgList(self.x), self.hist)
        self.pdf = ROOT.RooHistPdf("pdf" + self.name, "your pdf", ROOT.RooArgSet(self.x), self.datahist)
    
    def getDatahist(self):
        return self.datahist

        # xframe= self.x.frame()
        # self.pdf.plotOn(xframe)
        # xframe.Draw()
        # input()

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

    def getwidthdata(self, width):
        width = self.mass * (width)
        self.x.setBins(1000, "cache")
        m0 = ROOT.RooRealVar("m0" + self.name + str(width), "m0", 0)
        w = ROOT.RooRealVar("width" + self.name + str(width), "width", width)
        bw = ROOT.RooBreitWigner("bw" + self.name + str(width), "bw", self.x, m0, w)
        out = ROOT.RooFFTConvPdf("out" + self.name + str(width), "out", self.x, bw, self.pdf)
        data= out.generate(ROOT.RooArgSet(self.x), 1000000)
        outlist = []
        for i in range(0, 1000000-1):
            outlist.append(data.get(i).getRealValue("x1"))
        return np.array(outlist)
    
    # def test(self, width):
    #     width = self.mass * (width)
    #     self.x.setBins(1000, "cache")
    #     self.m0 = ROOT.RooRealVar("m0" + self.name + str(width), "m0", 0)
    #     self.w = ROOT.RooRealVar("width" + self.name + str(width), "width", width)

    #     self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(width), "lnm0", self.mass - width * 2, self.mass + width * 2)
    #     self.lnk = ROOT.RooRealVar("lnmk" + self.name + str(width), "lnmk", 1.1, 100)
    #     self.bw = ROOT.RooBreitWigner("bw" + self.name + str(width), "bw", self.x, self.m0, self.w)
    #     self.ln = ROOT.RooLognormal("ln" + self.name + str(width), "ln", self.x, self.lnm0, self.lnk)
    #     self.modifedbw = ROOT.RooProdPdf("modifedbw" + self.name + str(width), "modifedbw", ROOT.RooArgSet(self.bw, self.ln))
    #     out = ROOT.RooFFTConvPdf("out" + self.name + str(width), "out", self.x, self.modifedbw, self.pdf)
    #     return out

    def test(self, width):
        width = self.mass * (width)
        self.x.setBins(1000, "cache")
        self.m0 = ROOT.RooRealVar("m0" + self.name + str(width), "m0", 0)
        self.m0.setConstant()
        self.w = ROOT.RooRealVar("width" + self.name + str(width), "width", width)
        self.w.setConstant()

        self.lnm0 = ROOT.RooRealVar("lnm0" + self.name + str(width), "lnm0", self.mass - width * 2, self.mass + width * 2)
        self.lnk = ROOT.RooRealVar("lnmk" + self.name + str(width), "lnmk", 1.001, 100)
        #self.lnk = ROOT.RooRealVar("lnmk" + self.name + str(width), "lnmk", 1.1)
        self.bw = ROOT.RooBreitWigner("bw" + self.name + str(width), "bw", self.x, self.m0, self.w)

        self.shift = ROOT.RooRealVar("shift" + self.name + str(width), "shift", self.mass)
        self.shift.setConstant()
        self.xshift = ROOT.RooFormulaVar("xshift" + self.name + str(width), "@0+@1", ROOT.RooArgSet(self.x, self.shift))
        self.ln = ROOT.RooLognormal("ln" + self.name + str(width), "ln", self.xshift, self.lnm0, self.lnk)

        self.bfgmean = ROOT.RooRealVar("bfgmean" + self.name + str(width), "bfgmean", 0)
        self.bfgmean.setConstant()
        self.bfgsig1 = ROOT.RooRealVar("bfgsig1" + self.name + str(width), "sig1", 10)
        self.bfgsig2 = ROOT.RooRealVar("bfgsig2" + self.name + str(width), "sig2", 10)
        self.bfg = ROOT.RooBifurGauss("bfg" + self.name + str(width), "bfg", self.x, self.bfgmean, self.bfgsig1, self.bfgsig2)
        self.modifedbw = ROOT.RooProdPdf("modifedbw" + self.name + str(width), "modifedbw", ROOT.RooArgSet(self.bfg, self.ln))
        out = ROOT.RooFFTConvPdf("out" + self.name + str(width), "out", self.x, self.pdf, self.modifedbw)

        return out

def plot_mass(mass, bins):
    mass = str(mass)
    NWidthobj = loadroot("ntuples/" + mass + "w0.root")

    for each_width in [1, 5, 10, 20]:
        LWidthobj = loadroot("ntuples/" + mass + "w" + str(each_width) + ".root")
        center, height, error = LWidthobj.gethist()
        # for i in range(len(height)):
        #     if height[i] > 0.00001:
        #         print(center[i], height[i], error[i])
        datalist = NWidthobj.getwidthdata(each_width/100)
        histplot_withsub_raw([center, datalist], bins, weights=[np.array(height), np.array([1]*len(datalist))], usererror=[error, None], 
                             labels = ["MC", "smearing"], removenorm=True, filename="com"+mass+"w"+str(each_width),
                             title2=r"$\mathit{\sqrt{s}=13\:TeV}$", title3="mass = "+mass+"GeV, W = " + str(each_width) + "%", central="smearing")


def main():
    global glabalx
    glabalx = ROOT.RooRealVar("x1","mA", -2000, 4000.)
    mass = str(500)
    NWidthobj = loadroot("ntuples/" + mass + "w0.root")
    outpdf = NWidthobj.test(0.1)
    LWidthobj = loadroot("ntuples/" + mass + "w10.root")
    datahist = LWidthobj.getDatahist()
    #outpdf.fitTo(datahist)
    xframe = glabalx.frame()
    #datahist.plotOn(xframe)
    outpdf.plotOn(xframe)
    xframe.Draw()
    # lnm0 = ROOT.RooRealVar("lnm0", "lnm0", 1000)
    # lnk = ROOT.RooRealVar("lnmk", "lnmk", 1.1)
    # # bw = ROOT.RooBreitWigner("bw", "bw", glabalx, lnm0, w)
    # ln = ROOT.RooLognormal("ln", "ln", glabalx, lnm0, lnk)
    # xframe = glabalx.frame()
    # ln.plotOn(xframe)
    # xframe.Draw()
    input()

    exit(1)
    # b = loadroot("ntuples300w5.root")
    # b.loadfile()
    # # b.getwidthdata(0.1)
    # print(b.gethist())
    bins = linspace(0, 1000, 20)
    plot_mass(300, bins)
    bins = linspace(0, 1000, 20)
    plot_mass(500, bins)
    bins = linspace(400, 1600, 20)
    plot_mass(1000, bins)
    bins = linspace(1000, 3000, 20)
    plot_mass(2000, bins)

if __name__ == '__main__':
    main()
