from numpy.core.function_base import linspace
import ROOT
import numpy as np


x = ROOT.RooRealVar("x1","mA", -2000, 3000.)

path = "500w0.root"
ptl1 = []
ptl2 = []
ptb1 = []
ma = []
path = path
name = "500"
mass = 500


# load data and convert to pdf
hist = ROOT.TH1F("hist" + name, "", 1000, -2000, 4000.)
hist.Sumw2()
f = ROOT.TFile(path)
t1 = f.Get("data")
for entry in t1:
    if entry.ptl2/1000. < 20:
        continue
    if entry.ptl2/1000. < 27:
        continue
    if entry.ptb1/1000. < 45:
        continue
    hist.Fill(entry.mA/1000.)

datahist = ROOT.RooDataHist("datahist" + name, "datahist", ROOT.RooArgList(x), hist)
pdf = ROOT.RooHistPdf("pdf" + name, "your pdf", ROOT.RooArgSet(x), datahist)



# define RooBreitWigner pdf
width = mass * (0.2)
x.setBins(1000, "cache")
m0 = ROOT.RooRealVar("m0" + name + str(width), "m0", 0)
m0.setConstant()
w = ROOT.RooRealVar("width" + name + str(width), "width", width)
w.setConstant()

lnm0 = ROOT.RooRealVar("lnm0" + name + str(width), "lnm0", mass - width * 2, mass + width * 2)
lnk = ROOT.RooRealVar("lnmk" + name + str(width), "lnmk", 1.001, 100)
bw = ROOT.RooBreitWigner("bw" + name + str(width), "bw", x, m0, w)


# define RooLognormal pdf
shift = ROOT.RooRealVar("shift" + name + str(width), "shift", mass)
shift.setConstant()
alist = ROOT.RooArgList(x, shift)
xshift = ROOT.RooFormulaVar("xshift", "title", "@0+@1", alist)
ln = ROOT.RooLognormal("ln" + name + str(width), "ln", xshift, lnm0, lnk)

# calculate product of RooBreitWigner and RooLognormal
modifedbw = ROOT.RooProdPdf("modifedbw" + name + str(width), "modifedbw", ROOT.RooArgList(bw, ln))

# convolute pdf
out = ROOT.RooFFTConvPdf("out" + name + str(width), "out", x, pdf, modifedbw)


# make plot
c = ROOT.TCanvas("convolution", "convolution", 600, 600)
ROOT.gPad.SetLeftMargin(0.15)
xframe = x.frame(ROOT.RooFit.Title("test"))
xframe.GetYaxis().SetTitleOffset(1.4)
out.plotOn(xframe)
xframe.Draw()
c.SaveAs("test.png")
