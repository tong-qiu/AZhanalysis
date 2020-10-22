import os

import uproot

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib import rc
import sys
import zlib
import re
import copy

lib_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.append(lib_path)
sys.path.append('../package')

from package.events import *
from mlcut import *
from package.stackplot import *
from curveplot import *
from cutstring import *
import multiprocessing

def loadpickle(path):
    with open("pickles/" + path, 'rb') as f:
        out = pickle.load(f)
    return out

def eventpop(inputs, alias):
    for i in range(len(inputs)):
        if inputs[i].alias == alias:
            inputs.pop(i)
        break

def eventkepp(inputs, alias):
    poplist = []
    for i in range(len(inputs)):
        if inputs[i].alias != alias:
            poplist.append(i)
    return [x for i, x in enumerate(inputs) if i not in poplist]

def setcolour(inputs, colour):
    for each in inputs:
        each.colour = colour

def Zsplit(sample):
    all_sample_after = copy.deepcopy(sample)
    poplist = []
    for i in range(len(all_sample_after)):
        if all_sample_after[i].alias == "Zlljet":
            poplist.append(i)
    temzjet = []
    for i in poplist:
        temzjet.append(all_sample_after[i])
    zhf, zlf = splitszjetsamples(temzjet)
    zlf.colour = "lightblue"
    all_sample_after = [each for i, each in enumerate(all_sample_after) if i not in poplist]
    all_sample_after.append(zhf)
    all_sample_after.append(zlf)
    return all_sample_after

def counttotal(samples):
    for each in samples:
        if "Zlljet" == each.alias:
            total = len(each.weight)
    return total

def flow(new, old, ntag):

    t2 = r"$\mathit{\sqrt{s}=13\:TeV,139\:fb^{-1}}$"
    title3 = "mBBcr " + str(ntag) +" btags"
    bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1150, 1350, 1550, 1800]
    chi2, nod = stackplot(Zsplit(new), b'mVHres', bins, 1000.,
        xlabel=r"$m_{VH}[GeV]$", title3=title3, filename="mVH-new"+str(ntag)+"tag", print_height=False,
        title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True)
    chi2, nod = stackplot(Zsplit(old), b'mVHres', bins, 1000.,
        xlabel=r"$m_{VH}[GeV]$", title3=title3, filename="mVH-old"+str(ntag)+"tag", print_height=False,
        title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True)

    bins = np.array([2, 3, 4, 5, 6, 7, 8, 9])-0.001
    chi2, nod = stackplot(Zsplit(new), b"nJets", bins, 1,
        xlabel=r"number of jets", title3=title3, filename="njet-new"+str(ntag)+"tag", print_height=False,
        title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True)
    chi2, nod = stackplot(Zsplit(old), b"nJets", bins, 1,
        xlabel=r"number of jets", title3=title3, filename="njet-old"+str(ntag)+"tag", print_height=False,
        title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True)

    bins = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4.5])-0.001
    chi2, nod = stackplot(Zsplit(new), b"dEtaBB", bins, 1,
        xlabel=r"$\Delta\eta_{BB}$", title3=title3, filename="etabb-new"+str(ntag)+"tag", print_height=False,
        title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True)
    chi2, nod = stackplot(Zsplit(old), b"dEtaBB", bins, 1,
        xlabel=r"$\Delta\eta_{BB}$", title3=title3, filename="etabb-old"+str(ntag)+"tag", print_height=False,
        title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True)

    bins = [150, 200, 250, 300, 350, 400, 500, 600, 800, 1800]
    chi2, nod = stackplot(Zsplit(new), b"ptH", bins, 1000.,
        xlabel=r"$pT_{BB}[GeV]$", title3=title3, filename="ptH-new"+str(ntag)+"tag", print_height=False,
        title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True)
    chi2, nod = stackplot(Zsplit(old), b"ptH", bins, 1000.,
        xlabel=r"$pT_{BB}[GeV]$", title3=title3, filename="ptH-old"+str(ntag)+"tag", print_height=False,
        title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True)

    data = copy.deepcopy(old)
    eventpop(data, "Zlljet")
    for each in data:
        if each.alias != "data":
            each.weight *= -1 

    new = eventkepp(new, "Zlljet")
    setcolour(new, "r")
    old = eventkepp(old, "Zlljet")
    setcolour(old, "b")
    setcolour(data, "k")

    bins = [150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000, 1500, 1800]
    histplot_withsub([new, old, data], b'mVHres', bins, labels=["Sherpa v2.2.10", "Sherpa v2.2.1", "data"], scales=1000., removenorm=None, filename="mVH"+str(ntag)+"tag" , central="data", title3=title3, title2=t2, do_errorbar=True)
    histplot_withsub([new, old, data], b'mVHres', bins, labels=["Sherpa v2.2.10", "Sherpa v2.2.1", "data"], scales=1000., removenorm="data", filename="mVH"+str(ntag)+"tag_density", central="data", title3=title3, title2=t2, do_errorbar=True)
    bins = [150, 200, 250, 300, 350, 400, 500, 600, 800, 1800]
    histplot_withsub([new, old, data], b"ptH", bins, labels=["Sherpa v2.2.10", "Sherpa v2.2.1", "data"], scales=1000., xlabel=r"$pT_{BB}[GeV]$", removenorm=None, filename="ptbb"+str(ntag)+"tag" , central="data", title3=title3, title2=t2, do_errorbar=True)
    histplot_withsub([new, old, data], b"ptH", bins, labels=["Sherpa v2.2.10", "Sherpa v2.2.1", "data"], scales=1000., xlabel=r"$pT_{BB}[GeV]$",  removenorm="data", filename="ptbb"+str(ntag)+"tag_density", central="data", title3=title3, title2=t2, do_errorbar=True)
    bins = [150, 200, 250, 300, 350, 400, 500, 600, 800, 1800]
    histplot_withsub([new, old, data], b"pTV", bins, labels=["Sherpa v2.2.10", "Sherpa v2.2.1", "data"], scales=1000., xlabel=r"$pT_{V}[GeV]$", removenorm=None, filename="ptV"+str(ntag)+"tag" , central="data", title3=title3, title2=t2, do_errorbar=True)
    histplot_withsub([new, old, data], b"pTV", bins, labels=["Sherpa v2.2.10", "Sherpa v2.2.1", "data"], scales=1000., xlabel=r"$pT_{V}[GeV]$",  removenorm="data", filename="ptV"+str(ntag)+"tag_density", central="data", title3=title3, title2=t2, do_errorbar=True)
    bins = np.array([2, 3, 4, 5, 6, 7, 8, 9])-0.001
    histplot_withsub([new, old, data], b"nJets", bins, labels=["Sherpa v2.2.10", "Sherpa v2.2.1", "data"], scales=1, xlabel=r"number of jets", removenorm=None, filename="njet"+str(ntag)+"tag" , central="data", title3=title3, title2=t2, do_errorbar=True)
    histplot_withsub([new, old, data], b"nJets", bins, labels=["Sherpa v2.2.10", "Sherpa v2.2.1", "data"], scales=1, xlabel=r"number of jets",  removenorm="data", filename="njet"+str(ntag)+"tag_density", central="data", title3=title3, title2=t2, do_errorbar=True)
    bins = np.array([0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4.5])-0.001
    histplot_withsub([new, old, data], b"dEtaBB", bins, labels=["Sherpa v2.2.10", "Sherpa v2.2.1", "data"], scales=1, xlabel=r"$\Delta\eta_{BB}$", removenorm=None, filename="eta"+str(ntag)+"tag" , central="data", title3=title3, title2=t2, do_errorbar=True)
    histplot_withsub([new, old, data], b"dEtaBB", bins, labels=["Sherpa v2.2.10", "Sherpa v2.2.1", "data"], scales=1, xlabel=r"$\Delta\eta_{BB}$",  removenorm="data", filename="eta"+str(ntag)+"tag_density", central="data", title3=title3, title2=t2, do_errorbar=True)

sample1 = loadpickle("background-1.pickle")
samplenominal1 = loadpickle("background-nominal-1.pickle")
print(counttotal(sample1), counttotal(samplenominal1))
flow(sample1, samplenominal1, 1)

sample1 = loadpickle("background-2.pickle")
samplenominal1 = loadpickle("background-nominal-2.pickle")
print(counttotal(sample1), counttotal(samplenominal1))
flow(sample1, samplenominal1, 2)

sample1 = loadpickle("background-3.pickle")
samplenominal1 = loadpickle("background-nominal-3.pickle")
print(counttotal(sample1), counttotal(samplenominal1))
flow(sample1, samplenominal1, "3+")
print("end")