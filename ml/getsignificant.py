import os

import uproot

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib import rc
import sys
import zlib
import copy
import re
lib_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.append(lib_path)

from package.events import *
from mlcut import *
from package.stackplot import *
from curveplot import *
from cutstring import *
import multiprocessing

def get_signalid(samplenameininfo):
    masses = []
    ids = []
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
            if samplenameininfo == "ggA" or samplenameininfo == "bbA":
                for each in sample[2].split("_"):
                        if samplenameininfo in each:
                            masses.append(int(each.replace(samplenameininfo, "")))
                            break
            elif samplenameininfo == "HVT":
                masses.append(int(re.findall("\d+", sample[2].split("_")[-1])[0]))
        output = {}
    for eachmass, eachid in zip(masses, ids):
        output[eachid] = eachmass
    return output

def splitesamples(eventsobj, dicobj):
    output = []
    for eachkey in dicobj.keys():
        eventobj_tem = copy.deepcopy(eventsobj)
        eventobj_tem.cut(cut_mcid(eachkey))
        eventobj_tem.alias = dicobj[eachkey]
        output.append((dicobj[eachkey], eventobj_tem))
    return output

def pickleit(obj, path):
    outfile = open(path, 'wb')
    pickle.dump(obj, outfile)
    outfile.close()

def unpickleit(path):
    infile = open(path, 'rb')
    output = pickle.load(infile)
    infile.close()
    return output

def loadandcalculatesig(inputfile, binning):
    all_sample, title = unpickleit(inputfile)
    i_ggA = 0
    # pop sig
    for i in range(len(all_sample)):
        if all_sample[i].alias ==  "ggA":
            i_ggA = i
            break
    signal_sample = all_sample[i]
    all_sample.pop(i)

    # pop data
    for i in range(len(all_sample)):
        if "data" in all_sample[i].alias:
            i_ggA = i
            break
    all_sample.pop(i)

    allbkg = None
    for each in all_sample:
        if allbkg is None:
            allbkg = each
        else:
            allbkg = allbkg + each

    test = get_signalid("ggA")
    signals = splitesamples(signal_sample, test)
    output = []
    for each in signals:
        mass = each[0]
        sig = significant(allbkg, each[1], b"mVH", binning, scale=1000, logsig=True)
        output.append((mass, sig))
    output.sort(key=lambda a: a[0])
    return output


if __name__ == '__main__':
    binning1tagresolvedsr = [0.0, 190.0, 210.0, 230.0, 250.0, 270.0, 290.0, 310.0, 330.0, 350.0, 370.0, 390.0, 410.0, 430.0, 450.0, 470.0, 490.0,
    510.0, 530.0, 550.0, 570.0, 590.0, 610.0, 640.0, 670.0, 700.0, 730.0, 760.0, 790.0, 820.0, 880.0, 940.0, 1000.0, 1060.0,
    1120.0, 1180.0, 1240.0, 1300.0, 1360.0, 1420.0, 1480.0, 1540.0, 1600.0, 2000.0, 2400.0, 2800.0, 3200.0, 3600.0, 4320.0, 9000.0]
    print(loadandcalculatesig("../run/pickle/sr1tag.pickle", binning1tagresolvedsr))

    binning2tagresolvedsr = [0.0, 200.0, 220.0, 240.0, 260.0, 280.0, 300.0, 320.0, 340.0, 360.0, 380.0, 400.0, 420.0, 
    440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 630.0, 660.0, 690.0, 720.0, 750.0, 
    780.0, 810.0, 850.0, 910.0, 970.0, 1030.0, 1090.0, 1150.0, 1210.0, 1270.0, 1330.0, 1390.0, 1450.0, 
    1510.0, 1570.0, 1900.0, 2300.0, 2700.0, 3100.0, 9000.0]
    print(loadandcalculatesig("../run/pickle/sr2tag.pickle", binning2tagresolvedsr))
