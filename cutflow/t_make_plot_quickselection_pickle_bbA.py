import os

import uproot

import json
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib import rc
import sys
import zlib
import copy
lib_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.append(lib_path)

from easyloadall import *
from package.events import *
from package.cut import *
from package.stackplot import *
from curveplot import *
from cutstring import *
import multiprocessing
from package.loadnormfactor import *

def cut_ptl1(data):
    mask = data[b'ptL1'] > 27000
    return mask

def cut_ptl2tresolved(data):
    mask = data[b'ptL2'] > 20000
    return mask

def cut_ptl2tmerged(data):
    mask = data[b'ptL2'] > 25000
    return mask

def hvt300(data):
    mask = data[b'MCChannelNumber'] == 343509
    return mask
def hvt700(data):
    mask = data[b"MCChannelNumber"] == 343516
    return mask
def hvt2000(data):
    mask = data[b"MCChannelNumber"] == 302406
    return mask
def hvt5000(data):
    mask = data[b"MCChannelNumber"] == 302415
    return mask

def passtrigger(data):
    mask = data[b"passedTrigger"] == 1
    return mask

def leptonsf(data):
    mask = data[b'flavL1'] == data[b'flavL2']
    return mask

def osmuon(data):
    mask = np.logical_or(data[b'flavL1'] == 1, data[b'chargeL1'] != data[b'chargeL2'])
    return mask

def mll(data):
    lower_limit = [max(40, 87 - 0.03 * each) for each in (data[b"mVH"]/1000.)]
    higher_limit = 97 + 0.013 * data[b"mVH"]/1000.
    mask = lower_limit < data[b"mLL"]/1000.
    mask = np.logical_and(data[b"mLL"]/1000. < higher_limit, mask)
    return mask

@np.vectorize
def ptll1_cut(mvh):
    if mvh < 320.:
        return 0.
    return 20. + 9. * pow(mvh - 320., 0.5)

def ptvres(data):
    mask = data[b'pTV']/1000. > ptll1_cut(data[b'mVHres']/1000.)
    return mask

def ptvmer(data):
    mask = data[b'pTV']/1000. > ptll1_cut(data[b'mVHmerg']/1000.)
    return mask

def metht(data):
    mask = data[b"METHT"]/(1000.**0.5) < 1.15 + 8 * (10**(-3))*data[b"mVH"]/1000.
    return mask

def muoneta(data):
    mask = np.logical_or(data[b'flavL1'] == 1, abs(data[b'etaL1']) < 2.5)
    return mask

def resolvedleptonpt(data):
    mask = data[b'ptL1']/1000. > 27
    mask = np.logical_and(data[b'ptL2']/1000. > 20, mask) # essential
    return mask

def mergedleptonpt(data):
    mask = data[b'ptL1']/1000. > 27
    mask = np.logical_and(data[b'ptL2']/1000. > 25, mask) # essential
    return mask

def resolvedb1pt(data):
    mask = ((data[b'j1px']/1000.)**2 + (data[b'j1py']/1000.)**2)**0.5 > 45
    return mask

def mergedb1pt(data):
    mask = ((data[b'j1px']/1000.)**2 + (data[b'j1py']/1000.)**2)**0.5 > 250
    return mask

def resolvedmbb(data):
    mask = data[b'mBB']/1000. > 100
    mask = np.logical_and(data[b'mBB']/1000. < 145, mask)
    return mask

def mergedmbb(data):
    mask = data[b'mBB']/1000. > 75
    mask = np.logical_and(data[b'mBB']/1000. < 145, mask)
    return mask

def resolvednjet(data):
    mask = data[b"nSigJet"] >= 2
    return mask

def mergednjet(data):
    mask = data[b"nFatJets"] >= 1
    return mask

def ntrackjet(data):
    mask = data[b'nTrkjetsInFJ'] >= 2
    return mask

def resolved1btag(data):
    mask = data[b'nTags'] >= 1
    return mask

def resolved3pbtag(data):
    mask = data[b'nTags'] >= 3
    return mask

def resolved4pbtag(data):
    mask = data[b'nTags'] >= 4
    return mask

def merged1btag(data):
    mask = data[b'nbTagsInFJ'] >= 1
    return mask

def calculatecutcommonflow(samples, domerged):
    resolveoutput = []
    mergedoutput = []
    resolveoutputlen = []
    mergedoutputlen = []

    sum_tem = np.sum(samples.weight)
    resolveoutput.append(["all CxAODs", sum_tem])
    mergedoutput.append(["all CxAODs", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append(["all CxAODs", sum_temlen])
    mergedoutputlen.append(["all CxAODs", sum_temlen])

    samples.cut(passtrigger)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append(["pass triggers", sum_tem])
    mergedoutput.append(["pass triggers", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append(["pass triggers", sum_temlen])
    mergedoutputlen.append(["pass triggers", sum_temlen])

    samples.cut(leptonsf)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append(["lepton same flavor", sum_tem])
    mergedoutput.append(["lepton same flavor", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append(["lepton same flavor", sum_temlen])
    mergedoutputlen.append(["lepton same flavor", sum_temlen])

    samples.cut(cut_ptl1)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append(["ptl1 > 27 GeV", sum_tem])
    mergedoutput.append(["ptl1 > 27 GeV", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append(["ptl1 > 27 GeV", sum_temlen])
    mergedoutputlen.append(["ptl1 > 27 GeV", sum_temlen])

    samples.cut(osmuon)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append(["muon opposite charge", sum_tem])
    mergedoutput.append(["muon opposite charge", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append(["muon opposite charge", sum_temlen])
    mergedoutputlen.append(["muon opposite charge", sum_temlen])

    samples.cut(muoneta)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append(["muon eta", sum_tem])
    mergedoutput.append(["muon eta", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append(["muon eta", sum_temlen])
    mergedoutputlen.append(["muon eta", sum_temlen])

    samples.cut(mll)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append(["mll", sum_tem])
    mergedoutput.append(["mll", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append(["mll", sum_temlen])
    mergedoutputlen.append(["mll", sum_temlen])

    samples.cut(metht)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append(["METHT", sum_tem])
    mergedoutput.append(["METHT", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append(["METHT", sum_temlen])
    mergedoutputlen.append(["METHT", sum_temlen])


    samplesmerged = copy.deepcopy(samples)

    samples.cut(cut_ptl2tresolved)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append(["ptl2 > 20 GeV", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append(["ptl2 > 20 GeV", sum_temlen])

    samples.cut(ptvres)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append(["ptV", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append(["ptV", sum_temlen])

    samples.cut(resolvedb1pt)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append(["pTB1", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append(["pTB1", sum_temlen])

    samples.cut(resolvedmbb)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append(["mBB", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append(["mBB", sum_temlen])

    samples.cut(resolvednjet)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append([">= 2 signal jets", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append([">= 2 signal jets", sum_temlen])

    samples.cut(resolved1btag)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append([">= 1 b-tagged", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append([">= 1 b-tagged", sum_temlen])

    samples.cut(resolved3pbtag)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append([">= 3 b-tagged", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append([">= 3 b-tagged", sum_temlen])

    samples.cut(resolved4pbtag)
    sum_tem = np.sum(samples.weight)
    resolveoutput.append([">= 4 b-tagged", sum_tem])
    sum_temlen = len(samples.weight)
    resolveoutputlen.append([">= 4 b-tagged", sum_temlen])

    if not domerged:
        return [{"resolved": resolveoutput}, {"resolved": resolveoutputlen}]

    samplesmerged.cut(cut_ptl2tmerged)
    sum_tem = np.sum(samplesmerged.weight)
    mergedoutput.append(["ptl2 > 25 GeV", sum_tem])
    sum_temlen = len(samplesmerged.weight)
    mergedoutputlen.append(["ptl2 > 25 GeV", sum_temlen])

    samplesmerged.cut(ptvmer)
    sum_tem = np.sum(samplesmerged.weight)
    mergedoutput.append(["ptV", sum_tem])
    sum_temlen = len(samplesmerged.weight)
    mergedoutputlen.append(["ptV", sum_temlen])

    samplesmerged.cut(mergedb1pt)
    sum_tem = np.sum(samplesmerged.weight)
    mergedoutput.append(["pTB1", sum_tem])
    sum_temlen = len(samplesmerged.weight)
    mergedoutputlen.append(["pTB1", sum_temlen])

    samplesmerged.cut(mergedmbb)
    sum_tem = np.sum(samplesmerged.weight)
    mergedoutput.append(["mBB", sum_tem])
    sum_temlen = len(samplesmerged.weight)
    mergedoutputlen.append(["mBB", sum_temlen])

    samplesmerged.cut(mergednjet)
    sum_tem = np.sum(samplesmerged.weight)
    mergedoutput.append([">= 1 fat jets", sum_tem])
    sum_temlen = len(samplesmerged.weight)
    mergedoutputlen.append([">= 1 fat jets", sum_temlen])

    samplesmerged.cut(ntrackjet)
    sum_tem = np.sum(samplesmerged.weight)
    mergedoutput.append([">= 2 track jets", sum_tem])
    sum_temlen = len(samplesmerged.weight)
    mergedoutputlen.append([">= 2 track jets", sum_temlen])

    samplesmerged.cut(merged1btag)
    sum_tem = np.sum(samplesmerged.weight)
    mergedoutput.append([">= 1 b-tagged", sum_tem])
    sum_temlen = len(samplesmerged.weight)
    mergedoutputlen.append([">= 1 b-tagged", sum_temlen])

    samplesmerged.matacut(s_merged)
    sum_tem = np.sum(samplesmerged.weight)
    mergedoutput.append(["priority resolved", sum_tem])
    sum_temlen = len(samplesmerged.weight)
    mergedoutputlen.append(["priority resolved", sum_temlen])

    return [{"resolved": resolveoutput, "merged": mergedoutput}, {"resolved": resolveoutputlen, "merged": mergedoutputlen}]
def _stack_cxaod(sample_directory, each_names, each_alias, each_color, branches_list_data, debug, m_allsamples, matas=None):
    sample = load_CxAODs(sample_directory,each_names,branches_list_data, debug, 
                        colour=each_color,alias=each_alias,matanames=matas)
    if not sample:
        print("Warning: No "+each_alias+" samples found!")
    if sample:
        m_allsamples.append(sample)
    return 0

def pickleit(obj, path):
    outfile = open(path, 'wb')
    pickle.dump(obj, outfile)
    outfile.close()

def unpickleit(path):
    infile = open(path, 'rb')
    output = pickle.load(infile)
    infile.close()
    return output

if __name__ == "__main__":
    sample_directory = ["../sample/a/", "../sample/d/", "../sample/e/"]
    branches_list_data = [b"EventWeight", b'mVHres', b'mVHmerg', b'mVH', b'nTags', b'nSigJet', b"passedTrigger", b'flavL1', b'flavL2', b'chargeL1', b'chargeL2', b"mLL", b"METHT", b'pTV', b'j1px', b'j1py', b'mBB', b'ptL1', b'ptL2', b'nFatJets', b'etaL1', b'MCChannelNumber', b'nbTagsInFJ', b'nTrkjetsInFJ', b'nbTagsOutsideFJ']
    mc_Wlvjet = ["Wenu_Sh221", "WenuB_Sh221", "WenuC_Sh221", "WenuL_Sh221", "Wmunu_Sh221", "WmunuB_Sh221", "WmunuC_Sh221", "WmunuL_Sh221", "Wtaunu_Sh221", "WtaunuB_Sh221", "WtaunuC_Sh221", "WtaunuL_Sh221"]
    mc_Zlljet1 = ["Zee_Sh221", "ZeeB_Sh221"]
    mc_Zlljet2 = ["ZeeC_Sh221", "ZeeL_Sh221"]
    mc_Zlljet3 = ["Zmumu_Sh221", "ZmumuB_Sh221"]
    mc_Zlljet4 = ["ZmumuC_Sh221", "ZmumuL_Sh221"]
    mc_Zlljet5 = ["Ztautau_Sh221", "ZtautauB_Sh221", "ZtautauC_Sh221", "ZtautauL_Sh221","Znunu_Sh221", "ZnunuB_Sh221", "ZnunuC_Sh221", "ZnunuL_Sh221"]
    mc_Zlljet = ["Znunu_Sh221", "ZnunuB_Sh221", "ZnunuC_Sh221", "ZnunuL_Sh221"]
    mc_Zlljet = ["Zee_Sh221", "ZeeB_Sh221", "ZeeC_Sh221", "ZeeL_Sh221", "Zmumu_Sh221", "ZmumuB_Sh221", "ZmumuC_Sh221", "ZmumuL_Sh221", "Ztautau_Sh221", "ZtautauB_Sh221", "ZtautauC_Sh221", "ZtautauL_Sh221"]
    mc_tt_bar = [ "ttbar_nonallhad_PwPy8", "ttbar_allhad_PwPy8", "ttbar_dilep_PwPy8"]#"ttbar_nonallhad_PwPy8", , "ttbar_allhad_PwPy8"]#"ttbar_nonallhad_PwPy8"]#, "ttbar_allhad_PwPy8"]
    mc_singletop = ["stops_PwPy8", "stopt_PwPy8", "stopWt_dilep_PwPy8"] # "stopWt_PwPy8", 
    mc_Diboson = ["WqqWlv_Sh221", "WqqZll_Sh221", "WqqZvv_Sh221", "ZqqZll_Sh221", "ZqqZvv_Sh221", "WlvZqq_Sh221", "ggZqqZll_Sh222", "ggWqqWlv_Sh222"]
    #sm_Higgs = ["qqWlvHbbJ_PwPy8MINLO", "qqZllHbbJ_PwPy8MINLO", "qqZvvHbbJ_PwPy8MINLO", "ggZllHbb_PwPy8", "ggZvvHbb_PwPy8", "ggHbb_PwPy8NNLOPS"] 
    sm_Higgs = ["bbHinc_aMCatNLOPy8", "ggHinc_PwPy8", "ggZllHbb_PwPy8","ggZllHcc_PwPy8","ggZvvHbb_PwPy8","ggZvvHcc_PwPy8","qqWlvHbbJ_PwPy8MINLO","qqWlvHccJ_PwPy8MINLO","qqZllHbbJ_PwPy8MINLO","qqZllHccJ_PwPy8MINLO","qqZvvHbbJ_PwPy8MINLO","qqZvvHccJ_PwPy8MINLO"]
    #other = ["ggZqqZll_Sh222", "ggWqqWlv_Sh222"]#,"ttV_aMCatNLOPy8","ggWqqWlv_Sh222","ggZqqZvv_Sh222","stoptZ_MGPy8"]#[ "ttV_aMCatNLOPy8"]#"VV_fulllep_Sh222",
    data = ["data16", "data15", "data17", "data18"]
    dbls = [ "HVT", "bbA"]
    matas = ["Regime"]
    file_name_array = [mc_tt_bar]
    alias = ["ttbar"]
    colors = ['yellow']

    loadit = False
    output = {}
    outputlen = {}
    ttbar = []
    zjet = []
    dbl = []
    stop = []
    diboson = []
    if loadit:
        # _stack_cxaod(sample_directory, mc_Diboson, "diboson", 'yellow', branches_list_data, False, diboson, matas)
        # pickleit(diboson, "diboson")
        # _stack_cxaod(sample_directory, mc_singletop, "stop", 'yellow', branches_list_data, False, stop, matas)
        # pickleit(stop, "stop")
        _stack_cxaod(sample_directory, dbls, "dbl", 'yellow', branches_list_data, False, dbl, matas)
        pickleit(dbl, "dbl")
        _stack_cxaod(sample_directory, mc_tt_bar, "ttbar", 'yellow', branches_list_data, False, ttbar, matas)
        pickleit(ttbar, "ttbar")
        _stack_cxaod(sample_directory, mc_Zlljet, "zlljet", 'blue', branches_list_data, False, zjet, matas)
        pickleit(zjet, "zlljet")
        exit(1)
    # ttbar = unpickleit("ttbar")
    # outputtem = calculatecutcommonflow(ttbar[0], True)
    # output["ttbar"] = outputtem[0]
    # outputlen["ttbar"] = outputtem[1]

    # zjet = unpickleit("zlljet")
    # outputtem = calculatecutcommonflow(zjet[0], True)
    # output["Z+jets"] = outputtem[0]
    # outputlen["Z+jets"] = outputtem[1]

    # diboson = unpickleit("diboson")
    # outputtem = calculatecutcommonflow(diboson[0], False)
    # output["diboson"] = outputtem[0]
    # outputlen["diboson"] = outputtem[1]

    # stop = unpickleit("stop")
    # outputtem = calculatecutcommonflow(stop[0], False)
    # output["stop"] = outputtem[0]
    # outputlen["stop"] = outputtem[1]

    dbl_tem = unpickleit("dbl")
    dbl_tem1 = copy.deepcopy(dbl_tem)
    dbl_tem1[0].cut(hvt300)
    outputtem = calculatecutcommonflow(dbl_tem1[0], True)
    output["bbA400"] = outputtem[0]
    outputlen["bbA400"] = outputtem[1]

    dbl_tem1 = copy.deepcopy(dbl_tem)
    dbl_tem1[0].cut(hvt700)
    outputtem = calculatecutcommonflow(dbl_tem1[0], True)
    output["bbA600"] = outputtem[0]
    outputlen["bbA600"] = outputtem[1]

    # dbl_tem1 = copy.deepcopy(dbl_tem)
    # dbl_tem1[0].cut(hvt2000)
    # outputtem = calculatecutcommonflow(dbl_tem1[0], True)
    # output["HVT2000"] = outputtem[0]
    # outputlen["HVT2000"] = outputtem[1]


    # dbl_tem1 = copy.deepcopy(dbl_tem)
    # dbl_tem1[0].cut(hvt5000)
    # outputtem = calculatecutcommonflow(dbl_tem1[0], True)
    # output["HVT5000"] = outputtem[0]
    # outputlen["HVT5000"] = outputtem[1]

    with open("cutflow.json", "w") as f:
        f.write(json.dumps(output))

    # with open("cutflowlen.json", "w") as f:
    #     f.write(json.dumps(outputlen))
