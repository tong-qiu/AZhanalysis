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
# lib_path = os.path.abspath(os.path.join(__file__, '..', '..'))
# sys.path.append(lib_path)

from package.events import *
from package.cut import *
from package.stackplot import *
from curveplot import *
from cutstring import *
import multiprocessing
from package.loadnormfactor import *

def poly(x, argv):
    s = 0
    for i, each in enumerate(argv):
        s += x**i * each
    return s

def rescale(inputfile, cfg="C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/confignormonly.cfg", 
            txt="C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/GlobalFit_fitres_unconditionnal_mu0_normonly.txt"):
    rescaledic = loadnorm(cfg, txt)
    inputfile = copy.deepcopy(inputfile)
    for i in range(len(inputfile)):
        for each_key in rescaledic.keys():
            if 'ALL' in rescaledic[each_key]:
                factor = rescaledic[each_key]['ALL'] + 1
                mask = inputfile[i].mata["Sample"] == zlib.adler32(each_key.encode())
                if True in mask:
                    inputfile[i].rescale(factor, mask)
    return inputfile

def slopecorrection(inputfile, csv="../run/output/slopefit/" + "pTV-mbbcut-"+str(1)+"tagpolyfitresult.csv"):
    inputfile = copy.deepcopy(inputfile)
    p1s = []
    p2s = []
    bottom = 0
    middle = 0
    top = 0
    with open(csv) as f:
        for each in f:
            each_array = each.split(',')
            if top == 0:
                bottom = float(each_array[0])
                middle = float(each_array[1])
                top = float(each_array[2])
            elif not p1s:
                p1s = each_array[0:-1]
            else:
                p2s = each_array[0:-1]
    for i in range(len(p1s)):
        p1s[i] = float(p1s[i])
        p2s[i] = float(p2s[i])

    for i in range(len(inputfile)):
        mask1 = inputfile[i].data[b'pTV']/1000. < middle
        mask2 = inputfile[i].data[b'pTV']/1000. >= middle
        mask2 = np.logical_and(inputfile[i].data[b'pTV']/1000. < top, mask2)
        print("before", inputfile[i].weight.sum())
        if True in mask1:
            inputfile[i].weight[mask1] = inputfile[i].weight[mask1] * (poly(inputfile[i].data[b'pTV'][mask1]/1000., p1s))
        if True in mask2:
            inputfile[i].weight[mask2] = inputfile[i].weight[mask2] * (poly(inputfile[i].data[b'pTV'][mask2]/1000., p2s))
    return inputfile




def _stack_cxaod(sample_directory, each_names, each_alias, each_color, branches_list_data, debug, selection, m_allsamples, matas=None):
    sample = load_CxAODs(sample_directory,each_names,branches_list_data, debug, 
                        colour=each_color,alias=each_alias,matanames=matas)
    if not sample:
        print("Warning: No "+each_alias+" samples found!")
    if sample:
        if "highmbbcr" in selection:
            sample.matacut(s_mbbcr)
            sample.cut(cut_highmbb)
        if "lowmbbcr" in selection:
            sample.matacut(s_mbbcr)
            sample.cut(cut_highmbb)
        if "mbbcr" in selection:
            sample.matacut(s_mbbcr)
        if "sr" in selection:
            sample.matacut(s_sr)
        if "topemucr" in selection:
            sample.matacut(s_topemucr)

        if "resolved" in selection:
            sample.matacut(s_resolved)
        if "merged" in selection:
            sample.matacut(s_merged)

        if "1tag" in selection:
            sample.cut_parameter(cut_btag_is, 1)
        if "2tag" in selection:
            sample.cut_parameter(cut_btag_is, 2)
        if "3ptag" in selection:
            sample.cut_parameter(cut_btag_more, 3)
        m_allsamples.append(sample)
    return 0

def get_sample(selection, debug=False, sample_list=None, tag="run2"):
    sample_directory = ["../CxAOD31_01a/"]

    t2 = r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$"
    if tag == "a":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$"
        sample_directory = ["../sample/" + tag + "/"]
    if tag == "d":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,43.6\:fb^{-1}}$"
        sample_directory = ["../sample/" + tag + "/"]
    if tag == "e":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,58.5\:fb^{-1}}$"
        sample_directory = ["../sample/" + tag + "/"]
    if tag == "run2":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,139\:fb^{-1}}$"
        sample_directory = ["../sample/a/", "../sample/d/", "../sample/e/"]

    mc_Wlvjet = ["Wenu_Sh221", "WenuB_Sh221", "WenuC_Sh221", "WenuL_Sh221", "Wmunu_Sh221", "WmunuB_Sh221", "WmunuC_Sh221", "WmunuL_Sh221", "Wtaunu_Sh221", "WtaunuB_Sh221", "WtaunuC_Sh221", "WtaunuL_Sh221"]
    mc_Zlljet1 = ["Zee_Sh221", "ZeeB_Sh221"]
    mc_Zlljet2 = ["ZeeC_Sh221", "ZeeL_Sh221"]
    mc_Zlljet3 = ["Zmumu_Sh221", "ZmumuB_Sh221"]
    mc_Zlljet4 = ["ZmumuC_Sh221", "ZmumuL_Sh221"]
    mc_Zlljet5 = ["Ztautau_Sh221", "ZtautauB_Sh221", "ZtautauC_Sh221", "ZtautauL_Sh221","Znunu_Sh221", "ZnunuB_Sh221", "ZnunuC_Sh221", "ZnunuL_Sh221"]
    mc_tt_bar = [ "ttbar_nonallhad_PwPy8", "ttbar_allhad_PwPy8", "ttbar_dilep_PwPy8"]#"ttbar_nonallhad_PwPy8", , "ttbar_allhad_PwPy8"]#"ttbar_nonallhad_PwPy8"]#, "ttbar_allhad_PwPy8"]
    mc_singletop = ["stops_PwPy8", "stopt_PwPy8", "stopWt_dilep_PwPy8"] # "stopWt_PwPy8", 
    mc_Diboson = ["WqqWlv_Sh221", "WqqZll_Sh221", "WqqZvv_Sh221", "ZqqZll_Sh221", "ZqqZvv_Sh221", "WlvZqq_Sh221", "ggZqqZll_Sh222", "ggWqqWlv_Sh222"]
    #sm_Higgs = ["qqWlvHbbJ_PwPy8MINLO", "qqZllHbbJ_PwPy8MINLO", "qqZvvHbbJ_PwPy8MINLO", "ggZllHbb_PwPy8", "ggZvvHbb_PwPy8", "ggHbb_PwPy8NNLOPS"] 
    sm_Higgs = ["bbHinc_aMCatNLOPy8", "ggHinc_PwPy8", "ggZllHbb_PwPy8","ggZllHcc_PwPy8","ggZvvHbb_PwPy8","ggZvvHcc_PwPy8","qqWlvHbbJ_PwPy8MINLO","qqWlvHccJ_PwPy8MINLO","qqZllHbbJ_PwPy8MINLO","qqZllHccJ_PwPy8MINLO","qqZvvHbbJ_PwPy8MINLO","qqZvvHccJ_PwPy8MINLO"]
    #other = ["ggZqqZll_Sh222", "ggWqqWlv_Sh222"]#,"ttV_aMCatNLOPy8","ggWqqWlv_Sh222","ggZqqZvv_Sh222","stoptZ_MGPy8"]#[ "ttV_aMCatNLOPy8"]#"VV_fulllep_Sh222",
    data = ["data16", "data15", "data17", "data18"]
    bbA300 = [ "bbA300"]
    ggA300 = [ "ggA300"]

    file_name_array = []
    alias = []
    colors = []
    if sample_list is None:
        file_name_array = [data, mc_Diboson, mc_tt_bar,  mc_singletop, mc_Zlljet1, mc_Zlljet2, mc_Zlljet3, mc_Zlljet4, mc_Zlljet5, mc_Wlvjet, sm_Higgs]#, bbA300, ggA300]#, sm_ggHiggs, sm_qqHiggs]
        alias = ["data", "Diboson", "ttbar", "singletop", "Zlljet", "Zlljet", "Zlljet", "Zlljet", "Zlljet", "Wlvjet", "smHiggs"]#,'bbA300', 'ggA300']#, "sm_ggHiggs", "sm_qqHiggs"]
        colors = [None,   'g',       'yellow',     'tab:orange',   'royalblue', 'royalblue', 'royalblue', 'royalblue', 'royalblue', 'm',    'teal']#, 'k', 'dimgrey']
    else:
        if "zlljet" in sample_list:
            file_name_array += [mc_Zlljet1, mc_Zlljet2, mc_Zlljet3, mc_Zlljet4, mc_Zlljet5]
            alias += ["Zlljet", "Zlljet", "Zlljet", "Zlljet", "Zlljet"]
            colors += ['royalblue', 'royalblue', 'royalblue', 'royalblue', 'royalblue']
        if "data" in sample_list:
            file_name_array += [data]
            alias += ["data"]
            colors += [None]
    branches_list_data = [b"mBBres", b"EventWeight", b"pTV", b'mVH', b'nTags', b'j1px', b'j1py', b'j2px', b'j2py', b'ptH', b'ptHcorr', b'ptL1', b'ptL2', b'METHT']
    matas = ["Sample", "Description", "Regime"]
    branches_list_MC = branches_list_data

    processes = []
    manager = multiprocessing.Manager()
    all_sample = manager.list()
    for each_names, each_alias, each_color in zip(file_name_array,alias,colors):
        if "data" in each_alias:
            t = multiprocessing.Process(target=_stack_cxaod, args=(sample_directory, each_names, each_alias, each_color, branches_list_data, debug, selection, all_sample, matas))
        else:
            t = multiprocessing.Process(target=_stack_cxaod, args=(sample_directory, each_names, each_alias, each_color, branches_list_MC, debug, selection, all_sample, matas))
        processes.append(t)
        t.start()

    i = 0
    for each_process, each_alias in zip(processes, alias):
        i += 1
        print(i," Waiting for " + each_alias + "...")
        each_process.join()
        print(i, each_alias + " finished.")
    print("All done.")
    all_sample_after = [each for each in all_sample]

    new_data_list = []
    for each in all_sample_after:
        ifpass = False
        for i, each_new in enumerate(new_data_list):
            if each_new.alias == each.alias:
                new_data_list[i] = each_new + each
                print("Info: Merge ", each_new.alias)
                ifpass = True
                break
        if not ifpass:
            new_data_list.append(each)

    return (new_data_list, t2)

if __name__ == "__main__":
    with open("../run/output/slopefit/" + "pTV-mbbcut-"+str(1)+"tagpolyfitresult.csv") as f:
        pass
    mysample, t2 = get_sample(["mbbcr", "resolved", "1tag"])
    mysample = rescale(mysample)
    mysample = slopecorrection(mysample)
    title3="mBBcr " + str(2) +" btags"
    bins = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300]
    print(bins)
    chi2, nod = stackplot(mysample,b'pTV',bins,1000.,
            xlabel=r"$p_{TV}[GeV]$", title3=title3, filename="test1", print_height=True,
            title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
    bins = range(0,2000,30)
    chi2, nod = stackplot(mysample,b'mVH',bins,1000.,
            xlabel=r"$p_{TV}[GeV]$", title3=title3, filename="testmvh1", print_height=True,
            title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)