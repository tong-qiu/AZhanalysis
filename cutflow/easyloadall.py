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
lib_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.append(lib_path)

from package.events import *
from package.cut import *
from package.stackplot import *
from curveplot import *
from cutstring import *
import multiprocessing
from package.loadnormfactor import *


def _stack_cxaod(sample_directory, each_names, each_alias, each_color, branches_list_data, debug, m_allsamples, matas=None):
    sample = load_CxAODs(sample_directory,each_names,branches_list_data, debug, 
                        colour=each_color,alias=each_alias,matanames=matas)
    if not sample:
        print("Warning: No "+each_alias+" samples found!")
    if sample:
        m_allsamples.append(sample)
    return 0

def get_allsample(debug=False, sample_list=None, tag="run2"):
    sample_directory = ["../CxAOD31_01a/"]

    if tag == "a":
        sample_directory = ["../sample/" + tag + "/"]
    if tag == "d":
        sample_directory = ["../sample/" + tag + "/"]
    if tag == "e":
        sample_directory = ["../sample/" + tag + "/"]
    if tag == "run2":
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
    dbl = [ "HVT", "bbA"]


    file_name_array = []
    alias = []
    colors = []
    # if sample_list is None:
    #     file_name_array = [mc_Diboson, mc_tt_bar,  mc_singletop, mc_Zlljet1, mc_Zlljet2, mc_Zlljet3, mc_Zlljet4, mc_Zlljet5, mc_Wlvjet, dbl]#, bbA300, ggA300]#, sm_ggHiggs, sm_qqHiggs]
    #     alias = [ "Diboson", "ttbar", "singletop", "Zlljet", "Zlljet", "Zlljet", "Zlljet", "Zlljet", "Wlvjet", "dbl"]#,'bbA300', 'ggA300']#, "sm_ggHiggs", "sm_qqHiggs"]
    #     colors = ['g',       'yellow',     'tab:orange',   'royalblue', 'royalblue', 'royalblue', 'royalblue', 'royalblue', 'm', 'dimgrey']#, 'k']
    # else:
    #     if "zlljet" in sample_list:
    #         file_name_array += [mc_Zlljet1, mc_Zlljet2, mc_Zlljet3, mc_Zlljet4, mc_Zlljet5]
    #         alias += ["Zlljet", "Zlljet", "Zlljet", "Zlljet", "Zlljet"]
    #         colors += ['royalblue', 'royalblue', 'royalblue', 'royalblue', 'royalblue']
    file_name_array = []
    branches_list_data = [b"EventWeight", b'mVH', b'nTags', b'nSigJet', b"passedTrigger", b'flavL1', b'flavL2', b'chargeL1', b'chargeL2', b"mLL", b"METHT", b'pTV', b'j1px', b'j1py', b'mBB', b'ptL1', b'ptL2', b'nFatJets', b'etaL1', b'MCChannelNumber']#b'nbTagsInFJ',,  b'nTrkjetsInFJ', b'nbTagsOutsideFJ'
    matas = ["Sample", "Description", "Regime"]
    branches_list_MC = branches_list_data

    file_name_array = [mc_tt_bar]
    alias = ["ttbar"]
    colors = ['yellow']
    processes=[]
    manager = multiprocessing.Manager()
    all_sample = manager.list()
    for each_names, each_alias, each_color in zip(file_name_array,alias,colors):
        t = multiprocessing.Process(target=_stack_cxaod, args=(sample_directory, each_names, each_alias, each_color, branches_list_MC, False, all_sample, matas))
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

    return new_data_list
