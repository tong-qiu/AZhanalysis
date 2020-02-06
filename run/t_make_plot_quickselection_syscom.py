import os

import uproot

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib import rc
import sys
import zlib

lib_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.append(lib_path)

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

def stack_cxaod(sample_directory, each_names, each_alias, each_color, branches_list_data, debug, cut, m_allsamples, matas=None):
    sample = load_CxAODs(sample_directory,each_names,branches_list_data, debug, 
                        colour=each_color,alias=each_alias,matanames=matas)
    if not sample:
        print("Warning: No "+each_alias+" samples found!")
    if cut and sample:
        sample.matacut(s_sr)
        sample.matacut(s_merged)
        sample.cut_parameter(cut_btagin_is, 2)
        sample.cut_parameter(cut_btagout_more, 1)
        #sample.cut(srcut)
        #sample.cut(cut_btag)
        #sample.cut(cut_electron)
        #sample.more()
        #sample.cut(crtopcut)
        #sample.cut(crmbbcut)

        m_allsamples.append(sample)
    if not cut:
        m_allsamples.append(sample)

    #print(each_alias)
    return 0

if __name__ == '__main__':
    debug = False
    cut = True
    sample_directory = ["../reduced_easytree/"]
    tag = "run2"

    t2 = r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$"
    if tag == "a":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$"
        sample_directory = ["../sample/CxAOD32_06" + tag + "/"]
    if tag == "d":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,43.6\:fb^{-1}}$"
        sample_directory = ["../sample/CxAOD32_06" + tag + "/"]
    if tag == "e":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,58.5\:fb^{-1}}$"
        sample_directory = ["../sample/CxAOD32_06" + tag + "/"]
    if tag == "run2":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,138.2\:fb^{-1}}$"
        sample_directory = ["../sample/reduced_easytree/CxAOD32_14a/", "../sample/reduced_easytree/CxAOD32_14d/", "../sample/reduced_easytree/CxAOD32_14e/"]
        #sample_directory = ["../phi/a/", "../phi/d/", "../phi/e/"]
        #sample_directory = ["../phi/a/"]
    mc_Wlvjet = ["WenuB_Sh221", "WenuC_Sh221", "WenuL_Sh221", "WmunuB_Sh221", "WmunuC_Sh221", "WmunuL_Sh221", "WtaunuB_Sh221", "WtaunuC_Sh221", "WtaunuL_Sh221"]
    mc_Zlljet1 = ["Zee_Sh221", "ZeeB_Sh221"]
    mc_Zlljet2 = ["ZeeC_Sh221", "ZeeL_Sh221"]
    mc_Zlljet3 = ["Zmumu_Sh221", "ZmumuB_Sh221"]
    mc_Zlljet4 = ["ZmumuC_Sh221", "ZmumuL_Sh221"]
    mc_Zlljet5 = ["Ztautau_Sh221", "ZtautauB_Sh221", "ZtautauC_Sh221", "ZtautauL_Sh221","Znunu_Sh221", "ZnunuB_Sh221", "ZnunuC_Sh221", "ZnunuL_Sh221"]
    mc_tt_bar = [ "ttbar_nonallhad_PwPy8", "ttbar_allhad_PwPy8", "ttbar_dilep_PwPy8"]#"ttbar_nonallhad_PwPy8", , "ttbar_allhad_PwPy8"]#"ttbar_nonallhad_PwPy8"]#, "ttbar_allhad_PwPy8"]
    file_name_array = [mc_tt_bar]
    alias = ["mc_tt_bar"]
    colors = ['royalblue']

    mc_Wlvjet_mg = ["WenuB_MGPy8", "WenuC_MGPy8", "WenuL_MGPy8", "WmunuB_MGPy8", "WmunuC_MGPy8", "WmunuL_MGPy8", "WtaunuB_MGPy8", "WtaunuC_MGPy8", "WtaunuL_MGPy8"]
    mc_Zlljet1_mg = ["ZeeB_MGPy8"]
    mc_Zlljet2_mg = ["ZeeC_MGPy8", "ZeeL_MGPy8"]
    mc_Zlljet3_mg = ["ZmumuB_MGPy8"]
    mc_Zlljet4_mg = ["ZmumuC_MGPy8", "ZmumuL_MGPy8"]


    branches_list_data = [b"EventWeight", b'nbTagsInFJ', b'nBTrkJets', b'nbTagsOutsideFJ', b'MCEventWeight_nominal', b'MCEventWeight_FSRUp', b'MCEventWeight_FSRDown']
    matas = ["Sample", "Description", "Regime"] #"Regime",
    branches_list_MC = branches_list_data
    bins = range(100,1400,50)
    bins = range(100,145,5)
    #bins = np.linspace(100,140,16)
    #all_sample = []
    processes = []
    manager = multiprocessing.Manager()
    all_sample = manager.list()
    for each_names, each_alias, each_color in zip(file_name_array,alias,colors):
        if "data" in each_alias:
            t = multiprocessing.Process(target=stack_cxaod, args=(sample_directory, each_names, each_alias, each_color, branches_list_data, debug, cut, all_sample, matas))
        else:
            t = multiprocessing.Process(target=stack_cxaod, args=(sample_directory, each_names, each_alias, each_color, branches_list_MC, debug, cut, all_sample, matas))
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
    print(all_sample_after)

    total = all_sample_after[0].weight.sum()
    upu = (all_sample_after[0].weight*all_sample_after[0].data[b'MCEventWeight_FSRUp']/all_sample_after[0].data[b'MCEventWeight_nominal']).sum()
    downu = (all_sample_after[0].weight*all_sample_after[0].data[b'MCEventWeight_FSRDown']/all_sample_after[0].data[b'MCEventWeight_nominal']).sum()
    print(total)
    print(upu)
    print(downu)
    print(upu/total, downu/total)