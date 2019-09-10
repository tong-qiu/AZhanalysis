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
        sample.matacut(s_mbbcr)
        #sample.cut(cut_lowmbb)
        sample.cut(cut_highmbb)
        #sample.matacut(s_resolved)
        sample.cut_parameter(cut_btag_is, 2)
        #sample.cut(srcut)
        #sample.cut(cut_btag)
        #sample.cut(cut_muon)
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
    sample_directory = ["../CxAOD31_01a/"]
    tag = "run2"
    rescale = True
    slopecorrection = False

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
        sample_directory = ["../sample/CxAOD32_06a/", "../sample/CxAOD32_06d/", "../sample/CxAOD32_06e/"]
        #sample_directory = ["../sample/CxAOD32_06a/", "../sample/CxAOD32_06d/", "../sample/CxAOD32_06e/"]
        #sample_directory = ["../phi/a/", "../phi/d/", "../phi/e/"]
        #sample_directory = ["../phi/a/"]
    mc_Wlvjet = ["Wenu_Sh221", "WenuB_Sh221", "WenuC_Sh221", "WenuL_Sh221", "Wmunu_Sh221", "WmunuB_Sh221", "WmunuC_Sh221", "WmunuL_Sh221", "Wtaunu_Sh221", "WtaunuB_Sh221", "WtaunuC_Sh221", "WtaunuL_Sh221"]
    #mc_Zlljet = ["Zee_Sh221", "ZeeB_Sh221", "ZeeC_Sh221", "ZeeL_Sh221", "Zmumu_Sh221", "ZmumuB_Sh221", "ZmumuC_Sh221", "ZmumuL_Sh221", "Ztautau_Sh221", "ZtautauB_Sh221", "ZtautauC_Sh221", "ZtautauL_Sh221","Znunu_Sh221", "ZnunuB_Sh221", "ZnunuC_Sh221", "ZnunuL_Sh221"]
    # mc_Zlljet1 = ["ZeeB_MGPy8", "ZeeC_MGPy8"]
    # mc_Zlljet2 = ["ZeeL_MGPy8", "ZmumuB_MGPy8"]
    # mc_Zlljet3 = ["Zmumu_Sh221", "ZmumuB_MGPy8"]
    # mc_Zlljet4 = ["ZmumuC_MGPy8", "ZmumuL_MGPy8"]

    mc_Zlljet1 = ["Zee_Sh221", "ZeeB_Sh221"]
    mc_Zlljet2 = ["ZeeC_Sh221", "ZeeL_Sh221"]
    mc_Zlljet3 = ["Zmumu_Sh221", "ZmumuB_Sh221"]
    mc_Zlljet4 = ["ZmumuC_Sh221", "ZmumuL_Sh221"]
    mc_Zlljet5 = ["Ztautau_Sh221", "ZtautauB_Sh221", "ZtautauC_Sh221", "ZtautauL_Sh221","Znunu_Sh221", "ZnunuB_Sh221", "ZnunuC_Sh221", "ZnunuL_Sh221"]
    mc_tt_bar = [ "ttbar_nonallhad_PwPy8", "ttbar_allhad_PwPy8"]#"ttbar_nonallhad_PwPy8", , "ttbar_allhad_PwPy8"]#"ttbar_nonallhad_PwPy8"]#, "ttbar_allhad_PwPy8"]
    mc_singletop = ["stops_PwPy8", "stopt_PwPy8", "stopWt_PwPy8"]
    mc_Diboson = ["WqqWlv_Sh221", "WqqZll_Sh221", "WqqZvv_Sh221", "ZqqZll_Sh221", "ZqqZvv_Sh221", "WlvZqq_Sh221", "ggZqqZll_Sh222", "ggWqqWlv_Sh222"]
    #sm_Higgs = ["qqWlvHbbJ_PwPy8MINLO", "qqZllHbbJ_PwPy8MINLO", "qqZvvHbbJ_PwPy8MINLO", "ggZllHbb_PwPy8", "ggZvvHbb_PwPy8", "ggHbb_PwPy8NNLOPS"] 
    sm_Higgs = ["bbHinc_aMCatNLOPy8", "ggHinc_PwPy8", "ggZllHbb_PwPy8","ggZllHcc_PwPy8","ggZvvHbb_PwPy8","ggZvvHcc_PwPy8","qqWlvHbbJ_PwPy8MINLO","qqWlvHccJ_PwPy8MINLO","qqZllHbbJ_PwPy8MINLO","qqZllHccJ_PwPy8MINLO","qqZvvHbbJ_PwPy8MINLO","qqZvvHccJ_PwPy8MINLO"]
    #other = ["ggZqqZll_Sh222", "ggWqqWlv_Sh222"]#,"ttV_aMCatNLOPy8","ggWqqWlv_Sh222","ggZqqZvv_Sh222","stoptZ_MGPy8"]#[ "ttV_aMCatNLOPy8"]#"VV_fulllep_Sh222",
    data = ["data16", "data15", "data17", "data18"]
    bbA300 = [ "bbA300"]
    ggA300 = [ "ggA300"]
    file_name_array = [data, mc_Diboson, mc_tt_bar,  mc_singletop, mc_Zlljet1, mc_Zlljet2, mc_Zlljet3, mc_Zlljet4, mc_Zlljet5, mc_Wlvjet, sm_Higgs]#, bbA300, ggA300]#, sm_ggHiggs, sm_qqHiggs]
    alias = ["data", "Diboson", "ttbar", "singletop", "Zlljet", "Zlljet", "Zlljet", "Zlljet", "Zlljet", "Wlvjet", "smHiggs"]#,'bbA300', 'ggA300']#, "sm_ggHiggs", "sm_qqHiggs"]
    colors = [None,   'g',       'yellow',     'tab:orange',   'royalblue', 'royalblue', 'royalblue', 'royalblue', 'royalblue', 'm',    'teal']#, 'k', 'dimgrey']


    branches_list_data = [b"mBBres", b"EventWeight", b"pTV", b'mBB', b'mVH', b'nbJets', b'flavL1', b'flavL2']
    matas = ["Sample", "Description" ] #"Regime",
    branches_list_MC = branches_list_data
    bins = range(100,1400,50)
    bins = range(100,145,5)
    #bins = np.linspace(100,140,16)
    #all_sample = []
    rescaledic = None
    if rescale:
        rescaledic = loadnorm("C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/confignormonly.cfg",
        "C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/GlobalFit_fitres_unconditionnal_mu0_normonly.txt")
    if slopecorrection:
        p1s = []
        p2s = []
        bottom = 0
        middle = 0
        top = 0
        with open("output/slopefit/" + "pTV-mbbcut-2tagpolyfitresult.csv") as f:
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
        print(p1s,p2s)
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
    #print(rescaledic)
    all_sample_after = [each for each in all_sample]

    if rescale:
        print("Performing rescale...")
        for i in range(len(all_sample_after)):
            for each_key in rescaledic.keys():
                if 'ALL' in rescaledic[each_key]:
                    factor = rescaledic[each_key]['ALL'] + 1
                    mask = all_sample_after[i].mata["Sample"] == zlib.adler32(each_key.encode())
                    if True in mask:
                        all_sample_after[i].rescale(factor, mask)
    
    if slopecorrection:
        print("Performing slope correction...")
        for i in range(len(all_sample_after)):
            mask1 = all_sample_after[i].data[b'pTV']/1000. < middle
            mask2 = all_sample_after[i].data[b'pTV']/1000. >= middle
            mask2 = np.logical_and(all_sample_after[i].data[b'pTV']/1000. < top, mask2)
            print("before", all_sample_after[i].weight.sum())
            if True in mask1:
                all_sample_after[i].weight[mask1] = all_sample_after[i].weight[mask1] * (poly(all_sample_after[i].data[b'pTV'][mask1]/1000., p1s))
            if True in mask2:
                all_sample_after[i].weight[mask2] = all_sample_after[i].weight[mask2] * (poly(all_sample_after[i].data[b'pTV'][mask2]/1000., p2s))
            print("after", all_sample_after[i].weight.sum())

    bins = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300]
    

    title3="highmBBcr 2 btags"
    direct = "output/t_make_plot/"
    name = "-highmbbcut-2tag"
    if rescale:
        direct = "output/t_make_plot_rescale/"
    if slopecorrection:
        direct = "output/t_make_plot_rescale_slopecorrection/"
    #bins = [0, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950]
    chi2, nod = stackplot(all_sample_after,b'pTV',bins,1000.,
            xlabel=r"$p_{TV}[GeV]$", title3=title3, filename=direct + "pTV" + name, print_height=True,
            title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
    print("pTV", chi2, nod)
    bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1150, 1350, 1550, 1800]
    chi2, nod = stackplot(all_sample_after,b'mVH',bins,1000.,
            xlabel=r"$m_{VH}[GeV]$", title3=title3, filename=direct + "mVH" + name, print_height=True,
            title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
    print("mVH", chi2, nod)
    bins = range(30, 310, 10)
    stackplot(all_sample_after,b'mBBres',bins,1000.,
            xlabel=r"$m_{BB}[GeV]$", title3=title3, filename=direct + "mBBres" + name, print_height=True,
            title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
    # bins = range(0,1000,50)
    # stackplot(all_sample_after,b'pTB1',bins,1000.,
    #         xlabel=r"$pt_{b1}[GeV]$", title3="loose selection, 2 btags", filename="pTB1", print_height=True,
    #         title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True)
    # bins = range(5,80,3)
    # stackplot(all_sample_after,b'MET',bins,1000.,
    #         xlabel=r"$pt_{L2}[GeV]$", title3="muon same charge srcut, < 2 btag", filename="ptL2_3_test", print_height=True,
    #         title2=t2,auto_colour=False, limit_y = 0.5, upper_y=3.0, log_y=False)

    # bins = np.linspace(0, 3.15, num=20)
    # stackplot(all_sample_after,b'delphi1',bins,1.,
    #         xlabel=r"$delphi1$", title3="l same charge srcut, < 3 btag",title4 = "PTl2 < 20 GeV" ,filename="delphi1l", print_height=True,
    #         title2=t2,auto_colour=False, limit_y = 0.5, upper_y=3, log_y=False, blind = False) 
    # stackplot(all_sample_after,b'delphi2',bins,1.,
    #         xlabel=r"$delphi2$", title3="l same charge srcut, < 3 btag",title4 = "PTl2 < 20 GeV" ,filename="delphi2l", print_height=True,
    #         title2=t2,auto_colour=False, limit_y = 0.5, upper_y=3, log_y=False, blind = False) 
    bins = range(0,120,4)
    # stackplot(all_sample_after,b'MT',bins,1000.,
    #         xlabel=r"$MT[GeV]$", title3="mu same charge srcut",title4 = "PTl2 > 20 GeV" ,filename="MTmu_l1_h", print_height=True,
    #         title2=t2,auto_colour=False, limit_y = 0.5, upper_y=3, log_y=False, blind = False)

    # stackplot(all_sample_after,b'MT',bins,1000.,
    #         xlabel=r"$MT[GeV]$", title3="lep same charge srcut",title4 = "PTl2 < 20 GeV" ,filename="MTl_l2_full", print_height=True,
    #         title2=t2,auto_colour=False, limit_y = 0.5, upper_y=3, log_y=False, blind = False)

    # bins = range(0,120,4)
    # stackplot(all_sample_after,b'MT',bins,1000.,
    #         xlabel=r"$MT[GeV]$", title3="mu same charge srcut, < 2 btag",title4 = "PTl2 > 20 GeV" ,filename="MTmu_l1_h", print_height=True,
    #         title2=t2,auto_colour=False, limit_y = 0.5, upper_y=3, log_y=False, blind = False)

    # bins = range(0,120,3)
    # stackplot(all_sample_after,b'MT',bins,1000.,
    #         xlabel=r"$MT[GeV]$", title3="mu same charge srcut",title4="ptl2<15 GeV", filename="MT", print_height=True,
    #         title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=False)
    # bins = range(0,100,2)
    # stackplot(all_sample_after,b'PhiMET',bins,1000.,
    #         xlabel=r"$pt_{MT}[GeV]$", title3="muon same charge srcut, 0, 1 btag", filename="MT", print_height=True,
    #         title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=False)
    # bins = range(0,100,2)
    # stackplot(all_sample_after,b'PhiL1',bins,1000.,
    #         xlabel=r"$pt_{MT}[GeV]$", title3="muon same charge srcut, 0, 1 btag", filename="MT", print_height=True,
    #         title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=False)

'''
    #additional jobs
    wdata = 0
    zdata = 0
    bbhdata = 0
    for each in all_sample_after:
        if each.alias == "Wlvjet":
            wdata = each
        if each.alias == "Zlljet":
            zdata = each
        print(each.alias)
        if each.alias == "AZh":
            bbhdata = each
    bins = range(0,100,2)
    histplot([[wdata],[zdata],[bbhdata]], b'ptL2', bins, ["W+jets", "Z+jets","bbA"], 1000,norm=True,filename="misidbbA", title3="SR, 2 btags",  xlabel=r"$pt_{l2}[GeV]$")
    #histplot([[wdata],[zdata],[bbhdata]], b'ptL2', bins, ["W+jets", "Z+jets","ggA"], 1000,norm=True,filename="misid", title3="SR, 2 btags",  xlabel=r"$pt_{l2}[GeV]$")
'''