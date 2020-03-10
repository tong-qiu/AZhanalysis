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

def poly(x, argv):
    s = 0
    for i, each in enumerate(argv):
        s += x**i * each
    return s

def fitfunction(x, p0, p1, p2, p4):
    y = np.zeros(len(x))
    y += (p0 + p1 * 40 + p2 * 40**2) * (x <= 40)
    y += (p0 + p1 * x + p2 * x**2) * (x <= p4) * (x > 40)
    y += (p0 + p1 * p4 + p2 * p4**2) * (x > p4)
    return y

def fitfunction_real(x, p0, p1, p2, p4):
    if x < 40:
        return p0 + p1 * 40 + p2 * 40**2
    if x < p4:
        return p0 + p1 * x + p2 * x**2
    return p0 + p1 * p4 + p2 * p4**2

def autobin(data_list, bins, alias=None, variable=b"pTV"):
    new_data = None
    for each in data_list:
        if (alias is None and alias != "data") or each.alias == alias:
            if new_data is None:
                new_data = each
            else:
                new_data = new_data + each
    height, error = new_data.binned_weight_variation(variable, bins, scale=1000.)
    newbin = [bins[0]]

    sum_weight = 0
    sum_error = 0
    for i in range(len(height)):
        if i == len(height) - 1:
            newbin.append(bins[i+1])
            break
        sum_weight += height[i]
        sum_error += error[i]

        if (sum_error**0.5)/sum_weight < 0.35 and sum_weight > 3:
            newbin.append(bins[i+1])
            sum_weight = 0
            sum_error = 0
    return newbin

def autobin_withdata(data_list, bins, alias=None, variable=b"pTV"):
    data_list = copy.deepcopy(data_list)
    new_data = None
    for each in data_list:
        if (alias is None and alias != "data") or each.alias == alias:
            if new_data is None:
                new_data = each
            else:
                new_data = new_data + each
    height, error = new_data.binned_weight_variation(variable, bins, scale=1000.)

    new_data = None
    for each in data_list:
        if each.alias == "data":
            if new_data is None:
                new_data = each
            else:
                new_data = new_data + each
    height_data, error_data = new_data.binned_weight_variation(variable, bins, scale=1000.)

    newbin = [bins[0]]
    sum_weight = 0
    sum_error = 0
    sum_weight_data = 0
    sum_error_data = 0
    for i in range(len(height)):
        if i == len(height) - 1:
            newbin.append(bins[i+1])
            break
        sum_weight += height[i]
        sum_error += error[i]
        sum_weight_data += height_data[i]
        sum_error_data += error_data[i]

        if (sum_error**0.5)/sum_weight < 0.35 and sum_weight > 3 and (sum_error_data**0.5)/sum_weight_data < 0.35:
            newbin.append(bins[i+1])
            sum_weight = 0
            sum_error = 0
    return newbin


ntag = 2

def stack_cxaod(sample_directory, each_names, each_alias, each_color, branches_list_data, debug, cut, m_allsamples, matas=None):
    sample = load_CxAODs(sample_directory,each_names,branches_list_data, debug, 
                        colour=each_color,alias=each_alias,matanames=matas)
    if not sample:
        print("Warning: No "+each_alias+" samples found!")
    if cut and sample:
        sample.matacut(s_resolved)
        # sample.cut(cut_lowmbb)
        #sample.cut(cut_highmbb)
        sample.matacut(s_mbbcr)
        sample.cut_parameter(cut_btag_is, ntag)
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
    slopecorrection = True

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
    mc_tt_bar = [ "ttbar_nonallhad_PwPy8", "ttbar_allhad_PwPy8", "ttbar_dilep_PwPy8"]#"ttbar_nonallhad_PwPy8", , "ttbar_allhad_PwPy8"]#"ttbar_nonallhad_PwPy8"]#, "ttbar_allhad_PwPy8"]
    mc_singletop = ["stops_PwPy8", "stopt_PwPy8", "stopWt_PwPy8", "stopWt_dilep_PwPy8"]
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


    branches_list_data = [b"mBBres", b"EventWeight", b"pTV", b'mVH', b'nTags']
    matas = ["Sample", "Description", "Regime"]
    branches_list_MC = branches_list_data
    bins = range(100,1400,50)
    bins = range(20,200,5)
    #bins = np.linspace(100,140,16)
    #all_sample = []
    rescaledic = None
    if rescale:
        rescaledic = loadnorm("C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/confignormonly.cfg",
        "C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/GlobalFit_fitres_unconditionnal_mu0_normonly.txt")

    if slopecorrection:
        p1s = []
        with open("output/slopefit/" + "pTV-mbbcut-"+str(ntag)+"tagpolyfitresult.csv") as f:
            for each_line in f:
                p1s = each_line.split(',')
                for i in range(len(p1s)-1):
                    p1s[i] = float(p1s[i])
                print(p1s)
                break
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
        all_sample_after_backup = copy.deepcopy(all_sample_after)
        print("Performing slope correction...")
        for i in range(len(all_sample_after)):
            all_sample_after[i].weight = all_sample_after[i].weight * (fitfunction(all_sample_after[i].data[b'pTV']/1000., p1s[0], p1s[1], p1s[2], p1s[3]))
            print("after", all_sample_after[i].weight.sum())
    bins = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300]
    

    title3="mBBcr " + str(ntag) +" btags"
    direct = "output/t_make_plot/"
    name = "-mbbcut-" + str(ntag) +"tag"
    if rescale:
        direct = "output/t_make_plot_rescale/"
    if slopecorrection and rescale:
        direct = "output/t_make_plot_rescale_slopecorrection_new/"
    if slopecorrection and not rescale:
        direct = "output/t_make_plot_slopecorrection/"


    bins = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,  1050, 1100, 1150, 1200, 1250, 1300]
    bins = range(0,1400,20)
    bins = autobin_withdata(all_sample_after, bins, alias="Zlljet", variable=b"pTV")
    print(bins)
    chi2, nod = stackplot(all_sample_after,b'pTV',bins,1000.,
            xlabel=r"$p_{TV}[GeV]$", title3=title3, filename=direct + "pTV" + name, print_height=True,
            title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
    print("pTV", chi2, nod)
    bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1150, 1350, 1550, 1800]
    chi2, nod = stackplot(all_sample_after,b'mVH',bins,1000.,
            xlabel=r"$m_{VH}[GeV]$", title3=title3, filename=direct + "mVH" + name, print_height=True,
            title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
    print("mVH", chi2, nod)
    bins = range(20, 200, 1)
    stackplot(all_sample_after,b'mBBres',bins,1000.,
            xlabel=r"$m_{BB}[GeV]$", title3=title3, filename=direct + "mBB" + name, print_height=True,
            title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=False, printzpjets=True, chi2=True)


    if not slopecorrection:
        all_sample_after1 = copy.deepcopy(all_sample_after)
        all_sample_after2 = copy.deepcopy(all_sample_after)

        for i in range(len(all_sample_after)):
            all_sample_after1[i].cut(cut_lowmbb)
            all_sample_after2[i].cut(cut_highmbb)
        

        title3="lowmBBcr " + str(ntag) +" btags"
        name = "-lowmbbcut-" + str(ntag) +"tag"
        bins = range(0,1400,20)
        bins = autobin_withdata(all_sample_after1, bins, alias="Zlljet", variable=b"pTV")
        print(bins)
        chi2, nod = stackplot(all_sample_after1,b'pTV',bins,1000.,
                xlabel=r"$p_{TV}[GeV]$", title3=title3, filename=direct + "pTV" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
        print("pTV", chi2, nod)
        bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1150, 1350, 1550, 1800]
        chi2, nod = stackplot(all_sample_after1,b'mVH',bins,1000.,
                xlabel=r"$m_{VH}[GeV]$", title3=title3, filename=direct + "mVH" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
        print("mVH", chi2, nod)
        bins = range(20, 200, 1)
        stackplot(all_sample_after1,b'mBBres',bins,1000.,
                xlabel=r"$m_{BB}[GeV]$", title3=title3, filename=direct + "mBB" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=False, printzpjets=True, chi2=True)


        title3="highmBBcr " + str(ntag) +" btags"
        name = "-highmbbcut-" + str(ntag) +"tag"
        bins = range(0,1400,20)
        bins = autobin_withdata(all_sample_after2, bins, alias="Zlljet", variable=b"pTV")
        print(bins)
        chi2, nod = stackplot(all_sample_after2,b'pTV',bins,1000.,
                xlabel=r"$p_{TV}[GeV]$", title3=title3, filename=direct + "pTV" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
        print("pTV", chi2, nod)
        bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1150, 1350, 1550, 1800]
        chi2, nod = stackplot(all_sample_after2,b'mVH',bins,1000.,
                xlabel=r"$m_{VH}[GeV]$", title3=title3, filename=direct + "mVH" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
        print("mVH", chi2, nod)
        bins = range(20, 200, 1)
        stackplot(all_sample_after2,b'mBBres',bins,1000.,
                xlabel=r"$m_{BB}[GeV]$", title3=title3, filename=direct + "mBB" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=False, printzpjets=True, chi2=True)


    if slopecorrection:
        all_sample_after1 = copy.deepcopy(all_sample_after_backup)
        all_sample_after2 = copy.deepcopy(all_sample_after_backup)

        for i in range(len(all_sample_after)):
            all_sample_after1[i].cut(cut_lowmbb)
            all_sample_after2[i].cut(cut_highmbb)

        print("Performing slope correction for low...")
        p1s = []
        with open("output/slopefit/" + "pTV-highmbbcut-"+str(ntag)+"tagpolyfitresult.csv") as f:
            for each_line in f:
                p1s = each_line.split(',')
                for i in range(len(p1s)-1):
                    p1s[i] = float(p1s[i])
                print(p1s)
                break
        for i in range(len(all_sample_after1)):
            all_sample_after1[i].weight = all_sample_after1[i].weight * (fitfunction(all_sample_after1[i].data[b'pTV']/1000., p1s[0], p1s[1], p1s[2], p1s[3]))
            print("after", all_sample_after1[i].weight.sum())

        print("Performing slope correction for high...")
        p1s = []
        with open("output/slopefit/" + "pTV-lowmbbcut-"+str(ntag)+"tagpolyfitresult.csv") as f:
            for each_line in f:
                p1s = each_line.split(',')
                for i in range(len(p1s)-1):
                    p1s[i] = float(p1s[i])
                print(p1s)
                break
        for i in range(len(all_sample_after2)):
            all_sample_after2[i].weight = all_sample_after2[i].weight * (fitfunction(all_sample_after2[i].data[b'pTV']/1000., p1s[0], p1s[1], p1s[2], p1s[3]))
            print("after", all_sample_after2[i].weight.sum())

        title3="lowmBBcr " + str(ntag) +" btags"
        name = "-lowmbbcuthighcorrection-" + str(ntag) +"tag"
        bins = range(0,1400,20)
        bins = autobin_withdata(all_sample_after1, bins, alias="Zlljet", variable=b"pTV")
        print(bins)
        chi2, nod = stackplot(all_sample_after1,b'pTV',bins,1000.,
                xlabel=r"$p_{TV}[GeV]$", title3=title3, filename=direct + "pTV" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
        print("pTV", chi2, nod)
        bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1150, 1350, 1550, 1800]
        chi2, nod = stackplot(all_sample_after1,b'mVH',bins,1000.,
                xlabel=r"$m_{VH}[GeV]$", title3=title3, filename=direct + "mVH" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
        print("mVH", chi2, nod)
        bins = range(20, 200, 1)
        stackplot(all_sample_after1,b'mBBres',bins,1000.,
                xlabel=r"$m_{BB}[GeV]$", title3=title3, filename=direct + "mBB" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=False, printzpjets=True, chi2=True)


        title3="highmBBcr " + str(ntag) +" btags"
        name = "-highmbbcutlowcorrection-" + str(ntag) +"tag"
        bins = range(0,1400,20)
        bins = autobin_withdata(all_sample_after2, bins, alias="Zlljet", variable=b"pTV")
        print(bins)
        chi2, nod = stackplot(all_sample_after2,b'pTV',bins,1000.,
                xlabel=r"$p_{TV}[GeV]$", title3=title3, filename=direct + "pTV" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
        print("pTV", chi2, nod)
        bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1150, 1350, 1550, 1800]
        chi2, nod = stackplot(all_sample_after2,b'mVH',bins,1000.,
                xlabel=r"$m_{VH}[GeV]$", title3=title3, filename=direct + "mVH" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpjets=True, chi2=True)
        print("mVH", chi2, nod)
        bins = range(20, 200, 1)
        stackplot(all_sample_after2,b'mBBres',bins,1000.,
                xlabel=r"$m_{BB}[GeV]$", title3=title3, filename=direct + "mBB" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=False, printzpjets=True, chi2=True)