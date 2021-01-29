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
from cut import *
from package.stackplot import *
from curveplot import *
from cutstring import *
import multiprocessing
from package.loadnormfactor import *

def fitfunction_real(x0, p0, p1, p2):
    out = p0 + p1 * x0 + p2 * x0**2
    if out < 0:
        return 0
    if x0 > 1700:
        return fitfunction_real(1700, p0, p1, p2)
    return out

def fitfunction(x0, p0, p1, p2):
    y = np.zeros(len(x0))
    out = p0 + p1 * x0 + p2 * x0**2
    y += (out ) * (out > 0) * (x0 < 1700)
    y += fitfunction_real(1700, p0, p1, p2) * (x0 >= 1700)
    return y

def get_slope_correction(path):
    p1s = []
    with open(path) as f:
        for each_line in f:
            p1s = each_line.split(',')
            for i in range(len(p1s)-1):
                p1s[i] = float(p1s[i])
            print(p1s)
            break
    return p1s

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

# function to load easytrees and perform event selections
def stack_cxaod(sample_directory, each_names, each_alias, each_color, branches_list_data, debug, cut, m_allsamples, matas=None):
    sample = load_CxAODs(sample_directory,each_names,branches_list_data, debug, 
                        colour=each_color,alias=each_alias,matanames=matas)
    if not sample:
        print("Warning: No "+each_alias+" samples found!")
    if cut and sample:
        # choose event selections here

        # select resolved event using the "Regime" branch
        # This is a string-based event selection. The easytree branch which contains the string 
        # should be defined in the "matas" list. 
        # The event selection criterion is defined in ml/cutstring.py
        # sample.matacut(s_resolved)
        sample.matacut(s_merged)
        sample.matacut(s_mbbcr)
        sample.ptfj1()

        # select events with certain number of b-tagged jets
        # This is a value-based event selection. The easytree branch which contains the values 
        # should be defined in the "branches_list_data" list.
        # The event selection criterion is defined in ml/mlcut.py
        # ntag = 2
        # sample.cut_parameter(cut_btag_is, ntag)

        # other user defined event selection
        # sample.cut(cut_basic)

        m_allsamples.append(sample)
    if not cut:
        # sample.cut_parameter(cut_btag_more, 0)
        # sample.add_region()
        # sample.add_regime()
        # sample.mata = 0
        m_allsamples.append(sample)
    return 0

if __name__ == '__main__':
    # only load limited number of the events if debug
    debug = False
    # Do event selection?
    cut = True
    # save event after selection as root file?
    saveevent = True
    rescale = True
    slopecorrection = False
    tag = "a"
    # directory of the easytrees
    sample_directory = ["../sample/a/", "../sample/d/", "../sample/e/"]
    data = ["data16", "data15", "data17", "data18"]
    # Text on the plot. Delete if not needed.
    t2 = r"$\mathit{\sqrt{s}=13\:TeV,139\:fb^{-1}}$"


    # rescaledic = loadnorm("../fitconfig/normonly12tag.cfg", "../fitconfig/normonly12tag.txt")
    # #rescaledic = loadnorm("../fitconfig/hvtcr/config_cr.cfg", "../fitconfig/hvtcr/GlobalFit_fitres_conditionnal_mu0.txt")
    # print("Performing rescale...")
    # for each_key in rescaledic.keys():
    #     if 'ALL' in rescaledic[each_key]:
    #         factor = rescaledic[each_key]['ALL'][0]
    #         error = rescaledic[each_key]['ALL'][1]**2
    #         if '1pfat0pjet' in rescaledic[each_key]:
    #             factor += rescaledic[each_key]['1pfat0pjet'][0] - 1
    #             error = rescaledic[each_key]['1pfat0pjet'][1]**2
    #         print(each_key, factor, error**0.5)
    # exit(1)

    # if tag == "a":
    #     t2 = r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$"
    #     sample_directory = ["../sample/CxAOD32_06" + tag + "/"]
    #     data = ["data16", "data15"]
    # if tag == "d":
    #     t2 = r"$\mathit{\sqrt{s}=13\:TeV,44.3\:fb^{-1}}$"
    #     sample_directory = ["../sample/CxAOD32_06" + tag + "/"]
    #     data = ["data17"]
    # if tag == "e":
    #     t2 = r"$\mathit{\sqrt{s}=13\:TeV,59.9\:fb^{-1}}$"
    #     sample_directory = ["../sample/CxAOD32_06" + tag + "/"]
    #     data = ["data18"]
    # if tag == "run2":
    #     t2 = r"$\mathit{\sqrt{s}=13\:TeV,140.3\:fb^{-1}}$"
    #     sample_directory = ["../sample/CxAOD32_06a/", "../sample/CxAOD32_06d/", "../sample/CxAOD32_06e/"]
    #     data = ["data16", "data15", "data17", "data18"]

    # define root files of each backgorund here. 
    # background files. Do not change!
    # --------------------------------------------------------------------
    mc_Wlvjet = ["Wenu_Sh221", "WenuB_Sh221", "WenuC_Sh221", "WenuL_Sh221", "Wmunu_Sh221", "WmunuB_Sh221", "WmunuC_Sh221", "WmunuL_Sh221", "Wtaunu_Sh221", "WtaunuB_Sh221", "WtaunuC_Sh221", "WtaunuL_Sh221"]
    mc_Zlljet1 = ["Zee_Sh221", "ZeeB_Sh221"]
    mc_Zlljet2 = ["ZeeC_Sh221", "ZeeL_Sh221"]
    mc_Zlljet3 = ["Zmumu_Sh221", "ZmumuB_Sh221"]
    mc_Zlljet4 = ["ZmumuC_Sh221", "ZmumuL_Sh221"]
    mc_Zlljet5 = ["Ztautau_Sh221", "ZtautauB_Sh221", "ZtautauC_Sh221", "ZtautauL_Sh221", "Znunu_Sh221", "ZnunuB_Sh221", "ZnunuC_Sh221", "ZnunuL_Sh221"]
    mc_tt_bar = [ "ttbar_nonallhad_PwPy8", "ttbar_allhad_PwPy8", "ttbar_dilep_PwPy8"] # must include all three
    mc_singletop = ["stops_PwPy8", "stopt_PwPy8", "stopWt_PwPy8"]
    mc_Diboson = ["WqqWlv_Sh221", "WqqZll_Sh221", "WqqZvv_Sh221", "ZqqZll_Sh221", "ZqqZvv_Sh221", "WlvZqq_Sh221", "ggZqqZll_Sh222", "ggWqqWlv_Sh222"]
    sm_Higgs = ["bbHinc_aMCatNLOPy8", "ggHinc_PwPy8", "ggZllHbb_PwPy8","ggZllHcc_PwPy8","ggZvvHbb_PwPy8","ggZvvHcc_PwPy8","qqWlvHbbJ_PwPy8MINLO","qqWlvHccJ_PwPy8MINLO","qqZllHbbJ_PwPy8MINLO","qqZllHccJ_PwPy8MINLO","qqZvvHbbJ_PwPy8MINLO","qqZvvHccJ_PwPy8MINLO"]
    ttV = ["ttV_aMCatNLOPy8_alternative-0"]
    # --------------------------------------------------------------------

    # signal files
    HVT = ["HVT"]

    # signal/background to be loaded
    file_name_array = [data, mc_Diboson, mc_tt_bar,  mc_singletop, mc_Zlljet1, mc_Zlljet2, mc_Zlljet3, mc_Zlljet4, mc_Zlljet5, mc_Wlvjet, HVT, ttV, sm_Higgs]
    # file_name_array = [HVT]
    # choose a name for your backgrounds
    alias = ["data", "Diboson", "ttbar", "singletop", "Zlljet", "Zlljet", "Zlljet", "Zlljet", "Zlljet", "Wlvjet", "HVT", 'ttV', 'sm_Higgs']
    # alias = ["HVT"]
    # Colours of each background on the plot. Delete if not needed.
    colors = [None, 'g', 'yellow', 'tab:orange', 'royalblue', 'royalblue', 'royalblue', 'royalblue', 'royalblue', 'm', "r", 'dimgrey', 'teal']

    # Variables to load.
    branches_list_data = [b"EventWeight", b'mVH', b'nTags', b'nbTagsInFJ', b'nbTagsOutsideFJ', b'ptTrkJetsinFJ1', b'ptTrkJetsinFJ2', b'fj1px', b'fj1py']
    # Strings to load.
    matas = ["Regime", "Description", "Sample"]
    branches_list_MC = copy.deepcopy(branches_list_data)
    branches_list_MC.append(b'MCChannelNumber')

    if slopecorrection:
        p1s = get_slope_correction("ptfj1mbbcut-1tagpolyfitresult.csv")
    # Load samples and event selection
    # Do not change!
    # --------------------------------------------------------------------
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
    all_sample_after = new_data_list
    # --------------------------------------------------------------------
    # end of sample loading

    # all_sample_after is a list of event.Events objects. event.Events is a class which stores 
    # all the events of a background.
    # variables and methods of event.Events class you may need to use:

    # Event.data: A dictionary contains all the variables requested in the branches_list_data list.
    # {b'mBBres': numpy array, b'mBB': numpy array, ....}
    # Event.weight: a numpy array contains the weight of each event
    # Event.alias: a string which stores the type of the background
    # Event.__add__: merged two event.Events object. Example: merged = Eventsobject1 + Eventsobject2


    if rescale:
        rescaledic = loadnorm("../fitconfig/normonly12tag.cfg", "../fitconfig/normonly12tag.txt")
        #rescaledic = loadnorm("../fitconfig/hvtcr/config_cr.cfg", "../fitconfig/hvtcr/GlobalFit_fitres_conditionnal_mu0.txt")
        print("Performing rescale...")
        for i in range(len(all_sample_after)):
            for each_key in rescaledic.keys():
                if 'ALL' in rescaledic[each_key]:
                    factor = rescaledic[each_key]['ALL'] + 1
                    if '1pfat0pjet' in rescaledic[each_key]:
                        factor += rescaledic[each_key]['1pfat0pjet']
                    mask = all_sample_after[i].mata["Sample"] == zlib.adler32(each_key.encode())
                    if True in mask:
                        all_sample_after[i].rescale(factor, mask)

    # save as root files for MVA. delete if not needed.
    if saveevent:
        print("Saving events for MVA training...")
        backgroundlist = None
        datalist = []
        signallist = []
        backgroundlist = []
        for content in all_sample_after:
            if "data" in content.alias:
                datalist.append(content)
            elif "HVT" in content.alias:
                signallist.append(content)
            else:
                backgroundlist.append(content)
        
        test = get_signalid("HVT")
        out = splitesamples(signallist[0], test)


    correctionname = ""
    if slopecorrection:
        correctionname = "-correction"
        print("Performing slope correction...")
        for i in range(len(backgroundlist)):
            backgroundlist[i].weight = backgroundlist[i].weight * (fitfunction(backgroundlist[i].data[b'ptfj1']/1000., p1s[0], p1s[1], p1s[2]))



    datatest = copy.deepcopy(datalist)
    bkgtest = copy.deepcopy(backgroundlist)
    sumbkg = 0
    sumdata = 0
    for each in bkgtest + datatest:
        each.cut_parameter(cut_btagin_is, 1)
        # each.cut_parameter(cut_btagout_less, 0)
        if "data" not in each.alias:
            sumbkg += sum(each.weight)
        else:
            sumdata += len(each.weight)
    print("Making plots...")
    direct = ""
    name = "mbbcut-1tag" + correctionname
    bins = [200, 220, 240, 260, 290, 320, 350, 380, 410,  450, 550, 650, 750, 850, 1000, 1200, 1400, 1800]
    stackplot(bkgtest + datatest,b'ptfj1',bins,1000.,
    xlabel=r"$p_{T}^{FJ}[GeV]$", filename=direct + "ptfj1" + name, print_height=True,
    title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.5, log_y=True, title3="2 lep. 1 btag merged", chi2=True)
    bins = [260, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400, 1800, 2200]
    stackplot(bkgtest + datatest,b'mVH',bins,1000.,
    xlabel=r"$m_{VH}[GeV]$", filename=direct + "mvh" + name, print_height=True,
    title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.5, log_y=True, title3="2 lep. 1 btag merged", chi2=True)


    datatest = copy.deepcopy(datalist)
    bkgtest = copy.deepcopy(backgroundlist)
    sumbkg = 0
    sumdata = 0
    for each in bkgtest + datatest:
        each.cut_parameter(cut_btagin_is, 1)
        each.cut_parameter(cut_btagout_less, 0)
        if "data" not in each.alias:
            sumbkg += sum(each.weight)
        else:
            sumdata += len(each.weight)
    print("Making plots...")
    direct = ""
    name = "mbbcut-1tag-noadd" + correctionname
    bins = [200, 220, 240, 260, 290, 320, 350, 380, 410,  450, 550, 650, 750, 850, 1000, 1200, 1400, 1800]
    stackplot(bkgtest + datatest,b'ptfj1',bins,1000.,
    xlabel=r"$p_{T}^{FJ}[GeV]$", filename=direct + "ptfj1" + name, print_height=True,
    title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.5, log_y=True, title3="2 lep. 1 btag merged, no add.", chi2=True)
    bins = [260, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400, 1800, 2200]
    stackplot(bkgtest + datatest,b'mVH',bins,1000.,
    xlabel=r"$m_{VH}[GeV]$", filename=direct + "mvh" + name, print_height=True,
    title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.5, log_y=True, title3="2 lep. 1 btag merged, no add.", chi2=True)

    datatest = copy.deepcopy(datalist)
    bkgtest = copy.deepcopy(backgroundlist)
    sumbkg = 0
    sumdata = 0
    for each in bkgtest + datatest:
        each.cut_parameter(cut_btagin_is, 1)
        each.cut_parameter(cut_btagout_more, 1)
        if "data" not in each.alias:
            sumbkg += sum(each.weight)
        else:
            sumdata += len(each.weight)
    print("Making plots...")
    direct = ""
    name = "mbbcut-1tag-add" + correctionname
    bins = [200, 220, 240, 260, 290, 320, 350, 380, 410,  450, 550, 650, 750, 850, 1000, 1200, 1400, 1800]
    stackplot(bkgtest + datatest,b'ptfj1',bins,1000.,
    xlabel=r"$p_{T}^{FJ}[GeV]$", filename=direct + "ptfj1" + name, print_height=True,
    title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.5, log_y=True, title3="2 lep. 1 btag merged, with add.", chi2=True)
    bins = [260, 300, 400, 500, 600, 700, 800, 1000, 1200, 1400, 1800, 2200]
    stackplot(bkgtest + datatest,b'mVH',bins,1000.,
    xlabel=r"$m_{VH}[GeV]$", filename=direct + "mvh" + name, print_height=True,
    title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.5, log_y=True, title3="2 lep. 1 btag merged, with add.", chi2=True)


    if not slopecorrection:
        datatest = copy.deepcopy(datalist)
        bkgtest = copy.deepcopy(backgroundlist)
        sumbkg = 0
        sumdata = 0
        for each in bkgtest + datatest:
            each.cut_parameter(cut_btagin_is, 2)
            # each.cut_parameter(cut_btagout_less, 0)
            if "data" not in each.alias:
                sumbkg += sum(each.weight)
            else:
                sumdata += len(each.weight)
        print("Making plots...")
        direct = ""
        name = "mbbcut-2tag" + correctionname
        bins = [200, 250, 350, 450, 550, 650, 750, 850, 1000, 1200, 1400]
        stackplot(bkgtest + datatest,b'ptfj1',bins,1000.,
        xlabel=r"$pT_{FJ}[GeV]$", filename=direct + "ptfj1" + name, print_height=True,
        title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.5, log_y=True, title3="2 lep. 2 btag merged", chi2=True)


        datatest = copy.deepcopy(datalist)
        bkgtest = copy.deepcopy(backgroundlist)
        sumbkg = 0
        sumdata = 0
        for each in bkgtest + datatest:
            each.cut_parameter(cut_btagin_is, 2)
            each.cut_parameter(cut_btagout_less, 0)
            if "data" not in each.alias:
                sumbkg += sum(each.weight)
            else:
                sumdata += len(each.weight)
        print("Making plots...")
        direct = ""
        name = "mbbcut-2tag-noadd" + correctionname
        bins = [200, 250, 350, 450, 550, 650, 750, 850, 1000, 1200, 1400]
        stackplot(bkgtest + datatest,b'ptfj1',bins,1000.,
        xlabel=r"$pT_{FJ}[GeV]$", filename=direct + "ptfj1" + name, print_height=True,
        title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.5, log_y=True, title3="2 lep. 2 btag merged, no add.", chi2=True)

        datatest = copy.deepcopy(datalist)
        bkgtest = copy.deepcopy(backgroundlist)
        sumbkg = 0
        sumdata = 0
        for each in bkgtest + datatest:
            each.cut_parameter(cut_btagin_is, 2)
            each.cut_parameter(cut_btagout_more, 1)
            if "data" not in each.alias:
                sumbkg += sum(each.weight)
            else:
                sumdata += len(each.weight)
        print("Making plots...")
        direct = ""
        name = "mbbcut-2tag-add" + correctionname
        bins = [200, 250, 350, 450, 550, 650, 750, 850, 1000, 1200, 1400]
        stackplot(bkgtest + datatest,b'ptfj1',bins,1000.,
        xlabel=r"$pT_{FJ}[GeV]$", filename=direct + "ptfj1" + name, print_height=True,
        title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.5, log_y=True, title3="2 lep. 2 btag merged, with add.", chi2=True)