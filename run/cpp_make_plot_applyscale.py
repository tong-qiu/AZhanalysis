import os
import sys
import json
import copy
import numpy as np
lib_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.append(lib_path)
# from ..package import *
import package
from package.events import *
from package.stackplot import *
from termcolor import colored

def fake_data(bins, hist, variable, stat2, sys2, alias, color):
    is_fake_data = True
    new_data = []
    weight = []
    if "data" in alias:
        is_fake_data = False 
        for i, each in enumerate(hist):
            for j in range(int(each)):
                new_data.append((bins[i]+bins[i+1])/2.)
                weight.append(1)
    else:
        for i, each in enumerate(hist):
            new_data.append((bins[i]+bins[i+1])/2.)
            weight.append(each)
    new_data = {variable: np.array(new_data)}
    sample = Events(new_data, weight, alias=alias, colour=color, fake_data=is_fake_data)
    sample.fake_stat2_per_event = stat2
    sample.fake_sys2_per_event = sys2
    return sample

def load_hist(path, cut, variable):
    i = 0
    allsys_up = dict()
    allsys_down = dict()
    allsys_other = dict()
    # systematics, cut name, variable name, [height, sigma2]
    if os.path.isfile(path + ".json"):
        print(os.path.isfile(path + ".json"),path)
        with open(path + ".json", "r") as f:
            for each_line in f:
                if i%2 == 0:
                    key = each_line[0:len(each_line)-1]
                if i%2 == 1:
                    # cut name, variable name, [height, sigma2]
                    hist_tem = json.loads(each_line)
                    if variable not in hist_tem[cut]: 
                        return (False, False, False)
                    if key == "Nominal":
                        hist = np.array(hist_tem[cut][variable][0])
                        stat2 = np.array(hist_tem[cut][variable][1])
                    elif "up" == key[len(key)-2:len(key)]:
                        allsys_up[key[0:len(key)-2]] = np.array(hist_tem[cut][variable][0])
                    elif "down" == key[len(key)-4:len(key)]:
                        allsys_down[key[0:len(key)-4]] = np.array(hist_tem[cut][variable][0])
                    else:
                        allsys_other[key] = np.array(hist_tem[cut][variable][0])
                i += 1
    else:
        print("warning: cannot find " + path + ".json.")
        return (False, False, False)

    keys = list(allsys_up.keys())
    for each_key in keys:
        if each_key not in allsys_down:
            allsys_other[each_key] = allsys_up[each_key]
            allsys_up.pop(each_key)

    keys = list(allsys_down.keys())
    for each_key in keys:
        if each_key not in allsys_up:
            allsys_other[each_key] = allsys_down[each_key]
            allsys_down.pop(each_key)
    # calculate systematics
    sys2 = np.zeros(len(hist))
    for each_key in allsys_up.keys():
        tem_sys_1 = abs(allsys_up[each_key] - allsys_down[each_key])
        tem_sys_2 = abs(allsys_up[each_key] - hist)
        tem_sys_3 = abs(hist - allsys_down[each_key])
        tem_sys_1 = np.maximum(tem_sys_1,tem_sys_2)
        tem_sys = np.maximum(tem_sys_1,tem_sys_3)
        sys2 += (tem_sys/2.) ** 2
    return (hist, stat2, sys2)

def bincheck(orginal, userbin):
    for i, each in enumerate(userbin):
        if each in orginal:
            continue
        if each < orginal[0] or each > orginal[-1]:
            continue
        print(colored("Warning: possible bin migration! Original bin is: " + str(orginal[0]) + ":" + str(orginal[-1]) + ":" + str(orginal[1] - orginal[0]),'red'))
        return 0


if __name__ == '__main__':
    rebin_factor = 1
    filename = "a_mVH_SR_2tag2pjet-_.json"
    #filename = "e_mVH_mBBcr_2tag2pjet-_.json"
    sub_filename = filename.split('_')
    period = sub_filename[0]
    variable_name = sub_filename[1]
    region = sub_filename[2]
    btag = sub_filename[3].split("tag")[0].replace('p', '+')
    njet = sub_filename[3].split("tag")[1].split('jet')[0].replace('p', '+')

    t2 = r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$"
    if period == "a":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$"
    if period == "d":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,43.6\:fb^{-1}}$"
    if period == "e":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,58.5\:fb^{-1}}$"
    if period == "run2":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,138.2\:fb^{-1}}$"
    
    mc_Wlvjet = ["Wl", "Wcl", "Wbl", "Wbb", "Wbc", "Wcc"]
    mc_Zlljet = ["Zcc", "Zcl", "Zbl", "Zbc", "Zl", "Zbb"]
    mc_tt_bar = ["ttbar"]
    mc_singletop = ["stopWt", "stops", "stopt"]
    mc_Diboson = ["WZ", "WW", "ZZ", "ggZZ", "ggWW"]
    sm_Higgs = ["ggZllH125", "qqZllH125"]
    file_name_array = [mc_Diboson, mc_tt_bar,  mc_singletop, mc_Zlljet, mc_Wlvjet]#, sm_Higgs]
    alias = ["Diboson", "ttbar", "singletop", "Zlljet", "Wlvjet", "smHiggs"]
    colors = ['g',    'yellow', 'tab:orange','royalblue', 'm',     'r']

    all_sample = []
    binning = []
    with open("../../postreader/PlotTool_Root/jsonoutput/"+filename, 'r') as f:
        for i, each_line in enumerate(f):
            jsondic = json.loads(each_line)
            if i == 0:
                sys_done = False
                for each_array, each_alias, each_color in zip(file_name_array, alias, colors):
                    data = "nodata"
                    this_sample = []
                    for each_sample in each_array:
                        if each_sample not in jsondic:
                            print("warning: connot find " + each_sample)
                            continue
                        binning = jsondic[each_sample]["binning"]
                        content = jsondic[each_sample]["content"]
                        stat = np.array(jsondic[each_sample]["stat"])
                        sys = [0 for each in content]
                        if not sys_done:
                            sys = np.array(jsondic["nominal"]["syst"])**2
                            sys_done = True
                        if data == "nodata":
                            data = fake_data(binning,content,"mVH",stat**2,sys,each_alias,each_color)
                        else: 
                            data = data + fake_data(binning,content,"mVH",stat**2,sys,each_alias,each_color)
                    all_sample.append(data)
                #all_sample[0].fake_sys2_per_event = np.array(jsondic["nominal"]["syst"])**2
            if i == 1:
                binning = jsondic["nominal"]["binning"]
                content = jsondic["nominal"]["content"]
                stat = np.array(jsondic["nominal"]["stat"])
                sys = [0 for each in content]
                data = fake_data(binning,content,"mVH",stat**2,sys,"data",'k')
                all_sample.append(data)
        bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1150, 1350, 1550, 1800]
        bincheck(binning, bins)
        all_sample1 = []
        for each in all_sample:
            #print(type(each))
            if isinstance(each, str):
                continue
            all_sample1.append(each)
        stackplot(all_sample1,variable_name,bins,1.,
                xlabel=r"$m_{VH}[GeV]$", title3="2 lep.," + btag +" b-tag " + region, filename="output/cpp_make_plot/" + filename[0:-5], print_height=False,
                title2=t2,auto_colour=False, limit_y = 0.6, upper_y=2.6, sys=True, log_y=True)
    '''
    sample_list = {"Wl", "Wcl", "Wbl", "Wbb", "Wbc", "Wcc", "WZ", "WW",              "Zcc", "Zcl", "Zbl", "Zbc", "Zl", "Zbb", "ZZ", "stopWt", "stops", "stopt", "ttbar", "ggZllH125", "qqZllH125", "stopWt_dilep"}
    colors =      ['g', 'g',    'g',    'g',  'g',   'g',  'tab:orange','tab:orange', 'royalblue','royalblue','royalblue','royalblue','royalblue','royalblue','royalblue','yellow','yellow','yellow','yellow', 'yellow','yellow','yellow']
    with open("mVHsr.txt","r") as f:
        for each_line in f:
            sample = json.loads(each_line)
    nominals = []
    for each, color in zip(sample_list,colors):
        if each not in sample:
            print("Warning: cannot find ", each)
            continue
        binning = sample[each]["binning"]
        content = sample[each]["content"]
        stat = np.array(sample[each]["stat"])
        sys = [0 for each in content]
        data = fake_data(binning,content,"mVH",stat**2,sys,each,color)
        nominals.append(data)
    nominals[0].fake_sys2 = np.array(sample["nominal"]["syst"])**2

    with open("mVHsrdata.txt","r") as f:
        for each_line in f:
            sample = json.loads(each_line)
    binning = sample["nominal"]["binning"]
    content = sample["nominal"]["content"]
    stat = np.array(sample["nominal"]["stat"])
    sys = [0 for each in content]
    data = fake_data(binning,content,"mVH",stat**2,sys,"data",color)
    nominals.append(data)

    bins = range(0,2000,120)
    stackplot(nominals,'mVH',bins,1.,
            xlabel=r"$m_{VH}[GeV]$", title3="2 lep., 2 b-tag", filename="test", print_height=True,
            title2="t2",auto_colour=True, limit_y = 0.4, upper_y=2.6, sys=True, log_y=True)
    '''
