import uproot
import numpy as np
import package
from package.events import *
from package.stackplot import *
import pandas as pd
from package.curveplot import histplot, curveplot, histplot_withsub
from package.loadnormfactor import *

def fake_data(bins, hist, variable, stat2, sys2, alias, color, rescaledic=None, rescaledicalias=None):
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
    if rescaledic is not None and rescaledicalias is not None:
        if rescaledicalias not in rescaledic:
            print("Warning: " + rescaledicalias + "rescale is not applied.")
            return sample

        if "ALL" in rescaledic[rescaledicalias]:
            sample = sample * (1 + rescaledic[rescaledicalias]["ALL"])

    return sample

bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1150, 1350, 1550, 1800]
bins = range(0,900,30)
file = uproot.open("../sample/histo/run2_fit.root")
region = "_mBBcr_"
variable = "pTV"
btag = "1tag2pjet"
rescale = True
#systematics = ["SysMODEL_VHFJets_MadGraph", "SysMODEL_VhlJets_MadGraph", "SysMODEL_VlJets_MadGraph", "SysMODEL_ZHFJets_MadGraph", "SysMODEL_ZhlJets_MadGraph", "SysMODEL_ZlJets_MadGraph"]

systematics = ["SysMODEL_VHFJets_MadGraph", "SysMODEL_VhlJets_MadGraph", "SysMODEL_ZHFJets_MadGraph", "SysMODEL_ZhlJets_MadGraph"]
#systematics = ["SysMODEL_ZHFJets_MadGraph", "SysMODEL_ZhlJets_MadGraph"]
#systematics = ["SysMODEL_VHFJets_MadGraph", "SysMODEL_VhlJets_MadGraph"]
mc_Wlvjet = ["Wl", "Wcl", "Wbl", "Wbb", "Wbc", "Wcc"]
mc_Zlljet = ["Zcc", "Zcl", "Zbl", "Zbc", "Zl", "Zbb"]
mc_tt_bar = ["ttbar"]
mc_singletop = ["stopWt", "stops", "stopt"]
mc_Diboson = ["WZ", "WW", "ZZ", "ggZZ", "ggWW"]
sm_Higgs = ["ggZllH125", "qqZllH125"]
file_name_array = [mc_Diboson, mc_tt_bar,  mc_singletop, mc_Zlljet, mc_Wlvjet]#, sm_Higgs]
file_name_array = [num for elem in file_name_array for num in elem]
#alias = ["Diboson", "ttbar", "singletop", "Zlljet", "Wlvjet", "smHiggs"]
#colors = ['g',    'yellow', 'tab:orange','royalblue', 'm',     'r']

rescaledic = None
if rescale:
    rescaledic = loadnorm("C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/configLLBB_190517_HVT_PRSR_MCstat0_Prun1_finalNPtreatment_RFfixC0_2000.cfg",
    "C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/GlobalFit_fitres_unconditionnal_mu0.txt")
allhistname = file["Systematics"].keys()
allhistname_nominal = file.keys()
nominal_dic = {}

nominal = None
data = None
print("loading Nominal...")
# load nominal and data
for each_mc_name in file_name_array + ["data"]:
    thename = None
    for each_name in allhistname_nominal:
        each_name_str = each_name.decode("utf-8") 
        if each_mc_name not in each_name_str:
            continue
        if '_' + variable + ';' not in each_name_str:
            continue
        if region not in each_name_str:
            continue
        if btag not in each_name_str:
            continue
        thename = each_name
        break
    if thename is None:
        print("Warning: cannot find: nominal " + each_mc_name)
        continue

    histpd = file[thename].pandas()
    edge = []
    count = []
    error = []
    for row in histpd.head(10000).itertuples():
        edge.append(row.Index.left)
        count.append(row.count)
        error.append(row.variance)
    edge.append(row.Index.right)
    if each_mc_name == "data":
        data = fake_data(edge, count, variable, error, None, "data", 'k', rescaledic, each_mc_name)
        continue
    data_tem = fake_data(edge, count, variable, error, None, "sys", 'b', rescaledic, each_mc_name)
    nominal_dic[each_mc_name] = data_tem
    if nominal is None:
        nominal = fake_data(edge, count, variable, error, None, "Nominal", 'r', rescaledic, each_mc_name)
    else:
        nominal = nominal + fake_data(edge, count, variable, error, None, "Nominal", 'r', rescaledic, each_mc_name)
if data is None:
    print("Warning: no data found!")
# histplot_withsub([[nominal], [nominal],[data]], variable, bins,labels = ["nominal", "nominal2","data"] )
# exit(1)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("loading Systematics...")
datasys = None
alreadyfound = []
# load nominal and data
for each_mc_name in file_name_array:
    thename = None
    for each_name in allhistname:
        each_name_str = each_name.decode("utf-8")
        if each_mc_name not in each_name_str:
            continue
        if '_' + variable + '_' not in each_name_str:
            continue
        if region not in each_name_str:
            continue
        if btag not in each_name_str:
            continue
        findit = False
        for each_sys in systematics:
            if each_sys in each_name_str:
                findit = True
        if not findit:
            continue
        thename = each_name
        #print(thename)
        break
    if thename is None:
        if each_mc_name not in nominal_dic:
            #print("Warning: cannot find: nominal " + each_mc_name)
            continue
        print("info: cannot find sys " + each_mc_name + ", replacing by nominal.")
        if datasys is None:
            datasys = nominal_dic[each_mc_name]
        else:
            datasys = datasys + nominal_dic[each_mc_name]
        continue
    # if each_mc_name in alreadyfound:
    #     print("Error: double count of " + each_mc_name)
    #     exit(1)
    #     continue
    alreadyfound.append(each_mc_name)

    histpd = file["Systematics"][thename].pandas()
    edge = []
    count = []
    error = []
    for row in histpd.head(10000).itertuples():
        edge.append(row.Index.left)
        count.append(row.count)
        error.append(row.variance)
    edge.append(row.Index.right)
    if each_mc_name == "data":
        continue
    if datasys is None:
        datasys = fake_data(edge, count, variable, error, None, "sys", 'b', rescaledic, each_mc_name)
    else:
        datasys = datasys + fake_data(edge, count, variable, error, None, "sys", 'b', rescaledic, each_mc_name)

histplot_withsub([[nominal],[data], [datasys]], variable, bins,labels =["nominal","data", "sys"], xlabel=r"$p_{TV}[GeV]$" )

# height_nominal, sigma2_mominal = nominal.binned_weight_variation(variable,bins,1)
# height_data, sigma2_data = data.binned_weight_variation(variable,bins,1)
# height_datasys, sigma2_datasys = datasys.binned_weight_variation(variable,bins,1)

# height_nominal = np.array(height_nominal)
# height_data = np.array(height_data)
# height_datasys = np.array(height_datasys)

# height_nominal = height_nominal/height_data
# height_datasys = height_datasys/height_data
# height_data = height_data/height_data


# curveplot([height_nominal,height_datasys, height_data], [bins, bins, bins], labels=["nominal", "sys", "data"])