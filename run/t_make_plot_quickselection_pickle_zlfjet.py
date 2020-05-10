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

from package.easyload import *
from package.events import *
from package.cut import *
from package.stackplot import *
from curveplot import *
from cutstring import *
import multiprocessing
from package.loadnormfactor import *

def pickleit(obj, path):
    outfile = open(path, 'wb')
    pickle.dump(obj, outfile)
    outfile.close()

def unpickleit(path):
    infile = open(path, 'rb')
    output = pickle.load(infile)
    infile.close()
    return output

def poly(x, argv):
    s = 0
    for i, each in enumerate(argv):
        s += x**i * each
    return s

# def fitfunction_real(x, p0, p1, p2, p4):
#     if x < p4:
#         return p0 + p1 * x + p2 * x**2
#     return p0 + p1 * p4 + p2 * p4**2

# def fitfunction(x, p0, p1, p2, p4):
#     y = np.zeros(len(x))
#     y += (p0 + p1 * x + p2 * x**2) * (x <= p4)
#     y += (p0 + p1 * p4 + p2 * p4**2) * (x > p4)
#     return y
def fitfunction(x, p0, p1, p2, p3):
    y = np.zeros(len(x))
    y += p0 * (x <= p1)
    y += (p2 * (x - p1) + p0 ) * (x <= p3) * (x > p1)
    y +=  (p2 * (p3 - p1) + p0 ) * (x > p3)
    return y

def fitfunction_real(x, p0, p1, p2, p3):
    if x <= p1:
        return p0
    if x <=p3:
        return p2 * (x - p1) + p0
    else:
        return p2 * (p3 - p1) + p0

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

        if (sum_error**0.5)/sum_weight < 0.35 and sum_weight > 7:
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

if __name__ == '__main__':
    ntag = 1
    debug = False
    cut = True
    sample_directory = ["../CxAOD31_01a/"]
    tag = "run2"
    rescale = True
    slopecorrection = True
    loadeasytree = False
    region = "mbbcr"

    if loadeasytree:
        mysamplembbcr1tag, t2 = get_sample(["mbbcr", "resolved", "1tag"])
        pickleit((mysamplembbcr1tag, t2), "pickle/mbbcr1tag")
        mysamplembbcr2tag, t2 = get_sample(["mbbcr", "resolved", "2tag"])
        pickleit((mysamplembbcr2tag, t2), "pickle/mbbcr2tag")
        mysamplesr1tag, t2 = get_sample(["sr", "resolved", "1tag"])
        pickleit((mysamplesr1tag, t2), "pickle/sr1tag")
        mysamplesr2tag, t2 = get_sample(["sr", "resolved", "2tag"])
        pickleit((mysamplesr2tag, t2), "pickle/sr2tag")
        exit(1)
    else:
        if ntag == 1 and region == "mbbcr":
            mysamplembbcr1tag, t2mbbcr1tag = unpickleit("pickle/mbbcr1tag")
            all_sample = mysamplembbcr1tag
        if ntag == 2 and region == "mbbcr":
            mysamplembbcr2tag, t2mbbcr2tag = unpickleit("pickle/mbbcr2tag")
            all_sample = mysamplembbcr2tag
        if ntag == 1 and region == "sr":
            mysamplesr1tag, t2sr1tag = unpickleit("pickle/sr1tag")
            all_sample = mysamplesr1tag
        if ntag == 2 and region == "sr":
            mysamplesr2tag, t2sr2tag = unpickleit("pickle/sr2tag")
            all_sample = mysamplesr2tag
    t2 = r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$"
    if tag == "a":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$"
    if tag == "d":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,43.6\:fb^{-1}}$"
    if tag == "e":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,58.5\:fb^{-1}}$"
    if tag == "run2":
        t2 = r"$\mathit{\sqrt{s}=13\:TeV,139\:fb^{-1}}$"
    # for i in range(len(all_sample)):
    #     all_sample[i].pth()
    
    all_sample_after = [each for each in all_sample]
    rescaledic = None
    if rescale:
        # rescaledic = loadnorm("C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/confignormonly.cfg",
        # "C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/GlobalFit_fitres_unconditionnal_mu0_normonly.txt")
        rescaledic = loadnorm("../fitconfig/config_m2000.cfg", "../fitconfig/GlobalFit_fitres_conditionnal_mu1.txt")
    if slopecorrection:
        p1s = get_slope_correction("output/slopefit_zlf/" + "pTH-mbbcut-"+str(1)+"tagpolyfitresult.csv")


    if rescale:
        print("Performing rescale...")
        for i in range(len(all_sample_after)):
            for each_key in rescaledic.keys():
                if 'ALL' in rescaledic[each_key]:
                    factor = rescaledic[each_key]['ALL'] + 1
                    mask = all_sample_after[i].mata["Sample"] == zlib.adler32(each_key.encode())
                    if True in mask:
                        all_sample_after[i].rescale(factor, mask)

    # to be made a function
    #-------------------------------------------------------
    poplist = []
    for i in range(len(all_sample_after)):
        if all_sample_after[i].alias == "Zlljet":
            poplist.append(i)
    temzjet = []
    for i in poplist:
        temzjet.append(all_sample_after[i])
    zhf, zlf = splitszjetsamples(temzjet)
    zlf.colour = "lightblue"
    all_sample_after = [each for i, each in enumerate(all_sample_after) if i not in poplist]
    all_sample_after.append(zhf)
    all_sample_after.append(zlf)
    #-------------------------------------------------------

    backup = [each for each in all_sample_after]
    all_sample_after_beforecorrection = copy.deepcopy(all_sample_after)

    if slopecorrection:
        print("Performing slope correction...")
        for i in range(len(all_sample_after)):
            if all_sample_after[i].alias == "Z+lf":
                all_sample_after[i].weight = all_sample_after[i].weight * (fitfunction(all_sample_after[i].data[b'ptH']/1000., p1s[0], p1s[1], p1s[2], p1s[3]))
    bins = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300]
    

    title3="mBBcr " + str(ntag) +" btags"
    direct = "output/t_make_plot/"
    name = "-mbbcut-" + str(ntag) +"tag"
    if rescale:
        direct = "output/t_make_plot_rescale_zlf/"
    if slopecorrection and rescale:
        direct = "output/t_make_plot_rescale_slopecorrection_zlf/"
    if slopecorrection and not rescale:
        direct = "output/t_make_plot_slopecorrection_zlf/"


    bins = [0, 25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000,  1050, 1100, 1150, 1200, 1250, 1300]
    bins = range(0,1400,20)
    bins = autobin_withdata(all_sample_after_beforecorrection, bins, alias="Z+lf", variable=b"pTV")
    print(bins)
    chi2, nod = stackplot(all_sample_after,b'pTV',bins,1000.,
            xlabel=r"$p_{TV}[GeV]$", title3=title3, filename=direct + "pTV" + name, print_height=True,
            title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpljets=True, chi2=True)
    print("pTV", chi2, nod)
    bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1150, 1350, 1550, 1800]
    chi2, nod = stackplot(all_sample_after,b'mVH',bins,1000.,
            xlabel=r"$m_{VH}[GeV]$", title3=title3, filename=direct + "mVH" + name, print_height=True,
            title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpljets=True, chi2=True)
    print("mVH", chi2, nod)
    bins = range(20, 200, 1)
    stackplot(all_sample_after,b'mBBres',bins,1000.,
            xlabel=r"$m_{BB}[GeV]$", title3=title3, filename=direct + "mBB" + name, print_height=True,
            title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=False, printzpljets=True, chi2=True)
    bins = range(0,1400,20)
    bins = autobin_withdata(all_sample_after_beforecorrection, bins, alias="Z+lf", variable=b'ptH')
    print(bins)
    chi2, nod = stackplot(all_sample_after,b'ptH',bins,1000.,
            xlabel=r"$p_{TH}[GeV]$", title3=title3, filename=direct + "pTH" + name, print_height=True,
            title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpljets=True, chi2=True)



    if not slopecorrection or True:
        all_sample_after1 = copy.deepcopy(backup)
        all_sample_after2 = copy.deepcopy(backup)

        for i in range(len(all_sample_after)):
            all_sample_after1[i].cut(cut_lowmbb)
            all_sample_after2[i].cut(cut_highmbb)
        all_sample_after1_beforecorrection = copy.deepcopy(all_sample_after1)
        all_sample_after2_beforecorrection = copy.deepcopy(all_sample_after2)
        if slopecorrection:
            p1s = get_slope_correction("output/slopefit_zlf/" + "pTH-highmbbcut-"+str(1)+"tagpolyfitresult.csv")
            for i in range(len(all_sample_after1)):
                if all_sample_after[i].alias == "Z+lf":
                    all_sample_after[i].weight = all_sample_after[i].weight * (fitfunction(all_sample_after[i].data[b'ptH']/1000., p1s[0], p1s[1], p1s[2], p1s[3]))
            p1s = get_slope_correction("output/slopefit_zlf/" + "pTH-lowmbbcut-"+str(1)+"tagpolyfitresult.csv")
            for i in range(len(all_sample_after2)):
                if all_sample_after[i].alias == "Z+lf":
                    all_sample_after[i].weight = all_sample_after[i].weight * (fitfunction(all_sample_after[i].data[b'ptH']/1000., p1s[0], p1s[1], p1s[2], p1s[3]))
        title3="lowmBBcr " + str(ntag) +" btags"
        name = "-lowmbbcut-" + str(ntag) +"tag"
        bins = range(0,1400,20)
        bins = autobin_withdata(all_sample_after1_beforecorrection, bins, alias="Z+lf", variable=b"pTV")
        print(bins)
        chi2, nod = stackplot(all_sample_after1,b'pTV',bins,1000.,
                xlabel=r"$p_{TV}[GeV]$", title3=title3, filename=direct + "pTV" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpljets=True, chi2=True)
        print("pTV", chi2, nod)
        bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1150, 1350, 1550, 1800]
        chi2, nod = stackplot(all_sample_after1,b'mVH',bins,1000.,
                xlabel=r"$m_{VH}[GeV]$", title3=title3, filename=direct + "mVH" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpljets=True, chi2=True)
        print("mVH", chi2, nod)
        bins = range(20, 200, 1)
        stackplot(all_sample_after1,b'mBBres',bins,1000.,
                xlabel=r"$m_{BB}[GeV]$", title3=title3, filename=direct + "mBB" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=False, printzpljets=True, chi2=True)
        bins = range(0,1400,20)
        bins = autobin_withdata(all_sample_after1_beforecorrection, bins, alias="Z+lf", variable=b'ptH')
        print(bins)
        chi2, nod = stackplot(all_sample_after1,b'ptH',bins,1000.,
                xlabel=r"$p_{TH}[GeV]$", title3=title3, filename=direct + "pTH" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpljets=True, chi2=True)



        title3="highmBBcr " + str(ntag) +" btags"
        name = "-highmbbcut-" + str(ntag) +"tag"
        bins = range(0,1400,20)
        bins = autobin_withdata(all_sample_after2_beforecorrection, bins, alias="Z+lf", variable=b"pTV")
        print(bins)
        chi2, nod = stackplot(all_sample_after2,b'pTV',bins,1000.,
                xlabel=r"$p_{TV}[GeV]$", title3=title3, filename=direct + "pTV" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpljets=True, chi2=True)
        print("pTV", chi2, nod)
        bins = [50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 1000, 1150, 1350, 1550, 1800]
        chi2, nod = stackplot(all_sample_after2,b'mVH',bins,1000.,
                xlabel=r"$m_{VH}[GeV]$", title3=title3, filename=direct + "mVH" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpljets=True, chi2=True)
        print("mVH", chi2, nod)
        bins = range(20, 200, 1)
        stackplot(all_sample_after2,b'mBBres',bins,1000.,
                xlabel=r"$m_{BB}[GeV]$", title3=title3, filename=direct + "mBB" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=False, printzpljets=True, chi2=True)
        bins = range(0,1400,20)
        bins = autobin_withdata(all_sample_after2_beforecorrection, bins, alias="Z+lf", variable=b'ptH')
        print(bins)
        chi2, nod = stackplot(all_sample_after2,b'ptH',bins,1000.,
                xlabel=r"$p_{TH}[GeV]$", title3=title3, filename=direct + "pTH" + name, print_height=True,
                title2=t2,auto_colour=False, limit_y = 0.5, upper_y=2.0, log_y=True, printzpljets=True, chi2=True)