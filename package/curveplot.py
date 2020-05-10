import matplotlib as mpl
mpl.use('Agg')
import math
import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib import rc, rcParams
from events import *
import matplotlib.font_manager as font_manager
import matplotlib.patches as mpatches

def curveplot(x_list, y_list, error_list=[], labels=None, **kwargs):
    rcParams['mathtext.fontset'] = 'custom'
    rcParams['mathtext.it'] = 'DejaVu Sans:italic'
    rcParams['mathtext.bf'] = 'DejaVu Sans:italic:bold'
    # default label
    settings = {
        "xlabel" : r"$m_{Vh}[GeV]$",
        "ylabel": 'Number of Events',
        "title1": r"$\mathbf{ATLAS}$",# \newline Ptl next-leading, full cuts, 2 b-tags $",
        "title1_1": r"$\mathit{Internal}$",
        "title2": r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$",# Ptl next-leading, full cuts, 2 b-tags $",
        #"title3": r"$\mathbf{2\;lep.,2\;b-tag}$",
        "title3": "2 lep., 2 b-tag",
        "filename": "deltatest2",
        "log_y":False
        }
    for each_key in kwargs.items():
        settings[each_key[0]] = kwargs[each_key[0]]
    plt.figure(figsize=(10, 8))
    i = 0
    all_curve = []
    for each_x, each_y in zip(x_list, y_list):
        label = "none"
        if labels is not None:
            label = labels[i]
        cur = plt.plot(each_x, each_y, label=label)
        all_curve.append(cur)
        i += 1
    if labels is not None:
        plt.legend(loc='upper right', prop={'size': 20})
    ax = plt.gca()
    shift = 1.3
    ax.text(0.05, (1.55 - shift) / 1.7, settings['title1'], fontsize=25, transform=ax.transAxes)
    ax.text(0.227, (1.55 - shift) / 1.7, settings['title1_1'], fontsize=21, transform=ax.transAxes)
    ax.text(0.05, (1.40 - shift) / 1.7, settings['title2'], fontsize=23, transform=ax.transAxes)
    #ax1.text(0.05, 1.26 / 1.7, settings['title3'], fontsize=18, weight='bold', style='italic', transform=ax1.transAxes)
    #ax1.text(0.05, 1.12 / 1.7, settings["title4"], fontsize=18, weight='bold', style='italic', transform=ax1.transAxes)
    
    if len(error_list)!=0:
        for each_x, each_y, each_error in zip(x_list, y_list, error_list):
            plt.errorbar(each_x, each_y, yerr=each_error,  fmt='.')
    plt.tick_params(labelsize=16)
    plt.tick_params(labelsize=16)
    axes = plt.gca()
    axes.set_ylim([0.5,1.2])
    plt.ylabel(settings['ylabel'], fontsize=20)
    plt.xlabel(settings['xlabel'], fontsize=20)
    plt.savefig(settings['filename'] + '.pdf', bbox_inches='tight', pad_inches = 0)
    plt.show()

def histplot(data_lists, varible_to_plot, bins, labels = None, scales=1., **kwargs):
    # data_list = [[sample1, sample2, ...],[sample3, sample4, ...], ...]
    rcParams['mathtext.fontset'] = 'custom'
    #rcParams["font.family"] = "Times New Roman"
    rcParams['mathtext.it'] = 'DejaVu Sans:italic'
    rcParams['mathtext.bf'] = 'DejaVu Sans:italic:bold'
    #rcParams['axes.unicode_minus'] = True
    # default label
    settings = {
        "xlabel" : r"$m_{Vh}[GeV]$",
        "ylabel": 'Number of Events',
        "title1": r"$\mathbf{ATLAS}$",# \newline Ptl next-leading, full cuts, 2 b-tags $",
        "title1_1": r"$\mathit{Internal}$",
        "title2": r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$",# Ptl next-leading, full cuts, 2 b-tags $",
        #"title3": r"$\mathbf{2\;lep.,2\;b-tag}$",
        "title3": "2 lep., 2 b-tag",
        "filename": "deltatest2",
        "log_y":False,
        "norm":False,
        }
    for each_key in kwargs.items():
        settings[each_key[0]] = kwargs[each_key[0]]
    plt.figure(figsize=(10, 8))
    all_height = []
    all_sigma2 = []
    all_label = labels
    for each_sample_list in data_lists:
        height = []
        sigma2 = []
        for each_sample in each_sample_list:
            height_tem, sigma2_tem = each_sample.binned_weight_variation(varible_to_plot,bins,scales)
            if not height:
                height = np.array(height_tem)
                sigma2 = np.array(sigma2_tem)
            else:
                height += np.array(height_tem)
                sigma2 += np.array(sigma2_tem)
        all_height.append(height)
        all_sigma2.append(sigma2)
        if labels is None:
            if all_label is None:
                all_label = []
            all_label.append(None)
    bin_centre = []
    bin_width = []
    for i in range(len(bins)-1):
        bin_centre.append((bins[i]+bins[i+1])/2.)
        bin_width.append(bins[i+1]-bins[i])
    for each_height, each_sigma2, each_label, each_width in zip(all_height, all_sigma2, all_label, bin_width):
        plt.hist(bin_centre, bins, weights = each_height, label=each_label, density=settings["norm"], histtype=u'step')
        #print("here")
        #plt.errorbar(bin_centre, each_height, yerr=each_sigma2**0.5, label=each_label, fmt='.') # xerr=each_width/2.
    if len(all_height) > 1:
        plt.legend(loc='upper right',prop={'size': 25})

    plt.tick_params(labelsize=16)
    plt.tick_params(labelsize=16)
    plt.ylabel(settings['ylabel'], fontsize=20)
    plt.xlabel(settings['xlabel'], fontsize=20)
    font = font_manager.FontProperties(weight='bold',
                                       style='normal', size=18)
    plt.savefig(settings['filename'] + '.pdf')
    plt.show()

def histplot_withsub(data_lists, varible_to_plot, bins, labels = None, scales=1., **kwargs):
    # data_list = [[sample1, sample2, ...],[sample3, sample4, ...], ...]
    rcParams['mathtext.fontset'] = 'custom'
    #rcParams["font.family"] = "Times New Roman"
    rcParams['mathtext.it'] = 'DejaVu Sans:italic'
    rcParams['mathtext.bf'] = 'DejaVu Sans:italic:bold'
    #rcParams['axes.unicode_minus'] = True
    # default label
    settings = {
        "xlabel" : r"$m_{Vh}[GeV]$",
        "ylabel": 'Number of Events',
        "title1": r"$\mathbf{ATLAS}$",# \newline Ptl next-leading, full cuts, 2 b-tags $",
        "title1_1": r"$\mathit{Internal}$",
        "title2": r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$",# Ptl next-leading, full cuts, 2 b-tags $",
        #"title3": r"$\mathbf{2\;lep.,2\;b-tag}$",
        "title3": "2 lep., 2 b-tag",
        "filename": "deltatest2",
        "log_y":False,
        "norm":False,
        "central":"data",
        }
    for each_key in kwargs.items():
        settings[each_key[0]] = kwargs[each_key[0]]
    #plt.figure(figsize=(10, 8))
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios':[3, 1]}, figsize=(10, 10))
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.06)
    all_height = []
    all_sigma2 = []
    all_color = []
    all_label = labels
    for each_sample_list in data_lists:
        height = []
        sigma2 = []
        color = None
        for each_sample in each_sample_list:
            color = each_sample.colour
            height_tem, sigma2_tem = each_sample.binned_weight_variation(varible_to_plot,bins,scales)
            if not height:
                height = np.array(height_tem)
                sigma2 = np.array(sigma2_tem)
            else:
                height += np.array(height_tem)
                sigma2 += np.array(sigma2_tem)
        all_height.append(height)
        all_sigma2.append(sigma2)
        all_color.append(color)
        if labels is None:
            if all_label is None:
                all_label = []
            all_label.append(None)
    bin_centre = []
    bin_width = []
    for i in range(len(bins)-1):
        bin_centre.append((bins[i]+bins[i+1])/2.)
        bin_width.append(bins[i+1]-bins[i])
    for each_height, each_sigma2, each_label, each_width, each_color in zip(all_height, all_sigma2, all_label, bin_width, all_color):
        ax1.hist(bin_centre, bins, weights = each_height, label=each_label, density=settings["norm"], histtype=u'step', color=each_color)
        #plt.errorbar(bin_centre, each_height, yerr=each_sigma2**0.5, label=each_label, fmt='.') # xerr=each_width/2.
    # lower plot
    data_height = None
    for each_height, each_sigma2, each_label, each_width in zip(all_height, all_sigma2, all_label, bin_width):
        #if "data" in each_label:
        if settings["central"] in each_label:
            data_height = np.array(each_height)
            break
    for each_height, each_sigma2, each_label, each_width, each_color in zip(all_height, all_sigma2, all_label, bin_width, all_color):
        #if "data" in each_label:
        if settings["central"] in each_label:
            continue

        new_each_height = np.array(each_height)/data_height
        new_each_height[np.isnan(new_each_height)] = -1
        # if "sys" in each_label:
        #     print("here")
        #     new_each_height = new_each_height + 0.1
        ax2.hist(bin_centre, bins, weights=new_each_height-1, label=each_label, histtype=u'step', color=each_color)
    ax2.plot([bins[0], bins[np.size(bins)-1]], [0, 0], linestyle='--', color='k')
    #ax2.set_ylim([0.99, 1.01])
    if len(all_height) > 1:
        ax1.legend(loc='upper right',prop={'size': 25})
    ax1.get_xaxis().set_ticks([])
    ax1.tick_params(labelsize=16)
    ax1.tick_params(labelsize=16)
    ax1.set_ylabel(settings['ylabel'], fontsize=20)
    ax2.set_ylabel('mc/' + settings["central"] + "-1", fontsize=20)
    ax2.set_xlabel(settings['xlabel'], fontsize=20)
    font = font_manager.FontProperties(weight='bold',
                                       style='normal', size=18)
    plt.savefig(settings['filename'] + '.pdf')
    plt.close(fig)