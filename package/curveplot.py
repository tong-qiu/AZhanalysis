import matplotlib as mpl
mpl.use('Agg')
import matplotlib
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np
import pickle
from matplotlib import rc, rcParams
from events import *
import matplotlib.font_manager as font_manager
import matplotlib.patches as mpatches

def histplot_raw(datas, bins, labels, weights=None, removenorm=False, scale=1., **kwargs):
    settings = {
        "xlabel" : r"$m_{Vh}[GeV]$",
        "ylabel": 'Number of Events',
        "title1": r"ATLAS",# \newline Ptl next-leading, full cuts, 2 b-tags $",
        "title1_1": r"Internal",
        "title2": r"$\mathit{\sqrt{s}=13\:TeV,139\:fb^{-1}}$",# Ptl next-leading, full cuts, 2 b-tags $",
        #"title3": r"$\mathbf{2\;lep.,2\;b-tag}$",
        "title3": "2 lep., 2 b-tag",
        "filename": "deltatest2",
        "log_y":False,
        "norm":False,
        "upper_y": 1.5, 
        }
    for each_key in kwargs.items():
        settings[each_key[0]] = kwargs[each_key[0]]
    

    if weights is None:
        weights = []
        for each in datas:
            weights.append(np.ones(len(each)))
    
    if removenorm:
        for i in range(len(weights)):
            weights[i] = np.array(weights[i]) / np.sum(weights[i])
    sigmas = []
    weight_in_binses = []
    for i in range(len(datas)):
        event_location = np.digitize(datas[i]/scale, bins)
        sigma2 = []
        weight_in_bins = []
        for j in range(np.size(bins) - 1):
            bin_weight = weights[i][np.where(event_location == j+1)[0]]
            sigma2.append(np.sum(bin_weight**2.))
            weight_in_bins.append(np.sum(bin_weight))
        sigmas.append(np.array(sigma2)**0.5)
        weight_in_binses.append(np.array(weight_in_bins))

    colors = ['b', 'g', 'r', 'c', 'm', 'y']
    fig, ax = plt.subplots(figsize=(10,8))
    ax.hist(np.array(datas)/scale, bins, histtype='step', fill=False, color=colors[0:len(datas)], weights=weights)
    bins = np.array(bins)
    for i in range(len(datas)):
        bin_centre = (bins[0:-1] + bins[1:])/2
        ax.errorbar(bin_centre, weight_in_binses[i], xerr=0.0001, yerr=sigmas[i], fmt='.', color=colors[i], label=str(labels[i]))
    ax.legend(loc='upper right',prop={'size': 20}, frameon=False)
    
    ymin, ymax = ax.get_ylim()
    ax.set_ylim([0,ymax* settings["upper_y"]])
    ax.text(0.05, 1.55 / 1.7, settings['title1'], fontsize=25, transform=ax.transAxes, style='italic', fontweight='bold')
    ax.text(0.227, 1.55/ 1.7, settings['title1_1'], fontsize=25, transform=ax.transAxes)
    ax.text(0.05, 1.40 / 1.7, settings['title2'], fontsize=23, transform=ax.transAxes, style='italic', fontweight='bold')
    ax.text(0.05, 1.26 / 1.7, settings['title3'], fontsize=18, weight='bold', style='italic', transform=ax.transAxes)
    ax.set_ylabel(settings['ylabel'], fontsize=20)
    ax.set_xlabel(settings['xlabel'], fontsize=20)
    if settings['log_y']:
        ax.set_yscale('log')
        ax.set_ylim([0.1, 10**(math.log10(ymax) * settings["upper_y"])])
        ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10,numticks=100))
        ax.minorticks_on()

    fig.savefig(settings['filename'] + '.pdf', bbox_inches='tight', pad_inches = 0.25)

def curveplot(x_list, y_list, error_list=[], labels=None, **kwargs):
    # rcParams['mathtext.fontset'] = 'custom'
    # rcParams['mathtext.it'] = 'DejaVu Sans:italic'
    # rcParams['mathtext.bf'] = 'DejaVu Sans:italic:bold'
    # default label
    settings = {
        "xlabel" : r"$m_{Vh}[GeV]$",
        "ylabel": 'Number of Events',
        "title1": r"ATLAS",# \newline Ptl next-leading, full cuts, 2 b-tags $",
        "title1_1": r"Internal",
        "title2": r"$\sqrt{s}=13\:TeV,139\:fb^{-1}$",# Ptl next-leading, full cuts, 2 b-tags $",
        #"title3": r"$\mathbf{2\;lep.,2\;b-tag}$",
        "title3": "2 lep., 2 b-tag",
        "filename": "deltatest2",
        "log_y":False,
        "ylimit":[0.5,1.2],
        "xlimit":None,
        "yshift":1.3,
        "xshift":0,
        "verticleline":None,
        "verticlelinetext": "",
        "horizontalline":None,
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
        if len(y_list) == 1:
            cur = plt.plot(each_x, each_y, "k.", label=label)
            plt.plot(each_x, each_y, "k")
        else:
            cur = plt.plot(each_x, each_y, label=label)
        all_curve.append(cur)
        i += 1
    if labels is not None:
        plt.legend(loc='upper right', prop={'size': 20})
    ax = plt.gca()
    shift = settings["yshift"]
    xshift = settings["xshift"]
    ax.text(0.05 + xshift, (1.55 - shift) / 1.7, settings['title1'], fontsize=25, transform=ax.transAxes, style='italic', fontweight='bold')
    ax.text(0.227 + xshift, (1.55 - shift) / 1.7, settings['title1_1'], fontsize=25, transform=ax.transAxes)
    ax.text(0.05 + xshift, (1.40 - shift) / 1.7, settings['title2'], fontsize=23, transform=ax.transAxes, style='italic', fontweight='bold')
    ax.text(0.05 + xshift, (1.26  - shift) / 1.7, settings['title3'], fontsize=18, weight='bold', style='italic', transform=ax.transAxes)
    #ax1.text(0.05, 1.12 / 1.7, settings["title4"], fontsize=18, weight='bold', style='italic', transform=ax1.transAxes)
    
    if len(error_list)!=0:
        for each_x, each_y, each_error in zip(x_list, y_list, error_list):
            plt.errorbar(each_x, each_y, yerr=each_error,  fmt='.')
    plt.tick_params(labelsize=16)
    plt.tick_params(labelsize=16)
    axes = plt.gca()
    axes.set_ylim(settings["ylimit"])
    if settings["xlimit"] is not None:
        axes.set_xlim(settings["xlimit"])
    if settings["verticleline"] is not None:
        plt.plot([settings["verticleline"], settings["verticleline"]], [0,100000], ':', color='silver')
        plt.annotate(settings["verticlelinetext"], xy=(settings["verticleline"],10), xycoords='data', color='grey')
    if settings["horizontalline"] is not None:
        plt.plot([0,100000], [settings["horizontalline"], settings["horizontalline"]], ':', color='silver')
    plt.ylabel(settings['ylabel'], fontsize=20)
    plt.xlabel(settings['xlabel'], fontsize=20)
    plt.savefig(settings['filename'] + '.pdf', bbox_inches='tight', pad_inches = 0.25)
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
    plt.clf()

def histplot_withsub(data_lists, varible_to_plot, bins, labels = None, scales=1., removenorm = None, **kwargs,):
    # data_list = [[sample1, sample2, ...],[sample3, sample4, ...], ...]
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
        "central":"none",
        "upper_y": 1.5, 
        "do_errorbar": False
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
            height_tem, sigma2_tem = each_sample.binned_weight_variation(varible_to_plot, bins, scales)
            if len(height) == 0:
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
    
    if removenorm is not None:
        for i in range(len(all_label)):
            if all_label[i] == removenorm:
                sumweightnorm = sum(all_height[i])
                break
        for i in range(len(all_label)):
            if all_label == "data":
                continue
            all_height[i] = all_height[i] * sumweightnorm/sum(all_height[i])
            all_sigma2[i] = all_sigma2[i] * sumweightnorm/sum(all_height[i])
    handles = []
    for each_height, each_sigma2, each_label, each_width, each_color in zip(all_height, all_sigma2, all_label, bin_width, all_color):
        ax1.hist(bin_centre, bins, weights = each_height, label=each_label, density=settings["norm"], histtype=u'step', color=each_color)
        handles.append(mpl.lines.Line2D([], [], c=each_color, label=each_label))

        #ax1.errorbar(bin_centre, each_height, yerr=each_sigma2*0, label=each_label, fmt='.', color=each_color) # xerr=each_width/2.
    # lower plot
    data_height = None
    data_error = None
    for each_height, each_sigma2, each_label, each_width in zip(all_height, all_sigma2, all_label, bin_width):
        #if "data" in each_label:
        if settings["central"] in each_label:
            data_height = np.array(each_height)
            data_error = np.array(each_sigma2)**0.5
            break
    for each_height, each_sigma2, each_label, each_width, each_color in zip(all_height, all_sigma2, all_label, bin_width, all_color):
        #if "data" in each_label:
        if settings["central"] in each_label:
            continue
        if data_height is not None:
            new_each_height = np.array(each_height)/data_height
            new_each_sigma2 = each_sigma2**0.5/np.array(each_height)
        else:
            print("histplot_withsub: WARNING Cannot find central")
            new_each_sigma2 = each_sigma2**0.5
            new_each_height = np.array(each_height)
        new_each_height[np.isnan(new_each_height)] = -1
        # if "sys" in each_label:
        #     print("here")
        #     new_each_height = new_each_height + 0.1
        if settings['do_errorbar']:
            ax2.errorbar(bin_centre, new_each_height-1, yerr=new_each_sigma2, label=each_label, fmt='.', color=each_color)
        ax2.hist(bin_centre, bins, weights=new_each_height-1, label=each_label, histtype=u'step', color=each_color)
    ax2.plot([bins[0], bins[np.size(bins)-1]], [0, 0], linestyle='--', color='k')

    i = -1
    if settings['do_errorbar']:
        for each_x, each_error, each_y_mc in zip(bin_centre, data_error, data_height):
            i += 1
            # do not plot bins without event?
            if each_y_mc == 0:
                continue
            ax2.add_patch(
                mpl.patches.Rectangle(
                    (each_x - (bins[i+1]-bins[i]) / 2., - each_error/each_y_mc), # x, y
                    bins[i+1]-bins[i],        # width
                    each_error/each_y_mc*2,        # height
                    color='black', alpha=0.5,
                    hatch='/////', fill=False, linewidth=0,
                    ))

    ymin, ymax = ax1.get_ylim()
    ax1.set_ylim([0,ymax* settings["upper_y"]])
    ax1.text(0.05, 1.55 / 1.7, settings['title1'], fontsize=25, transform=ax1.transAxes, style='italic', fontweight='bold')
    ax1.text(0.227, 1.55/ 1.7, settings['title1_1'], fontsize=25, transform=ax1.transAxes)
    ax1.text(0.05, 1.40 / 1.7, settings['title2'], fontsize=23, transform=ax1.transAxes, style='italic', fontweight='bold')
    ax1.text(0.05, 1.26 / 1.7, settings['title3'], fontsize=18, weight='bold', style='italic', transform=ax1.transAxes)
    #ax2.set_ylim([0.99, 1.01])
    if len(all_height) > 1:
        ax1.legend(loc='upper right',prop={'size': 20}, frameon=False, handles=handles)
    ax1.get_xaxis().set_ticks([])
    ax1.tick_params(labelsize=16)
    ax1.tick_params(labelsize=16)
    ax1.set_ylabel(settings['ylabel'], fontsize=20)
    ax2.set_ylabel('mc/' + settings["central"] + "-1", fontsize=20)
    ax2.set_xlabel(settings['xlabel'], fontsize=20)
    font = font_manager.FontProperties(weight='bold',
                                       style='normal', size=18)
    plt.savefig(settings['filename'] + '.pdf', bbox_inches='tight', pad_inches = 0.1)
    plt.close(fig)