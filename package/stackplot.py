import matplotlib
import math
import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib import rc, rcParams
from events import *
import matplotlib.font_manager as font_manager
import matplotlib.patches as mpatches

def stackplot(data_list, varible_to_plot, bins, scales=1., **kwargs):
    #rc('text', usetex=True)
    #rc('font',family='Times New Roman')
    rcParams['mathtext.fontset'] = 'custom'
    #rcParams["font.family"] = "Times New Roman"
    rcParams['mathtext.it'] = 'DejaVu Sans:italic'
    rcParams['mathtext.bf'] = 'DejaVu Sans:italic:bold'
    #rcParams['axes.unicode_minus'] = True
    # default label
    settings = {
        "xlabel" : r"$m_{Vh}[GeV]$",
        "ylabelup": 'Number of Events',
        "ylabeldown": "Data/Simulation",
        "title1": r"$\mathbf{ATLAS}$",# \newline Ptl next-leading, full cuts, 2 b-tags $",
        "title1_1": r"$\mathit{Internal}$",
        #"title1_1": r"$\mathit{Working\:in\:progress}$",
        "title2": r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$",# Ptl next-leading, full cuts, 2 b-tags $",
        #"title3": r"$\mathbf{2\;lep.,2\;b-tag}$",
        "title3": "2 lep., 2 b-tag",
        "title4": "",
        "filename": "deltatest2",
        "auto_colour": True,
        "print_height": False,
        "limit_y": 0,
        "pickle": False,
        "upper_y": 1.7,
        "ncol":1,
        "log_y":False,
        "sys":False,
        "blind":False,
        }
    for each_key in kwargs.items():
        settings[each_key[0]] = kwargs[each_key[0]]

    # load data
    varibles_content = []
    weight = []
    alias = []
    sigma2 = []
    colours = []
    data_content = []
    weight_in_bin = []
    sys2 = []
    data_weight = []

    new_data_list = []
    for each in data_list:
        ifpass = False
        for i, each_new in enumerate(new_data_list):
            if each_new.alias == each.alias:
                new_data_list[i] = each_new + each
                print("Info: Merge ", each_new.alias)
                ifpass = True
                break
        if not ifpass:
            new_data_list.append(each)
    data_list = new_data_list

    def nevent(asample):
        thesum = sum(asample.weight)
        return thesum#len(asample.data[varible_to_plot])
    data_list.sort(key=nevent)
    for each in data_list:
        print(nevent(each))
        if nevent(each) < 0:
            print(each.alias + " is empty.")
            continue
        if "data" not in each.alias:#np.mean(each.weight) != 1:# is not None:
            if each.alias in alias and False:
                i = alias.index(each.alias)
                print("Info: Merge ", alias[i])
                weight[i] = np.append(weight[i], each.weight)
                varibles_content[i] = np.append(varibles_content[i], each.data[varible_to_plot]/scales)
                sigma2[i] = [sum(x) for x in zip(sigma2[i], each.variation(varible_to_plot, bins, scales))]#sigma2[i] + each.variation(varible_to_plot, bins, scales)
                weight_in_bin[i] = [sum(x) for x in zip(weight_in_bin[i], each.binned_weight(varible_to_plot, bins, scales))]#each.binned_weight(varible_to_plot, bins, scales)
                if settings["sys"] and len(each.systematics(varible_to_plot, bins, scales))>0:
                    sys2[i] = [sum(x) for x in zip(sys2[i], each.systematics(varible_to_plot, bins, scales) )]#np.append(sys2[i], each.systematics(varible_to_plot, bins, scales))
                    #raise ValueError("Unfinished work")
                continue
            weight.append(each.weight)
            varibles_content.append(each.data[varible_to_plot]/scales)
            alias.append(each.alias)
            sigma2.append(each.variation(varible_to_plot, bins, scales))
            weight_in_bin.append(each.binned_weight(varible_to_plot, bins, scales))
            colours.append(each.colour)
            if settings["sys"] and len(each.systematics(varible_to_plot, bins, scales))>0:
                sys2.append(each.systematics(varible_to_plot, bins, scales))
        else:
            if settings["blind"]:
                data_weight = [0 for each in range(len(each.data[varible_to_plot]))]
            else:
                data_weight = [1 for each in range(len(each.data[varible_to_plot]))]
            data_content = each.data[varible_to_plot]/scales
    # calculate uncertainty
    sigma2 = np.transpose(sigma2)
    if settings["sys"]:
        sys2 = np.transpose(sys2)
        error_mc = np.sqrt(np.sum(sigma2, (1))+np.sum(sys2, (1)))
        #error_mc = np.sqrt(np.sum(sys2, (1)))
    else:
        error_mc = np.sqrt(np.sum(sigma2, (1)))
    # overall plot setting
    fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios':[3, 1]}, figsize=(10, 10))
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.06)
    ax1.set_xlim([bins[0], bins[np.size(bins)-1]])
    ax2.set_xlim([bins[0], bins[np.size(bins)-1]])
    #plt.set_ylim([ymin,ymax])

    # upper plot
    bin_heights, bin_edges = np.histogram(data_content, bins=bins, weights=data_weight)
    bincenters = 0.5*(bin_edges[1:]+bin_edges[:-1])
    mean_std = np.sqrt(bin_heights)

    # #add text to upper plot
    # ax1.set_ylim([0, max(bin_heights)* settings["upper_y"]])
    # ax1.text(0.05, 1.55 / 1.7, settings['title1'], fontsize=25, transform=ax1.transAxes)
    # ax1.text(0.227, 1.55/ 1.7, settings['title1_1'], fontsize=21, transform=ax1.transAxes)
    # ax1.text(0.05, 1.40 / 1.7, settings['title2'], fontsize=23, transform=ax1.transAxes)
    # ax1.text(0.05, 1.26 / 1.7, settings['title3'], fontsize=18, weight='bold', style='italic', transform=ax1.transAxes)
    # ax1.text(0.05, 1.12 / 1.7, "PTl2 < 23 GeV", fontsize=18, weight='bold', style='italic', transform=ax1.transAxes)
    # !!! change plot title
    #ax1.text(min(bincenters) + (max(bincenters)-min(bincenters)) * 0.05, max(bin_heights)*1.55 / 1.7  * settings["upper_y"],
    #         settings['title1'], fontsize=25)
    #ax1.text(bins[0] + (bins[len(bins)-1] - bins[0]) * 0.05, max(bin_heights)*1.55 / 1.7  * settings["upper_y"],
    #         settings['title1'], fontsize=25)
    #ax1.text(min(bincenters) + (max(bincenters)-min(bincenters)) * 0.235, max(bin_heights)*1.55/ 1.7  * settings["upper_y"],
    #         settings['title1_1'], fontsize=21)
    #ax1.text(bins[0] + (bins[len(bins)-1] - bins[0]) * 0.227, max(bin_heights)*1.55/ 1.7  * settings["upper_y"],
    #         settings['title1_1'], fontsize=21)
    #ax1.text(bins[0] + (bins[len(bins)-1] - bins[0]) * 0.05, max(bin_heights)*1.40 / 1.7  * settings["upper_y"],
    #         settings['title2'], fontsize=23)
    #ax1.text(bins[0] + (bins[len(bins)-1] - bins[0]) * 0.05, max(bin_heights)*1.26 / 1.7  * settings["upper_y"],
    #         settings['title3'], fontsize=18, weight='bold', style='italic')
    
    # plot data
    error_patch = ax1.errorbar(bincenters, bin_heights, color='k',
                               linestyle="", yerr=mean_std, fmt='o', markersize='8')
    #ax1.plot(bincenters, bin_heights, marker=".", linestyle="", color="k")

    # plot MC
    if settings['auto_colour']:
        y_mcs = ax1.hist(varibles_content, bins, weights=weight,
                         stacked=True, label=alias)
    else:
        y_mcs = ax1.hist(varibles_content, bins, weights=weight,
                         stacked=True, label=alias, color=colours)
    ax1.get_xaxis().set_ticks([])
    font = font_manager.FontProperties(weight='bold',
                                       style='normal', size=18)

    #add text to upper plot
    #y_mcs = np.transpose(y_mcs)
    #y_mc_max = np.max(np.sum(y_mcs[0], axis=0))
    y_mc_max = np.max(y_mcs[0][int(np.size(y_mcs[0])/np.size(y_mcs[0][0]))-1])
    if type(y_mcs[0][int(np.size(y_mcs[0])/np.size(y_mcs[0][0]))-1]) is not np.ndarray:
        y_mc_max = np.max(y_mcs[0])

    ax1.set_ylim([0,max([y_mc_max, max(bin_heights)])* settings["upper_y"]])
    ax1.text(0.05, 1.55 / 1.7, settings['title1'], fontsize=25, transform=ax1.transAxes)
    ax1.text(0.227, 1.55/ 1.7, settings['title1_1'], fontsize=21, transform=ax1.transAxes)
    ax1.text(0.05, 1.40 / 1.7, settings['title2'], fontsize=23, transform=ax1.transAxes)
    ax1.text(0.05, 1.26 / 1.7, settings['title3'], fontsize=18, weight='bold', style='italic', transform=ax1.transAxes)
    ax1.text(0.05, 1.12 / 1.7, settings["title4"], fontsize=18, weight='bold', style='italic', transform=ax1.transAxes)

    # plot the error bars of MC
    y_mc = y_mcs[0][int(np.size(y_mcs[0])/np.size(y_mcs[0][0]))-1]
    if type(y_mc) is not np.ndarray:
        y_mc = y_mcs[0]
    i = -1
    for each_x, each_error, each_y_mc in zip(bincenters, error_mc, y_mc):
        i += 1
        # do not plot bins without event?
        if each_y_mc == 0:
            print("Warning: bin at " + str(each_x) +  " has no MC events")
            continue
        ax1.add_patch(
            matplotlib.patches.Rectangle(
                (each_x - (bins[i+1]-bins[i]) / 2., each_y_mc - each_error), # x, y
                bins[i+1]-bins[i],        # width
                each_error *2,        # height
                color='black', alpha=0.5,
                hatch='/////', fill=False, linewidth=0,
                ))

    if settings["print_height"]:
        with open(settings["filename"]+".csv", "w") as f:
            f.write("data" + "," + "MC" + ",\n")
            for each_data, each_mc in zip(bin_heights, y_mc):
                f.write(str(each_data) + "," + str(each_mc) + ",\n")

    # add legend
    sys_patch = mpatches.Patch(color='black', hatch='/////', fill=False, linewidth=0)
    ax1.legend([error_patch, sys_patch] + y_mcs[2], ["data", "uncertainty"] + alias, prop=font, frameon=False, ncol=settings["ncol"])
    if settings["log_y"]:
        ax1.set_yscale('log')
        ax1.set_ylim([0.1, 10**(math.log10(max(bin_heights)) * settings["upper_y"])])
        #locmaj = matplotlib.ticker.LogLocator(base=10, subs=np.arange(2, 10))
        #ax1.yaxis.set_major_locator(locmaj)
        ax1.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10,numticks=100))
        ax1.minorticks_on()
        #locmin = matplotlib.ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * .1,
        #                                    numticks=100)

    # make lower plot
    y_mc = y_mcs[0][int(np.size(y_mcs[0])/np.size(y_mcs[0][0]))-1]
    if type(y_mc) is not np.ndarray:
        y_mc = y_mcs[0]
    # plot error bar of data
    # remove bins without MC events?
    error_bar_center = []
    error_bar_size = []
    error_bar_bincentre = []
    for each_y, each_y_mc, each_menstd, each_error_bar_bincentre in zip(bin_heights, y_mc, mean_std, bincenters):
        error_bar_center.append(each_y / each_y_mc)
        error_bar_size.append(each_menstd / each_y_mc)
        error_bar_bincentre.append(each_error_bar_bincentre)

    ax2.errorbar(error_bar_bincentre, error_bar_center, color='k', linestyle="",
                 yerr=error_bar_size, fmt='o', markersize='8')

    # plot line at the centre
    ax2.plot([bins[0], bins[np.size(bins)-1]], [1, 1], linestyle='--', color='k')

    # plot the error bars of MC
    i = -1
    for each_x, each_error, each_y_mc in zip(bincenters, error_mc, y_mc):
        i += 1
        # do not plot bins without event?
        if each_y_mc == 0:
            continue
        ax2.add_patch(
            matplotlib.patches.Rectangle(
                (each_x - (bins[i+1]-bins[i]) / 2., - each_error/each_y_mc + 1), # x, y
                bins[i+1]-bins[i],        # width
                each_error/each_y_mc*2,        # height
                color='black', alpha=0.5,
                hatch='/////', fill=False, linewidth=0,
                ))
    # limit the range of y
    if settings["limit_y"] > 0:
        ax2.set_ylim([1 - settings["limit_y"], 1 + settings["limit_y"]])
        error_bar_center = np.array(error_bar_center)
        error_bar_bincentre = np.array(error_bar_bincentre)
        error_bar_bincentre_up = error_bar_bincentre[error_bar_center > 1 + settings["limit_y"]]
        error_bar_bincentre_down = error_bar_bincentre[1 - settings["limit_y"] > error_bar_center]
        #error_bar_center_up = error_bar_center[error_bar_center > 1 + settings["limit_y"]]
        #error_bar_center_down = error_bar_center[error_bar_center < 1 - settings["limit_y"]]
        ax2.plot(error_bar_bincentre_up,np.linspace(1 + settings["limit_y"]*0.85, 1 + settings["limit_y"]*0.85, len(error_bar_bincentre_up)), 'k^', markersize='8')
        ax2.plot(error_bar_bincentre_down,np.linspace(1 - settings["limit_y"]*0.85, 1 - settings["limit_y"]*0.85, len(error_bar_bincentre_down)), 'kv', markersize='8')
    # tick setting
    ax1.tick_params(labelsize=16)
    ax2.tick_params(labelsize=16)
    if not settings['log_y']:
        ax1.minorticks_on()
    ax2.minorticks_on()
    # add label
    ax1.set_ylabel(settings['ylabelup'], fontsize=20)
    ax2.set_ylabel(settings['ylabeldown'], fontsize=20,labelpad=20)
    ax2.set_xlabel(settings['xlabel'], fontsize=20)
    plt.savefig(settings['filename'] + '.pdf', bbox_inches='tight', pad_inches = 0, )#transparent = True)
    if settings["pickle"]:
        with open(settings['filename'] + '.pkl', 'wb') as file:
            pickle.dump(fig, file)
    plt.show()
