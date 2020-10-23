import matplotlib
import math
import matplotlib.pyplot as plt
import numpy as np
import pickle
from events import *
import matplotlib.font_manager as font_manager
import matplotlib.patches as mpatches
import copy
def nevent(asample):
    thesum = sum(asample.weight)
    return thesum#len(asample.data[varible_to_plot])

def stackplot(data_list, varible_to_plot, bins, scales=1., **kwargs):
    data_list = copy.deepcopy(data_list)
    # default label
    settings = {
        "xlabel" : r"$m_{Vh}[GeV]$",
        "ylabelup": 'Number of Events',
        "ylabeldown": "Data/Simulation",
        "title1": r"ATLAS",# \newline Ptl next-leading, full cuts, 2 b-tags $",
        "title1_1": r"Internal",
        #"title1_1": r"$\mathit{Working\:in\:progress}$",
        "title2": r"$\sqrt{s}=13\:TeV,36.1\:fb^{-1}$",# Ptl next-leading, full cuts, 2 b-tags $",
        #"title3": r"$\mathbf{2\;lep.,2\;b-tag}$",
        "title3": "2 lep., 2 b-tag",
        "title4": "",
        "filename": "Notsave",
        "auto_colour": True,
        "print_height": False,
        "limit_y": 0,
        "pickle": False,
        "upper_y": 1.7,
        "ncol":1,
        "log_y":False,
        "sys":False,
        "blind":False,
        "printzpjets":False,
        "printzpljets":False,
        "chi2":False,
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
    Zpjet = []
    if settings["printzpjets"] or settings["printzpljets"]:
        thealias = "Zlljet"
        if settings["printzpljets"]:
            thealias = "Z+lf"
        for each in data_list:
            if thealias == "Zlljet":
                if "Z+lf" in each.alias or "Z+hf" in each.alias:
                    Zpjet.append(each)
            else:
                if thealias in each.alias:
                    Zpjet.append(each)
        sigma2_z = []
        weight_in_bin_z = []
        sys2_z = []
        for each in Zpjet:
            if nevent(each) < 0:
                continue
            sigma2_z.append(each.variation(varible_to_plot, bins, scales))
            weight_in_bin_z.append(each.binned_weight(varible_to_plot, bins, scales))
            if settings["sys"] and len(each.systematics(varible_to_plot, bins, scales))>0:
                sys2_z.append(each.systematics(varible_to_plot, bins, scales))
        # calculate uncertainty
        sigma2_z = np.transpose(sigma2_z)
        if settings["sys"]:
            sys2_z = np.transpose(sys2_z)
            error_mc_z = np.sqrt(np.sum(sigma2_z, (1))+np.sum(sys2_z, (1)))
        else:
            error_mc_z = np.sqrt(np.sum(sigma2_z, (1)))


    data_list.sort(key=nevent)
    for each in data_list:
        print(nevent(each))
        if nevent(each) < 0:
            print(each.alias + " is empty.")
            continue
        if "data" not in each.alias:#np.mean(each.weight) != 1:# is not None:
            weight.append(each.weight)
            varibles_content.append(each.data[varible_to_plot]/scales)

            if each.alias == "Z+hf":
                alias.append("Z+(bb,bc,cc)")
            elif each.alias == "Z+lf":
                alias.append("Z+(bl,cl,l)")
            elif each.alias == "Wlvjet":
                alias.append("W+jets")
            elif each.alias == "smHiggs":
                alias.append("SM Higgs")
            elif each.alias == "ttbar":
                alias.append(r"$t\bar{t}$")
            elif each.alias == "singletop":
                alias.append("single top")
            else:
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
    ax1.text(0.05, 1.55 / 1.7, settings['title1'], fontsize=25, transform=ax1.transAxes, style='italic', fontweight='bold')
    ax1.text(0.227, 1.55/ 1.7, settings['title1_1'], fontsize=25, transform=ax1.transAxes)
    ax1.text(0.05, 1.40 / 1.7, settings['title2'], fontsize=23, transform=ax1.transAxes, style='italic', fontweight='bold')
    ax1.text(0.05, 1.26 / 1.7, settings['title3'], fontsize=18, weight='bold', style='italic', transform=ax1.transAxes)
    #ax1.text(0.05, 1.12 / 1.7, settings["title4"], fontsize=18, weight='bold', style='italic', transform=ax1.transAxes)

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

    if settings["print_height"] and settings['filename'] != "Notsave":
        if not (settings["printzpjets"] or settings["printzpljets"]):
            with open(settings["filename"]+".csv", "w") as f:
                #f.write("data" + "," + "MC" + ",\n")
                for each_data, each_mc, each_edges, each_errormc in zip(bin_heights, y_mc, bin_edges[0:-1], error_mc):
                    f.write(str(each_edges) + "," + str(each_data) + "," + str(each_mc) + "," + str(each_errormc) + ",\n")
                f.write(str(bin_edges[-1]) + ",")
        else:
            weight_z = None
            for each in weight_in_bin_z:
                if weight_z is None:
                    weight_z = np.array(each)
                else:
                    weight_z += np.array(each)
            with open(settings["filename"]+".csv", "w") as f:
                #f.write("data" + "," + "MC" + ",\n")
                for each_data, each_mc, each_edges, each_errormc, eachmc_z, eacherrormc_z in zip(bin_heights, y_mc, bin_edges[0:-1], error_mc, weight_z, error_mc_z):
                    f.write(str(each_edges) + "," + str(each_data) + "," + str(each_mc) + "," + str(each_errormc) + "," + str(eachmc_z) + "," + str(eacherrormc_z) + ",\n")
                f.write(str(bin_edges[-1]) + ",")


    # add legend
    sys_patch = mpatches.Patch(color='black', hatch='/////', fill=False, linewidth=0)
    ax1.legend([error_patch, sys_patch] + y_mcs[2], ["data", "uncertainty"] + alias, prop=font, frameon=False, ncol=settings["ncol"])
    if settings["log_y"]:
        ax1.set_yscale('log')
        #math.log10(max(bin_heights))
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
    for each_y, each_y_mc, each_menstd, each_error_bar_bincentre, each_error_mc in zip(bin_heights, y_mc, mean_std, bincenters, error_mc):
        error_bar_center.append(each_y / each_y_mc)
        # error_bar_size.append(np.sqrt(each_menstd**2 / each_y_mc**2 + each_menstd**4 / each_y_mc**4 * each_error_mc**2)) #np.sqrt(data_point/ mc_point**2  + data_point**2/  mc_point**4 * mc_error**2)
        error_bar_size.append(each_menstd / each_y_mc)
        error_bar_bincentre.append(each_error_bar_bincentre)

    ax2.errorbar(error_bar_bincentre, error_bar_center, color='k', linestyle="",
                 yerr=error_bar_size, fmt='o', markersize='8')
    
    # calculate chi2
    error_bar_center = np.array(error_bar_center)
    error_bar_size = np.array(error_bar_size)

    chi2 = 0
    nod = 0
    for each_y, each_sigma in zip(error_bar_center, error_bar_size):
        if math.isnan(each_y) or math.isnan(each_sigma) or each_sigma == 0:
            continue
        chi2 += ((each_y - 1)/each_sigma)**2
        nod += 1
    if settings["chi2"]:
        ax1.text(0.05, 1.12 / 1.7, "$\chi^2$/ndf: " + str(round(chi2/nod,3)), fontsize=18, weight='bold', style='italic', transform=ax1.transAxes)
    #chi2 = sum((error_bar_center/error_bar_size)**2)p
    #nod = len(error_bar_center)

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
    if settings['filename'] != "Notsave":
        plt.savefig(settings['filename'] + '.pdf', bbox_inches='tight', pad_inches = 0.1, )#transparent = True)
    if settings["pickle"] and settings['filename'] != "Notsave":
        with open(settings['filename'] + '.pkl', 'wb') as file:
            pickle.dump(fig, file)
    plt.show()
    return (chi2, nod)
