import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.patches as mpatches
import matplotlib.font_manager as font_manager

class sortingobj():
    def __init__(self, list, name):
        self.list = list
        self.name = name
    def eventsum(self):
        return sum(self.list)

def maketextstackplot(sampleheights, sampleerrors, dataheights, binlabels, samplenames, prefitheights, colours=[], **kwargs):
    # default label
    settings = {
        "ylabelup": 'Number of Events',
        "ylabeldown": "Data/Simulation",
        "title1": r"ATLAS",# \newline Ptl next-leading, full cuts, 2 b-tags $",
        "title1_1": r"Internal",
        #"title1_1": r"$\mathit{Working\:in\:progress}$",
        "title2": r"$\sqrt{s}=13\:TeV,139\:fb^{-1}$",# Ptl next-leading, full cuts, 2 b-tags $",
        #"title3": r"$\mathbf{2\;lep.,2\;b-tag}$",
        "title3": "",
        "title4": "",
        "filename": "Notsave",
        "auto_colour": True,
        "print_height": False,
        "limit_y": 0,
        "pickle": False,
        "upper_y": 1.7,
        "log_y":False,
        "blind":False,
        }
    for each_key in kwargs.items():
        settings[each_key[0]] = kwargs[each_key[0]]
    
    objlist = []
    for each_h, each_n in zip(sampleheights, samplenames):
        objlist.append(sortingobj(each_h, each_n))
    objlist = sorted(objlist, key= lambda x: x.eventsum())
    sampleheights = []
    samplenames = []
    for each in objlist:
        sampleheights.append(each.list)
        samplenames.append(each.name)


    samplelength = len(sampleheights[0])
    if samplelength != len(sampleerrors) or samplelength != len(dataheights) or samplelength != len(binlabels):
        print("Error: inputs must have the same lenght")
        exit(1)
    if len(samplenames) != len(sampleheights):
        print("Error: samples and sample names must have the same lenght")
        exit(1)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, gridspec_kw={'height_ratios':[4, 1, 1]}, figsize=(10, 12))
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.02)
    ax1.get_xaxis().set_ticks([])
    ax2.get_xaxis().set_ticks([])
    bins = np.array(range(0, samplelength + 1), dtype=np.float64)
    bincenter = (bins[0:-1] + bins[1:])/2

    # main plot
    if settings['auto_colour']:
        y_mcs = ax1.hist([bincenter.tolist()] * len(sampleheights), bins, weights=sampleheights,
                         stacked=True, label=samplenames)
    else:
        y_mcs = ax1.hist([bincenter.tolist()] * len(sampleheights), bins, weights=sampleheights,
                         stacked=True, label=samplenames, color=colours)
    error_patch = ax1.errorbar(bincenter, dataheights, color='k',
                               linestyle="", yerr=np.array(dataheights)**0.5, fmt='o', markersize='8')
    
    # make sure all plots have the same x axis range
    ax1.set_xlim([0, len(bincenter)])
    ax2.set_xlim([0, len(bincenter)])
    ax3.set_xlim([0, len(bincenter)])
    ax3.set_xticks(bincenter)
    # plot label
    ax3.set_xticklabels(binlabels, fontsize=14)
    plt.setp(ax3.get_xticklabels(), rotation=45, ha="right",
        rotation_mode="anchor")

    totalsampleheights = np.array([0.] * len(sampleheights[0]))
    for each in sampleheights:
        totalsampleheights += np.array(each)
    y_mc_max = max(totalsampleheights.tolist())


    # set range and plot label
    ax1.set_ylim([0,max([y_mc_max, max(dataheights)])* settings["upper_y"]])
    ax1.text(0.05, 1.55 / 1.7, settings['title1'], fontsize=20, transform=ax1.transAxes, style='italic', fontweight='bold')
    ax1.text(0.2, 1.55/ 1.7, settings['title1_1'], fontsize=20, transform=ax1.transAxes)
    ax1.text(0.05, 1.40 / 1.7, settings['title2'], fontsize=15, transform=ax1.transAxes, style='italic', fontweight='bold')
    ax1.text(0.05, 1.26 / 1.7, settings['title3'], fontsize=18, weight='bold', style='italic', transform=ax1.transAxes)
    if settings["log_y"]:
        ax1.set_yscale('log')
        ax1.set_ylim([0.1, 10**(math.log10(max([y_mc_max, max(dataheights)])) * settings["upper_y"])])
        ax1.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10,numticks=100))
        ax1.minorticks_on()

    # add legend
    font = font_manager.FontProperties(weight='bold', style='normal', size=12)
    sys_patch = mpatches.Patch(color='black', hatch='/////', fill=False, linewidth=0)
    legend_styles = [sys_patch] + y_mcs[2]
    legend_name = ["uncertainty"] + samplenames
    if not settings["blind"]:
        legend_styles = [error_patch, sys_patch] + y_mcs[2]
        legend_name = ["data", "uncertainty"] + samplenames
    ax1.legend(legend_styles, legend_name, prop=font, frameon=False, ncol=2)

    ax2height = []
    ax2error = []
    for each_mcheight, each_dataheight in zip(totalsampleheights, dataheights):
        ax2height.append(each_dataheight/each_mcheight)
        ax2error.append(each_dataheight**0.5/each_mcheight)
    ax2.plot([bins[0], bins[np.size(bins)-1]], [1, 1], linestyle='--', color='k')
    ax2.errorbar(bincenter, ax2height, color='k', linestyle="",
                 yerr=ax2error, fmt='o', markersize='8')

    # plot arrow
    if settings["limit_y"] > 0 and not settings["blind"]:
        ax2.set_ylim([1 - settings["limit_y"], 1 + settings["limit_y"]])
        error_bar_center = np.array(ax2height)
        error_bar_bincentre = np.array(bincenter)
        error_bar_bincentre_up = error_bar_bincentre[error_bar_center > 1 + settings["limit_y"]]
        error_bar_bincentre_down = error_bar_bincentre[1 - settings["limit_y"] > error_bar_center]
        ax2.plot(error_bar_bincentre_up,np.linspace(1 + settings["limit_y"]*0.85, 1 + settings["limit_y"]*0.85, len(error_bar_bincentre_up)), 'k^', markersize='8')
        ax2.plot(error_bar_bincentre_down,np.linspace(1 - settings["limit_y"]*0.85, 1 - settings["limit_y"]*0.85, len(error_bar_bincentre_down)), 'kv', markersize='8')
    elif settings["limit_y"] > 0:
        ax2.set_ylim([1 - settings["limit_y"], 1 + settings["limit_y"]])
    else:
        ax2.set_ylim([0.85, 1.15])

    # plot the error bars of MC
    i = -1
    for each_x, each_error, each_y_mc in zip(bincenter, sampleerrors, totalsampleheights):
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

    # plot the error bars of MC
    i = -1
    for each_x, each_error, each_y_mc in zip(bincenter, sampleerrors, totalsampleheights):
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

    # thirdplot
    ax3.plot([bins[0], bins[np.size(bins)-1]], [1, 1], linestyle='--', color='k')
    ax3height = []
    for each_postfit, each_prefit in zip(totalsampleheights, prefitheights):
        ax3height.append(each_prefit/each_postfit)
    ax3.hist(bincenter, bins=bins, weights=ax3height, histtype="step", color='k')
    ax3limit = max(np.abs(np.array(ax3height) - 1)) * 1.4
    ax3.set_ylim([1 - ax3limit, 1 + ax3limit])

    ax1.set_ylabel(settings['ylabelup'], fontsize=20)
    ax2.set_ylabel("data/postfit", fontsize=13)
    ax3.set_ylabel("prefit/postfit", fontsize=13)
    if settings['filename'] != "Notsave":
        plt.savefig(settings['filename'] + '.pdf', bbox_inches='tight', pad_inches = 0.1, )
    


def main():
    data = [200, 250, 230, 200, 250, 230]
    ttbar = [100, 225, 215, 100, 225, 215]
    zjet = [190, 20, 20, 190, 20, 20]
    error = [10, 7, 9, 10, 7, 9]
    prefitheights = [300, 250, 230, 310, 250, 230]
    binlabel = ["region1", "retion2", "region3", "region1", "retion2", "region3"]
    samples = [ttbar, zjet]
    samplenames = ["ttbar", "Z+jets"]
    maketextstackplot(samples, error, data, binlabel, samplenames, prefitheights, filename="test", upper_y=6)

if __name__ == "__main__":
    main()
