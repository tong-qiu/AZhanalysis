import uproot
import numpy as np
import matplotlib.pyplot as plt

class truthsample():
    def __init__(self, path):
        self.path = path
        self.tree = uproot.open(path)["data"]
        self.weight = None
    def get_weight(self):
        if self.weight is None:
            output = []
            varname = b'mcEventWeights'
            content = self.tree.arrays([varname])[varname]
            for each in content:
                output.append(each)
            return np.array(output)
        return self.weight
    def get_variable(self, varname):
        output = []
        content = self.tree.arrays([varname])[varname]
        for each in content:
            output.append(each)
        return np.array(output)

def histplot_raw(datas, bins, labels, weights=None, removenorm=False, scale=1., **kwargs):
    settings = {
        "xlabel" : r"$m_{Vh}[GeV]$",
        "ylabel": 'Number of Events',
        "title1": r"ATLAS",# \newline Ptl next-leading, full cuts, 2 b-tags $",
        "title1_1": r"Internal",
        "title2": r"$\mathit{\sqrt{s}=13\:TeV,36.1\:fb^{-1}}$",# Ptl next-leading, full cuts, 2 b-tags $",
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
    ax.hist(datas, bins, histtype='step', fill=False, color=colors[0:len(datas)], weights=weights)
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
    fig.savefig(settings['filename'] + '.pdf', bbox_inches='tight', pad_inches = 0.25)
    
def main():
    masses = [220, 260, 340, 380]
    filename = "ggA"
    variables = [b"mA", b"mVH", b"mll", b"nbjets", b"ncjets", b"nljets", b"ptl1", b"etal1", b"ptb1", b"etab1"]
    objcollection = []
    curvecollection = {}
    curvecollection_weight = {}
    for each in variables:
        curvecollection[each] = []
        curvecollection_weight[each] = []
    for each in masses:
        thisobj = truthsample(filename + str(each) + ".root")
        objcollection.append(thisobj)
        for each_var in variables:
            entries = thisobj.get_variable(each_var)
            if np.mean(entries) > 100:
                entries = entries/1000.
            curvecollection[each_var].append(entries)
            curvecollection_weight[each_var].append(thisobj.get_weight())
    for each_var in variables:
        low = np.min(np.array(curvecollection[each_var]).flatten())
        high = np.max(np.array(curvecollection[each_var]).flatten())
        bins = np.linspace(low, high, num=25)
        if each_var in [b"nbjets", b"ncjets", b"nljets"]:
            bins = np.array(range(0, 10)) - 0.0001 
        histplot_raw(curvecollection[each_var], bins, masses, curvecollection_weight[each_var], title2=filename, title3="",
                     xlabel=each_var.decode("utf-8"), filename = filename + each_var.decode("utf-8"))



if __name__ == "__main__":
    main()