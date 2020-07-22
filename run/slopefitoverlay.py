import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal
from copy import deepcopy


def poly(x, argv):
    s = 0
    for i, each in enumerate(argv):
        s += x**i * each
    return s

def fitfunction_1tag(x, p0, p1, p2, p3):
    y = np.zeros(len(x))
    y += p0 * (x <= p1)
    y += (p2 * (x - p1) + p0 ) * (x <= p3) * (x > p1)
    y +=  (p2 * (p3 - p1) + p0 ) * (x > p3)
    return y

def fitfunction_2tag(x, p0, p1, p2, p3):
    y = np.zeros(len(x))
    y += (p1 * x + p0 ) * (x <= p2)
    y += p3 * (x > p2)
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

if __name__ == '__main__':
    tag = 2
    p1s = get_slope_correction("output/slopefit/" + "pTH-highmbbcut-"+str(tag)+"tagpolyfitresult.csv")
    p2s = get_slope_correction("output/slopefit/" + "pTH-lowmbbcut-"+str(tag)+"tagpolyfitresult.csv")

    x = []
    y1 = []
    y2 = []
    labelshift=0.6

    x = np.linspace(0,900,1600)
    y1 = fitfunction_2tag(x, p1s[0], p1s[1], p1s[2], p1s[3])
    y2 = fitfunction_2tag(x, p2s[0], p2s[1], p2s[2], p2s[3])
    plt.plot(x, y1, 'g-', label = 'high sideband')
    line2 = plt.plot(x, y2, 'b-', label = 'low sideband')
    ax = plt.gca()
    ax.legend(loc=5, fontsize='xx-large',frameon=False)
    plt.xlabel(r"$p_{T}^{BB}$ [GeV]", fontsize=17)
    plt.ylabel("reweight factor", fontsize=17)
    # if limity:
    #     plt.ylim([-0.5,2.5])
    # #plt.yscale("log")
    
    # plt.text(0.05, 0.1 + labelshift, "$\chi^2$/ndf: " + "{:.5f}".format(chi2nod[0]), fontsize=15, transform=ax.transAxes)
    # # plt.text(0.05, 0.03 + labelshift, "green chi2/ndf: " + "{:.5f}".format(chi2nod[1]), fontsize=15, transform=ax.transAxes)
    title1 = r"ATLAS"
    title1_1 = r"Internal"
    title3 = "2 lep., " + str(tag) + " b-tag"
    plt.text(0.45, 0.3 + labelshift, title1, fontsize=25, transform=ax.transAxes, style='italic', fontweight='bold')
    plt.text(0.727, 0.3 + labelshift, title1_1, fontsize=25, transform=ax.transAxes)
    plt.text(0.45, 0.2 + labelshift, title3, fontsize=18, weight='bold', style='italic', transform=ax.transAxes)
    plt.savefig("output/slopefit/overlay " + str(tag) + ".pdf" ,bbox_inches='tight', pad_inches = 0)
    plt.show()