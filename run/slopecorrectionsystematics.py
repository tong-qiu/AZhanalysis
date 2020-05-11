import matplotlib.pyplot as plt
import numpy as np

def get_slope_correction(path):
    p1s = []
    with open(path) as f:
        for each_line in f:
            p1s = each_line.split(',')
            for i in range(len(p1s)-1):
                p1s[i] = float(p1s[i])
            break
    return p1s

def fitfunction(x, p0, p1, p2, p3):
    y = np.zeros(len(x))
    y += p0 * (x <= p1)
    y += (p2 * (x - p1) + p0 ) * (x <= p3) * (x > p1)
    y +=  (p2 * (p3 - p1) + p0 ) * (x > p3)
    return y

def sysfunction(x,nominal, high, low):
    up = []
    down = []
    ynominal = fitfunction(x, nominal[0], nominal[1], nominal[2], nominal[3])
    for each in x:
        pass



ntag = 2
labelshift = 0.5

pnominal = get_slope_correction("output/slopefit/" + "pTH-mbbcut-"+str(ntag)+"tagpolyfitresult.csv")
phigh = get_slope_correction("output/slopefit/" + "pTH-highmbbcut-"+str(ntag)+"tagpolyfitresult.csv")
plow = get_slope_correction("output/slopefit/" + "pTH-lowmbbcut-"+str(ntag)+"tagpolyfitresult.csv")

xnominal = np.linspace(0, 800, 800)
ynominal = fitfunction(xnominal, pnominal[0], pnominal[1], pnominal[2], pnominal[3])
yhigh = fitfunction(xnominal, phigh[0], phigh[1], phigh[2], phigh[3])
ylow = fitfunction(xnominal, plow[0], plow[1], plow[2], plow[3])


plt.plot(xnominal, ynominal, 'k-')
plt.plot(xnominal, yhigh, 'r-')
plt.plot(xnominal, ylow, 'b-')
plt.xlabel(r"$p_{TH}$ [GeV]", fontsize=17)
plt.ylabel("reweight factor", fontsize=17)
#plt.ylim([0.5,1.5])
#plt.yscale("log")
ax = plt.gca()
title1 = r"ATLAS"
title1_1 = r"Internal"
title3 = "2 lep., " + str(ntag) + " b-tag"
plt.text(0.35, 0.3 + labelshift, title1, fontsize=25, transform=ax.transAxes, style='italic', fontweight='bold')
plt.text(0.627, 0.3 + labelshift, title1_1, fontsize=25, transform=ax.transAxes)
plt.text(0.35, 0.2 + labelshift, title3, fontsize=18, weight='bold', style='italic', transform=ax.transAxes)
plt.savefig("output/slopefit/sys-" + str(ntag)  + "tag.pdf" ,bbox_inches='tight', pad_inches = 0)
plt.show()

