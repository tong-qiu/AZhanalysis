import os
import uproot
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.stats
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)

def z2p(z):
    return 0.5 * (1 - math.erf(z / np.sqrt(2)))
ifhvt = True
alldirs = os.listdir("OutputDir")
masses = []
ps = []
for eachdir in alldirs:
    if not os.path.exists("OutputDir/"+ eachdir + "/outputLimit/asymptotics/test_CL95.root"):
        continue
    f = uproot.open("OutputDir/"+ eachdir + "/outputLimit/asymptotics/test_CL95.root")["stats"]
    mass = f.array("mass")[0]
    pobs = f.array("pb_obs")[0]
    masses.append(mass)
    ps.append(pobs)
print(masses)
print(ps)
zipped_pairs = zip(masses, ps) 
ps = [x for _, x in sorted(zipped_pairs)] 
masses = sorted(masses)
print(masses)
print(ps)
plt.plot(masses, ps, ".-", color="navy")
for z in [0, 1, 2, 3]:
    x = np.linspace(masses[0], masses[-1], 1000)
    y = [z2p(z)]*len(x)
    plt.plot(x, y, ':', color='silver')
    plt.annotate(str(z) + r' $\sigma$', xy=(x[-70],y[1]), xycoords='data', color='grey')
if ifhvt:
    plt.xlabel(r"$m_{Z'}$ [GeV]", fontsize=17)
else:
    plt.xlabel(r"$m_{A}$ [GeV]", fontsize=17)
plt.ylabel("p-value", fontsize=17)
plt.ylim([z2p(3.1), 10])
plt.xlim([masses[0], masses[-1]])
plt.yscale("log")
ax = plt.gca()
ax.xaxis.set_minor_locator(MultipleLocator(50))
title1 = r"ATLAS"
title1_1 = r"Internal"
title3 = r"$\sqrt{s}=13\:TeV,139\:fb^{-1}$"
plt.text(0.05, 0.9, title1, fontsize=18, transform=ax.transAxes, style='italic', fontweight='bold')
plt.text(0.25, 0.9, title1_1, fontsize=18, transform=ax.transAxes)
plt.text(0.05, 0.82, title3, fontsize=14, weight='bold', style='italic', transform=ax.transAxes)
if ifhvt:
    plt.text(0.05, 0.75, "HVT Z', 2 lep.", fontsize=14, transform=ax.transAxes)
    plt.savefig("hvt.pdf" ,bbox_inches='tight', pad_inches = 0.02)
else:
    plt.text(0.05, 0.75, "2HDM ggA, 2 lep.", fontsize=14, transform=ax.transAxes)
    plt.savefig("ggA.pdf" ,bbox_inches='tight', pad_inches = 0.02)
plt.show()