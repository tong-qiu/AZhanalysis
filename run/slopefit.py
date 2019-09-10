import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from decimal import Decimal
from copy import deepcopy

def poly(x, *argv):
    s = 0
    for i, each in enumerate(argv):
        s += x**i * each
    return s

labelshift = 0
nbtag = 2
if nbtag == 1:
    #middle = 12
    middle = 6
    def fitfunction1(x, p0, p1, p2):# p5, p6):
        return poly(x, p0, p1, p2)# p5, p6)

    def fitfunction2(x, p0, p1):#, p2):# p5, p6):
        return poly(x, p0, p1)#, p2)# p5, p6)
if nbtag == 2:
    labelshift = 0.55
    #middle = 10
    middle = 5
    def fitfunction1(x, p0, p1, p2):# p5, p6):
        return poly(x, p0, p1, p2)# p5, p6)

    def fitfunction2(x, p0, p1):#, p2):# p5, p6):
        return poly(x, p0, p1)#, p2)# p5, p6)
data_point = []
mc_point = []
mc_error = []
bin_edge = []
mc_point_z = []
mc_error_z = []

filename = "pTV-mbbcut-" + str(nbtag) + "tag"
# load data
with open("output/t_make_plot_rescale/" + filename + ".csv") as f:
    for each_line in f:
        data_tem = each_line.split(",")
        
        if len(data_tem) > 3:
            bin_edge.append(float(data_tem[0]))
            data_point.append(float(data_tem[1]))
            mc_point.append(float(data_tem[2]))
            mc_error.append(float(data_tem[3]))
            mc_point_z.append(float(data_tem[4]))
            mc_error_z.append(float(data_tem[5]))
        else:
            bin_edge.append(float(data_tem[0]))

#convert to numpy
data_point = np.array(data_point)
mc_point = np.array(mc_point)
mc_error =  np.array(mc_error)
bin_edge = np.array(bin_edge)
mc_error_z =  np.array(mc_error_z)
mc_point_z = np.array(mc_point_z)
mc_error_other = np.sqrt(mc_error**2 - mc_error_z**2)

# calculate bin centre
bin_centre = (bin_edge[0:-1] + bin_edge[1:])/2

# pop zero
mask = data_point != 0
data_point = data_point[mask]
mc_point = mc_point[mask]
mc_error = mc_error[mask]
mc_point_z = mc_point_z[mask]
mc_error_z = mc_error_z[mask]
bin_centre = bin_centre[mask]
mc_error_other = mc_error_other[mask]

# calculate diff
mc_o = mc_point - mc_point_z
diff = (data_point - mc_o) /mc_point_z
diff_error = np.sqrt(data_point/ mc_point_z**2 + mc_error_other**2/ mc_point_z**2  + (data_point - mc_o)**2/  mc_point_z**4 * mc_error_z**2)

#fit
chi2nod = []
if middle > 0:
    popt1, pcov1 = curve_fit(fitfunction1, bin_centre[0:middle+1], diff[0:middle+1], sigma=diff_error[0:middle+1])
    r = diff[0:middle+1] - fitfunction1(bin_centre[0:middle+1], *popt1)
    chisq = sum((r / diff_error[0:middle+1]) ** 2)
    chi2nod.append(chisq/(-len(popt1) + len(bin_centre[0:middle+1])))
    print(chisq/(-len(popt1) + len(bin_centre[0:middle+1])))
diff_ = deepcopy(diff)
diff_error_ = deepcopy(diff_error)
if middle > 0:
    diff_[middle] = fitfunction1(bin_centre[middle], *popt1)
    diff_error_[middle] = 0.00000001
popt2, pcov2 = curve_fit(fitfunction2, bin_centre[middle:], diff_[middle:], sigma=diff_error_[middle:])
r = diff[middle:] - fitfunction2(bin_centre[middle:], *popt2)
chisq = sum((r / diff_error[middle:]) ** 2)
chi2nod.append(chisq/(-len(popt2) + len(bin_centre[middle:])))
print(chisq/(-len(popt2) + len(bin_centre[middle:])))
# print(bin_centre)
# print(diff)

#makeplot
plt.errorbar(bin_centre, diff, yerr=diff_error, fmt='o')
if middle > 0:
    xs = np.linspace(bin_centre[0], bin_centre[middle],100)
    xs = np.linspace(0, bin_centre[middle],100)
    plt.plot(xs, fitfunction1(xs, *popt1), 'r-')#, label='fit: p0=%5.3f, p1=%5.3f, p2=%5.3f, p3=%5.3f' % tuple(popt))
xs = np.linspace(bin_centre[middle], bin_centre[-1],100)
plt.plot(xs, fitfunction2(xs, *popt2), 'g-')
plt.xlabel("pTV(GeV)", fontsize=17)
plt.ylabel("reweight factor", fontsize=17)
#plt.ylim([diff.min(),10])
plt.yscale("log")
ax = plt.gca()
plt.text(0.05, 0.1 + labelshift, "red chi2/nod: " + "{:.5f}".format(chi2nod[0]), fontsize=15, transform=ax.transAxes)
plt.text(0.05, 0.03 + labelshift, "green chi2/nod: " + "{:.5f}".format(chi2nod[1]), fontsize=15, transform=ax.transAxes)
title1 = r"$\mathbf{ATLAS}$"
title1_1 = r"$\mathit{Internal}$"
title3 = "2 lep., " + str(nbtag) + " b-tag"
plt.text(0.05, 0.3 + labelshift, title1, fontsize=25, transform=ax.transAxes)
plt.text(0.327, 0.3 + labelshift, title1_1, fontsize=21, transform=ax.transAxes)
plt.text(0.05, 0.2 + labelshift, title3, fontsize=18, weight='bold', style='italic', transform=ax.transAxes)
plt.savefig("output/slopefit/" + filename + "polyfitresult.pdf" ,bbox_inches='tight', pad_inches = 0)
plt.show()

popt1 = popt1.tolist()
popt2 = popt2.tolist()

while len(popt1) < 5:
    popt1.append(0)
while len(popt2) < 5:
    popt2.append(0)
# print data
with open("output/slopefit/" + filename + "polyfitresult.csv", "w") as f:
    f.write(str(bin_centre[0]) + "," + str(bin_centre[middle]) + "," + str(bin_centre[-1]) + "\n")
    if middle > 0:
        for each in popt1:
            f.write(str(Decimal(repr(each))) + ',')
    else: 
        for each in popt2:
            f.write(str(0) + ',')
    f.write('\n')
    for each in popt2:
        f.write(str(Decimal(repr(each))) + ',')