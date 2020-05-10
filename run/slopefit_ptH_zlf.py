import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from decimal import Decimal
from copy import deepcopy

def fitfunction(x, p0, p1, p2, p3):
    y = np.zeros(len(x))
    y += p0 * (x <= p1)
    y += (p2 * (x - p1) + p0 ) * (x <= p3) * (x > p1)
    y +=  (p2 * (p3 - p1) + p0 ) * (x > p3)
    return y
def fitfunction_real(x, p0, p1, p2, p3):
    if x <= p1:
        return p0
    if x <=p3:
        return p2 * (x - p1) + p0
    else:
        return p2 * (p3 - p1) + p0
def fitfunction2(x, p0, p1, p2, p3, p4):
    y = np.zeros(len(x))
    y += p0 * (x <= p1)
    y += (p2 * (x - p1) + p0 ) * (x <= p3) * (x > p1)
    y +=  (p2 * (p3 - p1) + p0 ) * (x > p3)
    return y
def poly(x, *argv):
    s = 0
    for i, each in enumerate(argv):
        s += x**i * each
    return s

labelshift = 0
nbtag = 1
highlow = "low"
if nbtag == 1:
    labelshift = 0.55
    g1 = 1
    g2 = 0
    g3 = 0
    g4 = 450
    # def fitfunction1(x, p0, p1, p2):# p5, p6):
    #     return poly(x, p0, p1, p2)# p5, p6)

    # def fitfunction2(x, p0, p1):#, p2):# p5, p6):
    #     return poly(x, p0, p1)#, p2)# p5, p6)
if nbtag == 2:
    labelshift = 0.55
    if highlow == "low":
        plt.ylim(top=3)
        plt.ylim(bottom=-0.4)
    g1 = 1.1
    g2 = 180
    g3 = 0
    g4 = 600
    # def fitfunction1(x, p0, p1, p2):# p5, p6):
    #     return poly(x, p0, p1, p2)# p5, p6)

    # def fitfunction2(x, p0, p1):#, p2):# p5, p6):
    #     return poly(x, p0, p1)#, p2)# p5, p6)
# if nbtag == 1:
#     labelshift = 0.55
#     #middle = 12
#     middle = 10
#     def fitfunction1(x, p0, p1, p2):# p5, p6):
#         #return poly(x, p0, p1, p2)# p5, p6)
#         return poly(x, p0, p1)

#     def fitfunction2(x, p0, p1):#, p2):# p5, p6):
#         #return poly(x, p0, p1)#, p2)# p5, p6)
#         return poly(x, p0)
# if nbtag == 2:
#     labelshift = 0.55
#     #middle = 10
#     middle = 6
#     def fitfunction1(x, p0, p1, p2):# p5, p6):
#         return poly(x, p0, p1)# p5, p6)

#     def fitfunction2(x, p0, p1):#, p2):# p5, p6):
#         return poly(x, p0, p1)#, p2)# p5, p6)
data_point = []
mc_point = []
mc_error = []
bin_edge = []
mc_point_z = []
mc_error_z = []

filename = "pTH-" + highlow + "mbbcut-" + str(nbtag) + "tag"
print(filename)
# load data
with open("output/t_make_plot_rescale_zlf/" + filename + ".csv") as f:
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

# xnew = np.arange(0, 1300, 1)
# tck = interpolate.splrep(bin_centre, diff, k=5, w =1/diff_error, s = 40)
# ynew = interpolate.splev(xnew, tck, der=0)

# cs = interpolate.CubicSpline(bin_centre, diff, bc_type =((1, 0.0), (1, 0.0)))

# plt.figure()
# # plt.plot(xnew, ynew, 'k')
# plt.plot(xnew, cs(xnew), 'k')
# plt.errorbar(bin_centre, diff, yerr=diff_error, fmt='o')
# #plt.yscale("log")
# plt.title('Cubic-spline interpolation')
# plt.show()

#fit
chi2nod = []
upper = len(diff_error)
for i, each in enumerate(diff_error):
    if i == 0:
        continue
    if each > 0.08:
        upper = i
        break
print(bin_centre[upper])
popt1, pcov1 = curve_fit(fitfunction, bin_centre, diff, sigma=diff_error, p0=[1,0,0,bin_centre[upper]-50], bounds=((-np.inf, -np.inf, -np.inf, 0), (np.inf, np.inf, np.inf, bin_centre[upper])) )
print(popt1)
print(pcov1)
print(np.sqrt(np.diag(pcov1)))
r = diff - fitfunction(bin_centre, *popt1)
chisq = sum((r / diff_error) ** 2)
chi2nod.append(chisq/(-len(popt1) + len(bin_centre)))
print(chisq/(-len(popt1) + len(bin_centre)))


# makeplot
plt.errorbar(bin_centre, diff, yerr=diff_error, fmt='k.')

xs = np.linspace(0, bin_centre[-1],10000)
ys1 = []
ys2 = []
ys3 = []
xs1 = []
xs2 = []
xs3 = []
for each in xs:
    if each <= popt1[1]:
        xs1.append(each)
        ys1.append(fitfunction_real(each, popt1[0], popt1[1], popt1[2], popt1[3]))
    elif each <= popt1[3]:
        xs2.append(each)
        ys2.append(fitfunction_real(each, popt1[0], popt1[1], popt1[2], popt1[3]))
    else:
        xs3.append(each)
        ys3.append(fitfunction_real(each, popt1[0], popt1[1], popt1[2], popt1[3]))
plt.plot(xs1, ys1, 'g-')
plt.plot(xs2, ys2, 'r-')
plt.plot(xs3, ys3, 'b-')
plt.xlabel(r"$p_{TH}$ [GeV]", fontsize=17)
plt.ylabel("reweight factor", fontsize=17)
#plt.ylim([0.5,1.5])
#plt.yscale("log")
ax = plt.gca()
plt.text(0.05, 0.1 + labelshift, "$\chi^2$/ndf: " + "{:.5f}".format(chi2nod[0]), fontsize=15, transform=ax.transAxes)
# plt.text(0.05, 0.03 + labelshift, "green chi2/ndf: " + "{:.5f}".format(chi2nod[1]), fontsize=15, transform=ax.transAxes)
title1 = r"ATLAS"
title1_1 = r"Internal"
title3 = "2 lep., " + str(nbtag) + " b-tag"
plt.text(0.05, 0.3 + labelshift, title1, fontsize=25, transform=ax.transAxes, style='italic', fontweight='bold')
plt.text(0.327, 0.3 + labelshift, title1_1, fontsize=25, transform=ax.transAxes)
plt.text(0.05, 0.2 + labelshift, title3, fontsize=18, weight='bold', style='italic', transform=ax.transAxes)
plt.savefig("output/slopefit_zlf/" + filename + "polyfitresult.pdf" ,bbox_inches='tight', pad_inches = 0)
plt.show()

popt1 = popt1.tolist()
pcov1 = np.sqrt(np.diag(pcov1)).tolist()

# print data
with open("output/slopefit_zlf/" + filename + "polyfitresult.csv", "w") as f:
    for each in popt1:
        f.write(str(Decimal(repr(each))) + ',')
    f.write('\n')
    for each in pcov1:
        f.write(str(Decimal(repr(each))) + ',')