import uproot
import pickle
from package.easyload import *
from package.cut import *
from scipy.interpolate import *
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use( 'tkagg' )
import numpy as np
#import ROOT


loadeasytree = False
makeplotsbeforecorrection = True
def pickleit(obj, path):
    outfile = open(path, 'wb')
    pickle.dump(obj, outfile)
    outfile.close()

def unpickleit(path):
    infile = open(path, 'rb')
    output = pickle.load(infile)
    infile.close()
    return output

def findzll(inputsample):
    safeidx = 0
    for i in range(len(inputsample)):
        if inputsample[i].alias == "Zlljet":
            if safeidx > 0:
                print("error")
                exit(1)
            safeidx = i
    return inputsample[safeidx]

def autobin(data_list, bins, alias=None, variable=b"pTV"):
    new_data = None
    for each in data_list:
        if (alias is None and alias != "data") or each.alias == alias:
            if new_data is None:
                new_data = each
            else:
                new_data = new_data + each
    height, error = new_data.binned_weight_variation(variable, bins, scale=1000.)
    newbin = [bins[0]]

    sum_weight = 0
    sum_error = 0
    for i in range(len(height)):
        if i == len(height) - 1:
            newbin.append(bins[i+1])
            break
        sum_weight += height[i]
        sum_error += error[i]

        if (sum_error**0.5)/sum_weight < 0.35 and sum_weight > 3:
            newbin.append(bins[i+1])
            sum_weight = 0
            sum_error = 0
    return newbin

allregions = ["Zbb_1tag2pjet_0ptv_SR_mVH", "Zbb_2tag2pjet_0ptv_SR_mVH",
              "Zbc_1tag2pjet_0ptv_SR_mVH", "Zbc_2tag2pjet_0ptv_SR_mVH",
              "Zbl_1tag2pjet_0ptv_SR_mVH", "Zbl_2tag2pjet_0ptv_SR_mVH",
              "Zcc_1tag2pjet_0ptv_SR_mVH", "Zcc_2tag2pjet_0ptv_SR_mVH",
              "Zcl_1tag2pjet_0ptv_SR_mVH", "Zcl_2tag2pjet_0ptv_SR_mVH",
              "Zl_1tag2pjet_0ptv_SR_mVH", "Zl_2tag2pjet_0ptv_SR_mVH"]


if __name__ == "__main__":
    if loadeasytree:
        mysamplembbcr1tag, t2 = get_sample(["mbbcr", "resolved", "1tag"])
        pickleit((mysamplembbcr1tag, t2), "pickle/mbbcr1tag")
        mysamplembbcr2tag, t2 = get_sample(["mbbcr", "resolved", "2tag"])
        pickleit((mysamplembbcr2tag, t2), "pickle/mbbcr2tag")
        mysamplesr1tag, t2 = get_sample(["sr", "resolved", "1tag"])
        pickleit((mysamplesr1tag, t2), "pickle/sr1tag")
        mysamplesr2tag, t2 = get_sample(["sr", "resolved", "2tag"])
        pickleit((mysamplesr2tag, t2), "pickle/sr2tag")
        exit(1)
    else:
        mysamplembbcr1tag, t2mbbcr1tag = unpickleit("pickle/mbbcr1tag")
        mysamplembbcr1tag = rescale(mysamplembbcr1tag)
        mysamplembbcr2tag, t2mbbcr2tag = unpickleit("pickle/mbbcr2tag")
        mysamplembbcr2tag = rescale(mysamplembbcr2tag)
        mysamplesr1tag, t2sr1tag = unpickleit("pickle/sr1tag")
        mysamplesr1tag = rescale(mysamplesr1tag)
        mysamplesr2tag, t2sr2tag = unpickleit("pickle/sr2tag")
        mysamplesr2tag = rescale(mysamplesr2tag)

    tag = 2
    theone = mysamplembbcr2tag
    for i in range(len(theone)):
        #theone[i].cut(cut_lowmbb)
        theone[i].pth()
    
    theone_main = copy.deepcopy(theone)
    theone_high = copy.deepcopy(theone)
    theone_low = copy.deepcopy(theone)
    for i in range(len(theone_high)):
        theone_high[i].cut(cut_highmbb)
    for i in range(len(theone_low)):
        theone_low[i].cut(cut_lowmbb)

    if makeplotsbeforecorrection:
        bins = range(0,1000,50)
        title3="mBBcr " + str(tag) +" btags"
        chi2, nod = stackplot(theone_main, b'pTV', bins, 1000.,
                xlabel=r"$p_{TV} GeV]$", print_height=True, filename="test1",
                title2=t2mbbcr1tag, title3 = title3, auto_colour=False, limit_y=0.5, upper_y=2.0, log_y=False, printzpjets=True, chi2=True)
        exit(1)
    # theone_main = slopecorrection(theone_main, csv="../run/output/slopefit/" + "pTV-mbbcut-"+str(tag)+"tagpolyfitresult.csv")
    # theone_high = slopecorrection(theone_high, csv="../run/output/slopefit/" + "pTV-highmbbcut-"+str(tag)+"tagpolyfitresult.csv")
    # theone_low = slopecorrection(theone_low, csv="../run/output/slopefit/" + "pTV-lowmbbcut-"+str(tag)+"tagpolyfitresult.csv")
    bins = range(0,1000,50)
    title3="mBBcr " + str(tag) +" btags"
    chi2, nod = stackplot(theone_main, b'pTV', bins, 1000.,
            xlabel=r"$p_{TV} GeV]$", print_height=True, filename="test1",
            title2=t2mbbcr1tag, title3 = title3, auto_colour=False, limit_y=0.5, upper_y=2.0, log_y=False, printzpjets=True, chi2=True)
    exit(1)
    theone_main = findzll(theone_main)
    theone_high = findzll(theone_high)
    theone_low = findzll(theone_low)
    bins = range(0,3000,20)
    bins = autobin([theone_main], bins, alias="Zlljet", variable=b"mVH")
    bincentre = []
    for i in range(0, len(bins)-1):
        bincentre.append((bins[i] + bins[i+1])/2.)

    (theone_main_weight, theone_main_var) = theone_main.binned_weight_variation(b"mVH", bins, scale=1000.)
    (theone_high_weight, theone_high_var) = theone_high.binned_weight_variation(b"mVH", bins, scale=1000.)
    (theone_low_weight, theone_low_var) = theone_low.binned_weight_variation(b"mVH", bins, scale=1000.)

    theone_high_weight = np.array(theone_high_weight) / np.array(theone_main_weight) - 1
    theone_low_weight = np.array(theone_low_weight) / np.array(theone_main_weight) - 1
    theone_main_weight = np.array(theone_main_weight) * 0

    plt.plot(bincentre, theone_main_weight, label='nominal')
    plt.plot(bincentre, theone_high_weight, label='high')
    plt.plot(bincentre, theone_low_weight, label='low')
    plt.show()


    # mysample = rescale(mysample)
    # mysample = slopecorrection(mysample)



    # f = uproot.open("run2dbla.root")["Zbb_2tag1pfat0pjet_0ptv_mBBcr_noaddbjetsr_pTV"]


    # newhistofile = ROOT.TFile("test","recreate")
    # oldhistofile = ROOT.TFile("run2dbla.root","read")
    # oldhisto1tag = oldhistofile.Get("Zbb_2tag1pfat0pjet_0ptv_mBBcr_noaddbjetsr_pTV")