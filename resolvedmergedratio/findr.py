import os
import uproot

def get_error(var1, var2, nom1, nom2):
    output = 0
    for each_var1, each_var2 in zip(var1, var2):
        output += (abs(each_var1/each_var2 - nom1/nom2)/(nom1/nom2))**2
    return output**0.5

def getweight(infiles, sample, hists1, hists2):
    out1 = 0
    out2 = 0
    for each_file in infiles:
        f = uproot.open(each_file)
        for each_key in f.keys():
            each_keystring = each_key.decode("utf-8")
            if sample not in each_keystring:
                continue
            for each_region in hists1:
                if each_region in each_keystring:
                    #print(each_keystring)
                    out1 += f[each_key].allvalues.sum()
                    break
            for each_region in hists2:
                if each_region in each_keystring:
                    #print(each_keystring)
                    out2 += f[each_key].allvalues.sum()
                    break
    return (out1, out2)

def main():
    # region1 = ["1tag1pfat0pjet_0ptv_SR_noaddbjetsr_mVH", "2tag1pfat0pjet_0ptv_SR_noaddbjetsr_mVH"]
    # region2 = ["1tag2pjet_0ptv_SR_mVH", "2tag2pjet_0ptv_SR_mVH"]
    # region1 = ["1tag2pjet_0ptv_SR_mVH", "2tag2pjet_0ptv_SR_mVH"]
    # region2 = ["3ptag2pjet_0ptv_SR_mVH"]
    altttbar = "ttbar_"
    altsa1 = ["a/hist-ttbar_dilep_aMCatNLOPwPy8_alternative-0.root", "a/hist-ttbar_semilep_aMCatNLOPwPy8_alternative-0.root"]
    altsa2 = ["a/hist-ttbar_PwHwg7_DiLep_alternative-0.root", "a/hist-ttbar_PwHwg7_SingleLep_alternative-0.root"]
    # altsa1 = ["e/hist-ttbar_dilep_aMCatNLOPwPy8-0.root", "e/hist-ttbar_semilep_aMCatNLOPwPy8-0.root"]
    # altsa2 = ["e/hist-ttbar_PwHwg7_DiLep-0.root", "e/hist-ttbar_PwHwg7_SingleLep-0.root"]
    region1 = ["1tag2pjet", "2tag2pjet"]
    region2 = ["3ptag2pjet"]
    #region2 = ["1tag2pjet_0ptv_topemucr_mVH", "2tag2pjet_0ptv_topemucr_mVH"]
    nom1, nom2 = getweight(["run2dbla.root"], "ttbar", region1, region2)
    nom1ttbarB, nom2ttbarB = getweight(["run2dbla.root"], altttbar, region1, region2)
    print(nom1ttbarB/nom1, nom2ttbarB/nom2)

    #region2 = ["3tag2pjet_0ptv_SR_mVH", "4ptag2pjet_0ptv_SR_mVH"]
    region2 = ["3tag2pjet", "4ptag2pjet"]
    nom1, nom2 = getweight(altsa1, "ttbar", region1, region2)
    nom1ttbarB, nom2ttbarB = getweight(altsa1, altttbar, region1, region2)
    print(nom1ttbarB/nom1, nom2ttbarB/nom2)

    nom1, nom2 = getweight(altsa2, "ttbar", region1, region2)
    nom1ttbarB, nom2ttbarB = getweight(altsa2, altttbar, region1, region2)
    print(nom1ttbarB/nom1, nom2ttbarB/nom2)

    nom1, nom2 = getweight(["d/hist-ttbar_nonallhad_PwPy8-0.root"], "ttbar", region1, region2)
    nom1ttbarB, nom2ttbarB = getweight(["d/hist-ttbar_nonallhad_PwPy8-0.root"], altttbar, region1, region2)
    print(nom1ttbarB/nom1, nom2ttbarB/nom2)

if __name__ == "__main__":
    main()
