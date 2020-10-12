import os
import uproot

def get_error(var1, var2, nom1, nom2):
    output = 0
    for each_var1, each_var2 in zip(var1, var2):
        output += (abs(each_var1/each_var2 - nom1/nom2)/(nom1/nom2))**2
    return output**0.5

def getweight(infile, sample, hists1, hists2):
    f = uproot.open(infile)
    out1 = 0
    out2 = 0
    for each_key in f.keys():
        each_keystring = each_key.decode("utf-8")
        if sample not in each_keystring:
            continue
        for each_region in hists1:
            if each_region in each_keystring:
                print(each_keystring)
                out1 += f[each_key].allvalues.sum()
                break
        for each_region in hists2:
            if each_region in each_keystring:
                print(each_keystring)
                out2 += f[each_key].allvalues.sum()
                break
    return (out1, out2)

def main():
    region1 = ["1tag1pfat0pjet_0ptv_SR_noaddbjetsr_mVH", "2tag1pfat0pjet_0ptv_SR_noaddbjetsr_mVH"]
    region2 = ["1tag2pjet_0ptv_SR_mVH", "2tag2pjet_0ptv_SR_mVH"]
    region1 = ["1tag2pjet_0ptv_SR_mVH", "2tag2pjet_0ptv_SR_mVH"]
    region2 = ["1tag2pjet_0ptv_topemucr_mVH", "2tag2pjet_0ptv_topemucr_mVH"]
    nom1, nom2 = getweight("run2dbl.root", "ttbar", region1, region2)
    altdir = ["a"]
    altsamples = ["PwHwg7", "aMCat"]
    alt = [[0, 0] for each in range(len(altsamples))]
    for eachdir in altdir:
        allfiles = os.listdir(eachdir)
        for each_infile in allfiles:
            alttem1, alttem2 = getweight(eachdir + "/" + each_infile, "ttbar", region1, region2)
            for i, each_altsample in enumerate(altsamples):
                if each_altsample in each_infile:
                    alt[i][0] += alttem1
                    alt[i][1] += alttem2
    print(nom1, nom2)
    print(alt)
    var1 = [i[0] for i in alt]
    var2 = [i[1] for i in alt]
    print(var1, var2)
    print(get_error(var1, var2, nom1, nom2))
    # var1 = (9000, 9200)
    # var2 = (1000, 800)
    # nom1 = 9100
    # nom2 = 2000
    # print(get_error(var2, var1, nom2, nom1))
    # print(get_error(var1, var2, nom1, nom2))

if __name__ == "__main__":
    main()
