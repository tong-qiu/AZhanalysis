import pickle
import ROOT
from ctypes import *
import decimal
decimal.getcontext().rounding = "ROUND_HALF_UP"
import copy
from textstackplot import maketextstackplot
import numpy as np

bbA = False
class Histo:
    def __init__(self, tag, region):
        self.tag = tag
        self.region = region
        self.data = 0
        self.totalbkg = 0
        self.data_error = 0
        self.totalbkg_error = 0
        self.bkgs = {}

# \begin{table}[t]
#     \begin{center}
#     \begin{footnotesize}
#     \begin{tabular}{l|ccc|cc}
#     \hline\hline
#       & \multicolumn{3}{c|}{ Resolved }  & \multicolumn{2}{c}{ Merged }  \\

nresolved = 0
nmerged = 0

def writetitle(f, histocollection):
    global nresolved
    global nmerged
    nresolved0 = 0
    nmerged0 = 0
    nresolved2 = 0
    nmerged2 = 0
    for each in histocollection.keys():
        if "llbb" in each:
            if "1pfat" in each:
                nmerged2 += 1
            else:
                nresolved2 += 1
        elif "vvbb" in each:
            if "1pfat" in each:
                nmerged0 += 1
            else:
                nresolved0 += 1
    nresolved = max(nresolved0, nresolved2)
    nmerged = max(nmerged0, nmerged2)

    f.write(r"\begin{table}[t]"+ "\n")
    f.write(r"    \begin{center}" + "\n")
    f.write(r"    \begin{footnotesize}"+ "\n")
    f.write(r"    \begin{tabular}{l|" + "c" * nresolved + "|" + "c" * nmerged + "}"+ "\n")
    f.write(r"    \hline\hline"+ "\n")
    f.write(r"     & \multicolumn{" + str(nresolved) + "}{c|}{ Resolved }  & \multicolumn{" + str(nmerged) + "}{c}{ Merged }  " + r"\\"+ "\n")

def getrank(name):
    out = 0
    if "1pfat" in name:
        out += 10000
    if "topemucr" in name:
        out += 1000
    if "mBBcr" in name:
        out += 2000
    if "_withadd" in name:
        out += 1000
    if "1tag" in name or "1ptag" in name:
        out += 1
    if "2tag" in name or "2ptag" in name:
        out += 2
    if "3tag" in name or "3ptag" in name:
        out += 3
    if "4tag" in name or "4ptag" in name:
        out += 4
    return out

def formatregion(name):
    out = ""
    if "1tag" in name:
        out += "1 $b$-tag"
    if "2tag" in name:
        out += "2 $b$-tag"
    if "3tag" in name:
        out += "3 $b$-tag"
    if "3ptag" in name:
        out += "3+ $b$-tag"
    if "4ptag" in name:
        out += "4+ $b$-tag"
    if "1ptag" in name:
        out += "1+2 $b$-tag"
    if "_withadd" in name:
        out += " add. $b$-tag"
    if "SR" in name:
        out += " SR"
    if "mBBcr" in name:
        out += " mBB CR"
    if "topemucr" in name:
        out += " top CR"
    return out

def getnumber(value, error):
    #np.format_float_positional(x, precision=2, unique=True, fractional=False, trim='-')
    if abs(value) < 0.0000001 and abs(error) < 0.0000001:
        return "0"
    return str(value) + r"$\pm$" + str(error)

def rankbkg(name):
    if name in ("top", "topLF"):
        return 0
    if name == "ttbarHF":
        return 1
    if name in ("Zhf"):
        return 2
    if name in ("Zclbl", "Zlf"):
        return 3
    if name in ("Zl"):
        return 4
    if name in ("Whf", "Wb"):
        return 5
    if name in ("Wclbl", "Wc"):
        return 6
    if name in ("Wl"):
        return 7
    if name in ("ttV"):
        return 8
    if name in ("Diboson", "VV"):
        return 9
    if name in ("VH125", "SMVH"):
        return 10
    return 1000

def formatbkg(name):
    if name == "ttbarHF":
        return "ttbar+HF"
    if name in ("top", "topLF"):
        if bbA:
            return "ttbar+LF, stop, ttH, ttV"
        return "ttbar, stop, ttH, ttV"
    if name in ("Wclbl", "Wc"):
        return "W+(bl, cl)"
    if name in ("Whf", "Wb"):
        return "W+(bb, bc, cc)"
    if name in ("Wl"):
        return "W+l"
    if name in ("Zclbl", "Zlf"):
        return "Z+(bl, cl)"
    if name in ("Zhf"):
        return "Z+(bb, bc, cc)"
    if name in ("Zl"):
        return "Z+l"
    if name in ("VH125", "SMVH"):
        return "SM VH"
    if name in ("Diboson", "VV"):
        return "Diboson"
    return name

def writemain(f, bkgname, histocollection, lepton):
    outregionname = []
    outdata = []
    outsample = {}
    outsampleerror = []
    for each in bkgname:
        outsample[formatbkg(each)] = []

    label = "llbb"
    if lepton == 0:
        label = "vvbb"
    selectedtags = [each for each in histocollection.keys() if label in each]
    selectedtags = sorted(selectedtags, key=getrank)
    if(len(selectedtags)==0):
        return (outregionname, outdata, outsample, outsampleerror)
    f.write("         " + str(lepton) + "-lepton " )


    nr = nresolved
    nm = nmerged
    for each in selectedtags:
        if "1pfat" in each:
            nm -= 1
        else:
            nr -= 1
    #print(selectedtags)
    #print(nm, nr, nresolved, nmerged)

    beginmerged = False
    for each_tag in selectedtags:
        if not beginmerged and "1pfat" in each_tag:
            beginmerged = True
            f.write(nr * " & ")
        f.write(" & " + formatregion(each_tag))
    f.write(nm * " & ")
    f.write(r" \\" +  "\n")
    f.write("         " + r"\hline" + "\n")

    # write number for each background
    for each_bkg in bkgname:
        # check empty background name
        value = 0
        error = 0
        for each_tag in selectedtags:
            if each_bkg in histocollection[each_tag].bkgs:
                value += histocollection[each_tag].bkgs[each_bkg][0]
                error += histocollection[each_tag].bkgs[each_bkg][1]

        if value == 0 and error == 0:
            continue

        # begin writing
        f.write("         " + formatbkg(each_bkg) + "   ")
        beginmerged = False
        for each_tag in selectedtags:
            if not beginmerged and "1pfat" in each_tag:
                beginmerged = True
                f.write(nr * " & ")
            value = 0
            error = 0
            if each_bkg in histocollection[each_tag].bkgs:
                value = histocollection[each_tag].bkgs[each_bkg][0]
                error = histocollection[each_tag].bkgs[each_bkg][1]
                # outsample[formatbkg(each_bkg)].append(float(value))
                if each_tag not in outregionname:
                    outregionname.append(each_tag)
            outsample[formatbkg(each_bkg)].append(float(value))
            f.write(" & " + getnumber(value, error))
        f.write(nm * " & ")
        f.write(r" \\" +  "\n")
    f.write("         " + r"\hline" + "\n")

    f.write("         " + "Total" + "   ")
    beginmerged = False
    for each_tag in selectedtags:
        if not beginmerged and "1pfat" in each_tag:
            beginmerged = True
            f.write(nr * " & ")
        value = histocollection[each_tag].totalbkg
        error = histocollection[each_tag].totalbkg_error
        outsampleerror.append(float(error))
        f.write(" & " + getnumber(value, error))
    f.write(nm * " & ")
    f.write(r" \\" +  "\n")

    f.write("         " + "Data" + "   ")
    beginmerged = False
    for each_tag in selectedtags:
        if not beginmerged and "1pfat" in each_tag:
            beginmerged = True
            f.write(nr * " & ")
        value = histocollection[each_tag].data
        error = histocollection[each_tag].data_error
        outdata.append(float(value))
        f.write(" & " + getnumber(value, error))
    f.write(nm * " & ")
    f.write(r" \\" +  "\n")
    f.write(r"         \hline\hline"+ "\n")

    for i in range(len(outregionname)):
        outstring = ""
        if "llbb" in outregionname[i]:
            outstring += "2L "
        if "vvbb" in outregionname[i]:
            outstring += "0L "
        if "1pfat" in outregionname[i]:
            outstring += "Merged "
        else:
            outstring += "Resolved "
        outstring += formatregion(outregionname[i])
        outstring.replace("$", "")
        outregionname[i] = outstring
    
    # print(outregionname)
    # print(outdata)
    # print(outsample)
    # print(outsampleerror)
    return (outregionname, outdata, outsample, outsampleerror)
    

def writeend(f):
    f.write(r"    \end{tabular}"+ "\n")
    f.write(r"    \end{footnotesize}"+ "\n")
    f.write(r"    \end{center}"+ "\n")
    f.write(r"\end{table}"+ "\n")

def loadandmaketable(filename, outname):
    with open(filename, 'rb') as f:
        inputfile = pickle.load(f)
    bkgname = set()
    tags = set()
    histocollection = {}
    # for each_plot in inputfile.keys():
    #     if "cr" in each_plot:
    #         print(inputfile[each_plot]['hsig'])
    #         cserror = c_double(0.)
    #         print (inputfile[each_plot]['hsig'].IntegralAndError(0, 99999, cserror, "width"))
    #         print (cserror.value)
    #exit(1)
    for each_plot in inputfile.keys():
        tag = each_plot.split("_")[0]
        if "noaddbjetsr" in each_plot:
            tag += "_noadd"
        if "topaddbjetcr" in each_plot:
            tag += "_withadd"
        tags.add(tag.replace("Ch", ""))
        region = each_plot.split("_")[2]
        print(region)
        if each_plot.split("_")[3] != "mVH":
            region += "_" + each_plot.split("_")[3]
        histotem = Histo(tag, region)
        cerror = c_double(0.)
        decimal.Decimal("0.645").quantize(decimal.Decimal("0.00"))
        histotem.data = decimal.Decimal(inputfile[each_plot]['hdata'].IntegralAndError(0, 99999, cerror, "width")).quantize(decimal.Decimal("0"))
        #print(inputfile[each_plot]['hdata'].IntegralAndError(0, (inputfile[each_plot]['hdata'].GetNbinsX()), cerror, "width"), cerror.value)
        #print(decimal.Decimal(inputfile[each_plot]['hdata'].IntegralAndError(0, 99999, cerror, "width")).quantize(decimal.Decimal("0")), decimal.Decimal(cerror.value).quantize(decimal.Decimal("0")))
        histotem.data_error = decimal.Decimal(cerror.value).quantize(decimal.Decimal("0"))
        histotem.totalbkg = decimal.Decimal(inputfile[each_plot]['hbkg'].IntegralAndError(0, 99999, cerror, "width")).quantize(decimal.Decimal("0.0"))
        histotem.totalbkg_error = decimal.Decimal(cerror.value).quantize(decimal.Decimal("0.0"))
        for each_bkg in inputfile[each_plot]['hbkgs']:
            bkg_tem = each_bkg.GetName().split("_")[3]
            total_tem = each_bkg.IntegralAndError(0, 99999, cerror, "width")
            # if bkg_tem in histotem.bkgs:
            #     exit(1)
            histotem.bkgs[bkg_tem] = [decimal.Decimal(total_tem).quantize(decimal.Decimal("0.0")), decimal.Decimal(cerror.value).quantize(decimal.Decimal("0.0"))]
            bkgname.add(bkg_tem)

        # print(histotem.bkgs)
        # if bbA:
        #     if "topLF" in histotem.bkgs and "top" in  histotem.bkgs:
        #         print("here")
        #         totaltop = decimal.Decimal(float(histotem.bkgs["topLF"][0]) + float(histotem.bkgs["top"][0])).quantize(decimal.Decimal("0.0"))
        #         totalerror = decimal.Decimal(float(histotem.bkgs["topLF"][1])**2. + float(histotem.bkgs["top"][1])**2.)**0.5.quantize(decimal.Decimal("0.0"))
        #         histotem.bkgs["top"] = [totaltop, totalerror]
        #         histotem.bkgs.pop("topLF")
        #         bkgname.remove("topLF")
        histocollection[histotem.region + "---" + histotem.tag] = copy.deepcopy(histotem)
        del histotem

    #print(bkgname)
    # print(inputfile.keys())
    # print(tags)

    with open(outname, "w") as f:
        # print(sorted(bkgname, key=rankbkg))
        writetitle(f, histocollection)
        L0outregionname, L0outdata, L0outsample, L0outsampleerror = writemain(f, sorted(bkgname, key=rankbkg), histocollection, 0)
        L2outregionname, L2outdata, L2outsample, L2outsampleerror = writemain(f, sorted(bkgname, key=rankbkg), histocollection, 2)
        writeend(f)

        regionname = L0outregionname + L2outregionname
        outdata = L0outdata + L2outdata
        sampleerror = L0outsampleerror + L2outsampleerror

        for each_key in L0outsample.keys():
            if len(L0outsample[each_key]) == 0:
                L0outsample[each_key] = [0.] * len(L0outdata)
        for each_key in L2outsample.keys():
            if len(L2outsample[each_key]) == 0:
                L2outsample[each_key] = [0.] * len(L2outdata)

        samplename = []
        samplevalue = []
        # print(L2outsample.keys())
        # print(L0outsample.keys())
        for each_key in L2outsample.keys():
            samplename.append(each_key)
            if len(L0outdata) == 0:
                samplevalue.append(L2outsample[each_key])
            else:
                samplevalue.append(L0outsample[each_key] + L2outsample[each_key])
                # print(L0outsample[each_key], L2outsample[each_key], each_key)
        return (samplevalue, sampleerror, outdata, regionname, samplename)

def main():
    global bbA
    bbA = False
    samplevalue, sampleerror, outdata, regionname, samplename = loadandmaketable("postfit.pickle", "postouttable.tex")
    presamplevalue, presampleerror, preoutdata, preregionname, presamplename = loadandmaketable("prefit.pickle", "preouttable.tex")
    # print(presamplevalue)
    prefitheights = np.array([0.] * len(outdata))
    for each in presamplevalue:
        # print(prefitheights)
        # print(each)
        prefitheights += np.array(each)
    # print(samplename)
    maketextstackplot(samplevalue, sampleerror, outdata, regionname, samplename, prefitheights, filename="test", log_y=True)

if __name__ == "__main__":
    main()