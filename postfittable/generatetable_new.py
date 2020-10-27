import pickle
import ROOT
from ctypes import *
import decimal
decimal.getcontext().rounding = "ROUND_HALF_UP"
import copy

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
        #return "ttbar+LF, stop, ttH"
        return "ttbar, stop, ttH"
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
    label = "llbb"
    if lepton == 0:
        label = "vvbb"
    selectedtags = [each for each in histocollection.keys() if label in each]
    selectedtags = sorted(selectedtags, key=getrank)
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
        for each_tag in selectedtags:
            value = 0
            error = 0
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
        f.write(" & " + getnumber(value, error))
    f.write(nm * " & ")
    f.write(r" \\" +  "\n")
    f.write(r"         \hline\hline"+ "\n")
    

def writeend(f):
    f.write(r"    \end{tabular}"+ "\n")
    f.write(r"    \end{footnotesize}"+ "\n")
    f.write(r"    \end{center}"+ "\n")
    f.write(r"\end{table}"+ "\n")

def main():
    with open("prefit.pickle", 'rb') as f:
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
        histocollection[histotem.region + "---" + histotem.tag] = copy.deepcopy(histotem)
        del histotem

    #print(bkgname)
    # print(inputfile.keys())
    # print(tags)


    with open("outtable.tex", "w") as f:
        writetitle(f, histocollection)
        writemain(f, sorted(bkgname, key=rankbkg), histocollection, 0)
        writemain(f, sorted(bkgname, key=rankbkg), histocollection, 2)
        writeend(f)

if __name__ == "__main__":
    main()