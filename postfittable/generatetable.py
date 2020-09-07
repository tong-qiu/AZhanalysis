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

def main():
    with open("prefit", 'rb') as f:
        inputfile = pickle.load(f)
    bkgname = set()
    tags = set()
    histocollection = []
    for each_plot in inputfile.keys():
        if "cr" in each_plot:
            print(inputfile[each_plot]['hsig'])
            cserror = c_double(0.)
            print (inputfile[each_plot]['hsig'].IntegralAndError(0, 99999, cserror, "width"))
            print (cserror.value)
    exit(1)
    for each_plot in inputfile.keys():
        tag = each_plot.split("_")[0]
        tags.add(tag)
        region = each_plot.split("_")[2]
        if each_plot.split("_")[3] != "mVH":
            region += "_" + each_plot.split("_")[3]
        histotem = Histo(tag, region)
        cerror = c_double(0.)
        decimal.Decimal("0.645").quantize(decimal.Decimal("0.00"))
        histotem.data = decimal.Decimal(inputfile[each_plot]['hdata'].IntegralAndError(0, 99999, cerror, "width")).quantize(decimal.Decimal("0"))
        print(inputfile[each_plot]['hdata'].IntegralAndError(0, (inputfile[each_plot]['hdata'].GetNbinsX()), cerror, "width"), cerror.value)
        print(decimal.Decimal(inputfile[each_plot]['hdata'].IntegralAndError(0, 99999, cerror, "width")).quantize(decimal.Decimal("0")), decimal.Decimal(cerror.value).quantize(decimal.Decimal("0")))
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
        histocollection.append(copy.deepcopy(histotem))
        del histotem

    print(bkgname)
    print(inputfile.keys())

    namemapper = [["top", "top"], ["Diboson", "Diboson"], ['Zl','Zl'], ['Zlf','Zhl'], ['Zhf','Zhf'], ['Wl','Wl'], ['Wc','Whl'], ['Wb','Whf'], ['SMVH', 'SM VH']]
    tagsegion = ['llbb1tag2pjetSR', 'llbb2tag2pjetSR', 'llbb3tag2pjetSR', 'llbb4ptag2pjetSR', 'llbb1ptag2pjettopemucr', 'llbb1tag1pfat0pjetSR_noaddbjetsr', 'llbb2tag1pfat0pjetSR_noaddbjetsr', "llbb1ptag1pfat0pjetSR_topaddbjetcr"]
    with open("outtable.tex", "w") as f:
        f.write(r"\begin{table}[t]"+ "\n")
        f.write(r"    \begin{center}" + "\n")
        f.write(r"    \begin{footnotesize}"+ "\n")
        f.write(r"    \begin{tabular}{l|ccccc|cccc}"+ "\n")
        f.write(r"    \hline\hline"+ "\n")
        f.write(r"      & \multicolumn{5}{c|}{ Resolved }  & \multicolumn{4}{c}{ Merged }  \\"+ "\n")
        f.write(r"    2-lepton &1 $b$-tag &2 $b$-tag &3 $b$-tag & 4+ $b$-tag & 1+2 $b$-tag topemucr & 1 $b$-tag &2 $b$-tag &\multicolumn{2}{c}{1+2 $b$-tag add. $b$-tag} \\ "+ "\n")
        f.write(r"    \hline"+ "\n")

        for each_mapper in namemapper:
            f.write("        " + each_mapper[1] + "    ")
            for each_tag in tagsegion:
                findtag = False
                for each_histo in histocollection:
                    if each_tag != each_histo.tag + each_histo.region:
                        continue
                    findtag = True
                    if each_tag != tagsegion[-1]:
                        if each_mapper[0] not in each_histo.bkgs:
                            f.write("  & 0    ")
                        else:
                            f.write("  & " + str(each_histo.bkgs[each_mapper[0]][0]) + r"$\pm$" + str(each_histo.bkgs[each_mapper[0]][1]))
                    else:
                        if each_mapper[0] not in each_histo.bkgs:
                            f.write(r"  & \multicolumn{2}{c}{0}")
                        else:
                            f.write(r"  & \multicolumn{2}{c}{ " + str(each_histo.bkgs[each_mapper[0]][0]) + r"$\pm$" + str(each_histo.bkgs[each_mapper[0]][1])+ r"}")
                if not findtag:
                    if each_tag != tagsegion[-1]:
                        f.write("  & 0    ")
                    else:
                        f.write(r"  & \multicolumn{2}{c}{0}")
            f.write(r" \\" + "\n")

        f.write(r"    \hline" + "\n")

        for each_mapper in ["Total", "Data"]:
            f.write("        " + each_mapper + "    ")
            for each_tag in tagsegion:
                findtag = False
                for each_histo in histocollection:
                    if each_tag != each_histo.tag + each_histo.region:
                        continue
                    findtag = True
                    if each_mapper == "Total":
                        if each_tag != tagsegion[-1]:
                            f.write("  & " + str(each_histo.totalbkg) + r"$\pm$" + str(each_histo.totalbkg_error))
                        else:
                            f.write(r"  & \multicolumn{2}{c}{ " + str(each_histo.totalbkg) + r"$\pm$" + str(each_histo.totalbkg_error) + r"}")
                    if each_mapper == "Data":
                        if each_tag != tagsegion[-1]:
                            f.write("  & " + str(each_histo.data) + r"$\pm$" + str(each_histo.data_error))
                        else:
                            f.write(r"  & \multicolumn{2}{c}{ " + str(each_histo.data) + r"$\pm$" + str(each_histo.data_error) + r"}")
                if not findtag:
                    if each_tag != tagsegion[-1]:
                        f.write("  & 0    ")
                    else:
                        f.write(r"  & \multicolumn{2}{c}{0}")
            f.write(r" \\" + "\n")


        f.write(r"    \hline\hline"+ "\n")
        f.write(r"    \end{tabular}"+ "\n")
        f.write(r"    \end{footnotesize}"+ "\n")
        f.write(r"    \end{center}"+ "\n")
        f.write(r"\end{table}"+ "\n")


    # ttbar     & xx$\pm$xx & xx$\pm$xx &xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & \multicolumn{2}{c}{xx$\pm$xx} \\
    # Single top & xx$\pm$xx & xx$\pm$xx &xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & \multicolumn{2}{c}{xx$\pm$xx} \\
    # Diboson    & xx$\pm$xx & xx$\pm$xx &xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & \multicolumn{2}{c}{xx$\pm$xx} \\
    # Zl        & xx$\pm$xx & xx$\pm$xx &xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & \multicolumn{2}{c}{xx$\pm$xx} \\
    # Zhl       & xx$\pm$xx & xx$\pm$xx &xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & \multicolumn{2}{c}{xx$\pm$xx} \\
    # Zhf       & xx$\pm$xx & xx$\pm$xx &xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & \multicolumn{2}{c}{xx$\pm$xx} \\
    # Wl        & xx$\pm$xx & xx$\pm$xx &xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & \multicolumn{2}{c}{xx$\pm$xx} \\
    # Whl       & xx$\pm$xx & xx$\pm$xx &xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & \multicolumn{2}{c}{xx$\pm$xx} \\
    # Whf       & xx$\pm$xx & xx$\pm$xx &xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & \multicolumn{2}{c}{xx$\pm$xx} \\
    # SM $Vh$    & xx$\pm$xx & xx$\pm$xx &xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & \multicolumn{2}{c}{xx$\pm$xx} \\
    # ttbar $h$ & xx$\pm$xx & xx$\pm$xx &xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & \multicolumn{2}{c}{xx$\pm$xx} \\
    # ttbar $V$ & xx$\pm$xx & xx$\pm$xx &xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & \multicolumn{2}{c}{xx$\pm$xx} \\
    # \hline
    # Total      & xx$\pm$xx & xx$\pm$xx &xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & xx$\pm$xx & \multicolumn{2}{c}{xx$\pm$xx} \\
    # Data       & xxxx & xxxx & xxxx & xxxx & xxxx & xxxx & \multicolumn{2}{c}{xxxx} \\
if __name__ == "__main__":
    main()