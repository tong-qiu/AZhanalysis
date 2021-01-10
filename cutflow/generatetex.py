import json
import os
def getcutname(inputs):
    # print(inputs)
    returnlist = []
    for each in inputs:
        temname = each
        if each == "mll":
            temname = r"$m_{ll}$"
        if each == "mBB":
            temname = r"$m_{BB}$"
        if each == "ptV":
            temname = r"$P_{T}^{ll}$"
        if each == "pTB1":
            temname = r"$P_{T}^{B1}$"
        if r">=" in each:
            temname = "at least" + temname[2:]
        returnlist.append(temname)
    return returnlist

def createlist(cuts, cutstex, rmkey, samplename, datadic):
    tem_samplename = []
    for eachsamplekey in samplename:
        if rmkey in datadic[eachsamplekey]:
            tem_samplename.append(eachsamplekey)
    samplename = tem_samplename
            
    with open("table.tex", "a") as f:
        f.write(r"\begin{table}[H]" + "\n")
        f.write(r"    \begin{center}" + "\n")
        f.write(r"        \begin{footnotesize}" + "\n")
        f.write(r"            \begin{tabular}{l|" + " c|" * len(samplename) + "}" + "\n")
        f.write(r"                \hline\hline" + "\n")

        f.write(r"                    cut  ")
        for each in samplename:
            theoutput = each
            if each == "ttbar":
                theoutput = r"$t\bar{t}$"
            if "HVT" in each:
                theoutput = each.replace("HVT", r"Z$^\prime$ ") + " GeV"

            f.write("   & " + theoutput)
        f.write(r"    \\" + "\n")

        f.write(r"                \hline\hline" + "\n")

        i = 0
        for eachselectionkey, eachkeytex in zip(cuts, cutstex):
            i += 1
            f.write(r"                    " + eachkeytex + "  ")
            for eachsamplekey in samplename:
                if rmkey not in datadic[eachsamplekey]:
                    print(eachsamplekey)
                    continue

                for eachvaluelist in datadic[eachsamplekey][rmkey]:
                    if eachvaluelist[0] == eachselectionkey:
                        f.write("   & " + str(int(eachvaluelist[1])))
                        break
            if i != len(cuts):
                f.write(r" \\" + " \n")
        f.write(r"    \\" + "\n")

        f.write(r"                \hline\hline" + "\n")
        f.write(r"            \end{tabular}" + "\n")
        f.write(r"        \end{footnotesize}" + "\n")
        f.write(r"    \end{center}" + "\n")
        f.write(r"    \caption{The cutflow table in " + rmkey + " regime. The values are the number of unweighted MC events after each cut.} " + "\n")
        f.write(r"    {\label{tab:2leventselectiontable" + rmkey + r"}}" + "\n")
        f.write(r"\end{table}" + "\n\n\n\n\n")

def main():
    with open("result/cutflowlen.json", 'r') as f:
        for each in f:
            datadic = json.loads(each)
            break
    #print(datadic)
    os.remove("table.tex")

    samplename = [*datadic.keys()]
    resolvedcut = [each[0] for each in datadic[samplename[0]]['resolved']]
    mergedcut = [each[0] for each in datadic[samplename[0]]['merged']]
    print(resolvedcut)
    resolvedcuttex = getcutname(resolvedcut)
    mergedcuttex = getcutname(mergedcut)
    print(resolvedcut)
    createlist(resolvedcut, resolvedcuttex, "resolved", samplename, datadic)
    createlist(mergedcut, mergedcuttex, "merged", samplename, datadic)



if __name__ == "__main__":
    main()
