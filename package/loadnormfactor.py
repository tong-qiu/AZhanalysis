import re
import yaml

def loadnorm(configpath, normpath):
    sysdic = {}
    with open(configpath) as f:
        for eachline in f:
            eachline = eachline.replace(" ", "")
            eachline = eachline.replace("\n", "")
            eachline = eachline.split("#")[0]
            if len(eachline) < 1:
                continue
            ss1 = eachline.split("=")
            if "FloatSyst" not in ss1[0] and "NormSyst" not in ss1[0] and "Luminosity" not in ss1[0]:
                continue
            ss2 = ss1[1].split(":")
            #print(ss2)
            if ss2[2] in sysdic:
                sysdic[ss2[2]].append([ss2[0], ss2[1], float(ss2[3]), float(ss2[4])])
            else:
                sysdic[ss2[2]] = [[ss2[0], ss2[1], float(ss2[3]), float(ss2[4])]]
    #print(sysdic)
    # key1: sysname, key2: region, key3: value
    norm_table = {}
    norm_table_witherror = {}
    with open(normpath) as f:
        for eachline in f:
            if "CORRELATION_MATRIX" in eachline:
                break
            eachline = eachline.replace("\n", "")
            eachline = eachline.replace("\\", "")
            eachline = eachline.replace(" ", "")
            match = re.search(r"(.+)\&\$([\-e0-9.]+)\^\{([e\-\+0-9.]+)\}\_\{([e\+\-0-9.]+)", eachline)
            if match is not None and match.group(1) in sysdic:
                for each in sysdic[match.group(1)]:
                    print(each)
                    if "norm" in match.group(1):
                        sys_tem = float(match.group(2)) - 1.
                        sys_tem_error = float(match.group(3))
                    else:
                        sys_tem = float(match.group(2)) * float(each[2])
                        sys_tem_error = float(match.group(3)) * float(each[2])
                        # sys_tem = 0.
                        # if float(match.group(2)) > 0:
                        #     sys_tem = float(match.group(2)) * float(each[2])
                        # if float(match.group(2)) < 0:
                        #     sys_tem = float(match.group(2)) * float(each[2])

                    if each[1] not in norm_table:
                        norm_table[each[1]] = {}
                        norm_table[each[1]][each[0]] = sys_tem
                    else:
                        if each[0] in norm_table[each[1]]:
                            norm_table[each[1]][each[0]] += sys_tem
                        else:
                            norm_table[each[1]][each[0]] = sys_tem

                    if each[1] not in norm_table_witherror:
                        norm_table_witherror[each[1]] = {}
                        norm_table_witherror[each[1]][each[0]] = [sys_tem + 1, sys_tem_error]
                    else:
                        if each[0] in norm_table_witherror[each[1]]:
                            norm_table_witherror[each[1]][each[0]][0] += sys_tem
                            norm_table_witherror[each[1]][each[0]][1] = (norm_table_witherror[each[1]][each[0]][1]**2 + sys_tem_error**2)**0.5
                        else:
                            norm_table_witherror[each[1]][each[0]] = [sys_tem + 1, sys_tem_error]
        print(norm_table_witherror)
        #print(norm_table)
        return norm_table

                # print(match.group(1))
                # print(match.group(2))
                # print(match.group(3))
                # print(float(match.group(4)))
                #break

def generatejson(inputs):
    resolved = {}
    merged = {}
    for eachbkg in inputs.keys():
        for eachreg in inputs[eachbkg].keys():
            if eachreg == "ALL":
                resolved[eachbkg] = inputs[eachbkg][eachreg] + 1
            if eachbkg not in merged:
                merged[eachbkg] = inputs[eachbkg][eachreg] + 1
            else:
                merged[eachbkg] += inputs[eachbkg][eachreg]
    dp = {"resolved": resolved, "merged": merged}
    print(yaml.dump(dp))


if __name__ == "__main__":
    # loadnorm("C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/configLLBB_190517_HVT_PRSR_MCstat0_Prun1_finalNPtreatment_RFfixC0_2000.cfg",
    #  "C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/GlobalFit_fitres_unconditionnal_mu0.txt")
    # rescaledic = loadnorm("C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/confignormonly.cfg",
    # "C:/Users/qiutt/Desktop/postreader/PlotTool_Root/jsonoutput/GlobalFit_fitres_unconditionnal_mu0_normonly.txt")
    # rescaledic = loadnorm("C:/Users/qiutt/Desktop/postreader/pullandcorrelation/2hdm_norm/config.cfg",
    # "C:/Users/qiutt/Desktop/postreader/pullandcorrelation/2hdm_norm/GlobalFit_fitres_unconditionnal_mu0.txt")
    rescaledic = loadnorm("../2HDMNPs/config_m300.cfg",
    "../2HDMNPs/GlobalFit_fitres_conditionnal_mu0.txt")
    generatejson(rescaledic)
