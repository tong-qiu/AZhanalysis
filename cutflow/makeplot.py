import matplotlib.pyplot as plt
import numpy as np
import json
from itertools import cycle

def barplot(height, tick, filename, legend):
    x = range(len(tick)+1)
    cycol = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k','yellowgreen'])
    for each, eachlegend in zip(height, legend):
        each1 = [each[0]] + each
        each1 = np.array(each1)/each1[0]
        plt.step(x, each1, color=next(cycol), label=eachlegend)
    plt.legend(prop={'size': 10}, frameon=False, ncol=2)
    x = range(len(tick))
    plt.xlim((0,len(tick)))
    plt.ylim((0,1.4))
    plt.xticks(np.array(x)+0.5, tick, fontsize=10)
    plt.xticks(rotation=45, rotation_mode="anchor", ha='right',)
    plt.ylabel("Relative yield after applying given cuts", fontsize=10)
    plt.savefig(filename + '.pdf', bbox_inches='tight', pad_inches = 0.2,)
    plt.show()

if __name__ == "__main__":
    # tick = ["asjfosdf fsf", "fsdfjoid fdsf", "jfodjfoisjfi fdsa", "fjdfjod"]
    # y = ([1000, 300, 100,10], [1300, 400, 30,40]) 
    # legend = ["121","1232"]
    # colors = ["r", "blue"]
    # barplot(y,tick,colors,"test", legend)
    # samples = {}
    # exit(1)
    with open("result/cutflowlen.json", 'r') as f:
        for each in f:
            samples = json.loads(each)
            break
    legend = []
    legend_merged = []
    resolvedheight = []
    resolvedcut = []
    mergedheight = []
    mergedcut = []
    for eachsample in samples.keys():
        legend.append(eachsample)
        resolved = samples[eachsample]['resolved']
        resolvedheight_tem = []
        resolvedcut = []
        for each in resolved:
            resolvedcut.append(each[0])
            resolvedheight_tem.append(each[1])
        resolvedheight.append(resolvedheight_tem)
        if "merged" in samples[eachsample]:
            legend_merged.append(eachsample)
            mergedheight_tem = []
            mergedcut = []
            for each1 in samples[eachsample]['merged']:
                mergedcut.append(each1[0])
                mergedheight_tem.append(each1[1])
            mergedheight.append(mergedheight_tem)
    barplot(resolvedheight,resolvedcut,"resolved", legend)
    barplot(mergedheight,mergedcut,"merged", legend_merged)