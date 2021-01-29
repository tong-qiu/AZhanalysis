import ROOT
import matplotlib.pyplot as plt
import numpy as np

gen = ROOT.TRandom2()
datas = []
for each in range(0, 100000):
    datas.append(gen.BreitWigner (100, 10))

fig, ax = plt.subplots(figsize=(10,8))
ax.hist(np.array(datas), np.linspace(0,200,1000), histtype='step', fill=False)
fig.savefig("test2.pdf")
