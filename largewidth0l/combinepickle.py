import pickle
import os


paras = {}
allfiles = os.listdir("pickle")
for each in allfiles:
    mass = int(each.replace("fitvalues", "").replace(".pickle", ""))
    with open("pickle/"  + each,'rb') as f:
        pvalues = pickle.load(f)
        paras[mass] = pvalues
with open('fitvalues.pickle', 'wb') as f:
    pickle.dump(paras, f)