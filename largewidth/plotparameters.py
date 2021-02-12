import pickle
import numpy as np
import matplotlib.pyplot as plt

with open("fitvalues_merged.pickle",'rb') as f:
    values = pickle.load(f)
    print(values)

mass = []
x = []
for each_mass in sorted(values.keys()):
    mass.append(each_mass)
    x.append(values[each_mass][20][0][0])

print(values)
plt.plot(mass, x)
plt.show()

