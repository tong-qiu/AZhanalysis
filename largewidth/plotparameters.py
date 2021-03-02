import pickle
import numpy as np
import matplotlib.pyplot as plt

with open("fitvalues.pickle",'rb') as f:
    values = pickle.load(f)
    print(values[360])

mass = []
x = []
for each_mass in sorted(values.keys()):
    mass.append(each_mass)
    x.append(values[each_mass][20][0][0])
    # print(each_mass, values[each_mass][20][3][0])
#print(values)
plt.plot(mass, x)
plt.show()

