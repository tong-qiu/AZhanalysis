import matplotlib.pyplot as plt
import ROOT
import math
def bw(x, p0, p1):
    # gamma = ((p1**2) * (p1**2 + p0**2))**0.5
    # k = 2 * 2**0.5 + p1 * p0 * gamma / (3.14*(p1**2 + gamma)**0.5)

    return 1/((x - p1)**2 + 0.25*p0*p0)

# for m in [100, 300, 500, 1000]:
x = range(-2000, 2000)
y = []
for mi in x:
    y.append(ROOT.Math.lognormal_pdf(mi, abs(math.log(100)), math.log(1.1), -2000))
    # y.append(bw(mi, 100, m))
plt.plot(x, y)
plt.show()
