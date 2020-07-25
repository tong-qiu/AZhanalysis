import ROOT
import math
def fitfunction_root(x, par):
    # source: https://github.com/root-project/root/blob/master/roofit/roofit/src/RooBukinPdf.cxx
    Xp = par[0]
    sigp = par[1]
    xi = par[2]
    rho1 = par[3]
    rho2 = par[4]
    ap = par[5]
    consts = 2*math.sqrt(2*math.log(2.0))
    r1=0
    r2=0
    r3=0
    r4=0
    r5=0
    hp=0
    x1 = 0
    x2 = 0
    fit_result = 0

    hp=sigp*consts
    r3=ROOT.TMath.Log(2.)
    r4=math.sqrt(xi*xi+1)
    r1=xi/r4

    if abs(xi) > math.exp(-6.):
        r5=xi/ROOT.TMath.Log(r4+xi)
    else:
        r5=1

    x1 = Xp + (hp / 2) * (r1-1)
    x2 = Xp + (hp / 2) * (r1+1)

    #--- Left Side
    if x[0] < x1:
        r2=rho1*(x[0] - x1) * (x[0] - x1)/(Xp - x1)/(Xp - x1)-r3 + 4 * r3 * (x[0] - x1)/hp * r5 * r4/(r4 - xi)/(r4 - xi)
    #--- Center
    elif x[0] < x2:
        if abs(xi) > math.exp(-6.):
            r2=ROOT.TMath.Log(1 + 4 * xi * r4 * (x[0] - Xp)/hp)/ROOT.TMath.Log(1 + 2*xi*(xi - r4))
            r2=-r3*r2*r2
        else:
            r2=-4*r3*(x[0] - Xp)*(x[0] - Xp)/hp/hp
    #--- Right Side
    else:
        r2=rho2*(x[0] - x2)*(x[0] - x2)/(Xp - x2)/(Xp - x2)-r3 - 4 * r3 * (x[0] - x2)/hp * r5 * r4/(r4 + xi)/(r4 + xi)

    if abs(r2) > 100:
        fit_result = 0
    else:
        #---- Normalize the result
        fit_result = math.exp(r2)
    return fit_result *ap