import logging

import ROOT


import numpy as np

ltx = ROOT.TLatex()

textsize = 0.06

# start_y = 0.92
start_y = 0.87
# start_y = 0.85
# step_y = 0.05
step_y = 0.075

# x_pos = 0.225
x_pos = 0.2
# x_pos = 0.125
# x_pos = 0.1
y_pos = start_y

logger = logging.getLogger(__name__)


def asymptotic_sensitivity(s, b, sigma):
    """
    Asymptotic sensitivity formula as specified in https://cds.cern.ch/record/2643488?

    Args:
        s: Number of signal events
        b: Number of background events
        sigma: Background uncertainty

    Returns:
        Z: The significance as calculated with formula (1) or (2) from the paper
    """
    if s == 0:
        raise ValueError("Function ill-defined for 0 expected signal events.")
    if b == 0:
        logger.debug("0 background events expected: Will use n_background = 1e-6 for calculation")
        b = 1e-6
    n = s + b
    if sigma == 0:
        sig = np.sqrt(2 * (n * np.log(n / b) - (n - b)))
        if n >= b:
            Z = sig
        else:
            Z = -sig
    else:
        term1 = n * (b + sigma ** 2)
        term1 /= b ** 2 + n * sigma ** 2
        term1 = np.log(term1)
        term1 *= n
        term2 = sigma ** 2 * (n - b)
        term2 /= b * (b + sigma ** 2)
        term2 += 1
        term2 = np.log(term2)
        term2 *= b ** 2 / sigma ** 2
        sig = np.sqrt(2 * (term1 - term2))
        if n >= b:
            Z = sig
        else:
            Z = -sig
    logger.debug(
        "asymptotic_sensitivity: s = %.2f, b = %.2f, sigma = %.2f; Z = %.2f", s, b, sigma, Z,
    )
    if np.isnan(Z):
        return 0
    return Z


def remove_ownership(obj):
    """
    Remove ROOT ownership of an object

    If you apply this method with e.g. a canvas you will be able to return it from
    a function without losing the plot contents.
    """
    for prim in obj.GetListOfPrimitives():
        if isinstance(prim, ROOT.TPad):
            ROOT.SetOwnership(prim, False)
            remove_ownership(prim)
        ROOT.SetOwnership(prim, False)


def ATLASSimulationInternal(x=None, y=None, size=textsize):
    """Draw on ROOT.gPad"""
    if not x:
        global x_pos
        x = x_pos
    if not y:
        global y_pos
        y = y_pos
        y_pos -= step_y * size / size
    ltx.SetTextSize(size)
    ltx.SetTextFont(42)
    ltx.DrawLatexNDC(x, y, "#bf{#it{ATLAS}} Simulation Internal")
    y_pos -= step_y


def ATLASInternal(x=None, y=None, size=textsize):
    """Draw on ROOT.gPad"""
    if not x:
        global x_pos
        x = x_pos
    if not y:
        global y_pos
        y = y_pos
        y_pos -= step_y * size / size
    ltx.SetTextSize(size)
    ltx.SetTextFont(42)
    # print('DrawLatexNDC({}, {}, #bf{{#it{{ATLAS}}}} Internal)'.format(x, y))
    ltx.DrawLatexNDC(x, y, "#bf{#it{ATLAS}} Internal")


def ATLASWIP():
    """Draw on ROOT.gPad"""
    global y_pos
    ltx.SetTextSize(textsize)
    ltx.SetTextFont(42)
    ltx.DrawLatexNDC(x_pos, y_pos, "#bf{#it{ATLAS}} Work In Progress")
    y_pos -= step_y


def ATLASSimWIP(x=None, y=None, size=textsize):
    """Draw on ROOT.gPad"""
    if not x:
        global x_pos
        x = x_pos
    if not y:
        global y_pos
        y = y_pos
        y_pos -= step_y * size / size
    ltx.SetTextSize(size)
    ltx.SetTextFont(42)
    ltx.DrawLatexNDC(x, y, "#bf{#it{ATLAS}} Simulation Work In Progress")


def Energy(energy=13, x=None, y=None, size=textsize):
    """Draw on ROOT.gPad"""
    if not x:
        global x_pos
        x = x_pos
    if not y:
        global y_pos
        y = y_pos
        y_pos -= step_y * size / size
    ltx.DrawLatexNDC(x, y, "#sqrt{{s}} = {} TeV".format(energy))
    y_pos -= step_y


def EnergyAndLumi(energy=13, lumi=36.1, x=None, y=None, size=textsize):
    """Draw on ROOT.gPad"""
    global y_pos
    ltx.SetTextSize(size)
    if not x:
        x = x_pos
    if not y:
        y = y_pos
        y_pos -= step_y * size / size
    # ltx.DrawLatexNDC(
    #     x, y, "#sqrt{{s}} = {} TeV, #scale[0.65]{{#int}} L = {} fb^{{-1}}".format(energy, lumi),
    # )
    ltx.DrawLatexNDC(
        x, y, "#sqrt{{s}} = {} TeV, {} fb^{{-1}}".format(energy, lumi),
    )


def DrawText(text, x=None, y=None, size=textsize):
    global y_pos
    if not x:
        x = x_pos
    if not y:
        y = y_pos
        y_pos -= step_y * size / textsize
    ltx.SetTextSize(size)
    obj = ltx.DrawLatexNDC(x, y, text)
    ltx.SetTextSize(textsize)
    return obj


def set_x_pos(pos):
    global x_pos
    x_pos = pos


def set_y_pos(pos):
    global y_pos
    y_pos = pos


def reset_y_pos():
    global y_pos
    y_pos = start_y
