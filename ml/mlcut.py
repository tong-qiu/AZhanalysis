import numpy as np

def cut_btag_more(data, b):
    mask = data[b'nbJets'] > b
    return mask

def cut_btag_is(data, b):
    #mask = data[b"nSigJets"] >= 2
    mask = data[b'nbJets'] == b
    return mask

def cut_btag(data):
    #mask = data[b"nSigJets"] >= 2
    mask = data[b'nbJets'] == 0
    return mask

def cut_basic(data):
    # two singal jets
    mask = data[b"nSigJets"] >= 2 # essential
    mask = np.logical_and(data[b'ptL1']/1000. > 27, mask) # essential
    mask = np.logical_and(data[b'ptL2']/1000. > 20, mask) # essential
    return mask

@np.vectorize
def ptll_cut(mvh):
    if mvh < 320.:
        return 0.
    return 20. + 9. * pow(mvh - 320., 0.5)

def srcut(data):
    mask = data[b"nSigJets"] >= 2 # essential
    #mask = np.logical_and(data[b'nbJets'] < 3, mask)
    #mask = np.logical_and(data[b'nbJets'] == 2, mask)
    mask = np.logical_and(data[b"passedTrigger"] == 1, mask)
    mask = np.logical_and(data[b'flavL1'] == data[b'flavL2'], mask)
    mask = np.logical_and(np.logical_or(data[b'flavL1'] == 1, data[b'chargeL1'] != data[b'chargeL2']), mask)
    #mask =np.logical_and(np.logical_or(data[b'flavL1'] == 1, abs(data[b'etaL1']) < 2.5), mask)
    # cut on mll
    lower_limit = [max(40, 87 - 0.03 * each) for each in (data[b"mVH"]/1000.)]
    higher_limit = 97 + 0.013 * data[b"mVH"]/1000.
    mask = np.logical_and(lower_limit < data[b"mLL"]/1000., mask) # essential
    mask = np.logical_and(data[b"mLL"]/1000. < higher_limit, mask)

    mask = np.logical_and(data[b"METHT"]/(1000.**0.5) < 1.15 + 8 * (10**(-3))*data[b"mVH"]/1000., mask)
    mask = np.logical_and(data[b'pTV']/1000. > ptll_cut(data[b'mVH']/1000.), mask)
    mask = np.logical_and(data[b'pTB1']/1000. > 45, mask) # essential
    mask = np.logical_and(data[b'mBBres']/1000. < 145, mask)
    mask = np.logical_and(100 < data[b'mBBres']/1000., mask)
    mask = np.logical_and(data[b'ptL1']/1000. > 27, mask) # essential
    mask = np.logical_and(data[b'ptL2']/1000. > 20, mask) # essential
    #mask = np.logical_and(data[b'ptL2']/1000. > 20, mask)
    return mask

def crmbbcut(data):
    mask = data[b"nSigJets"] >= 2 # essential
    #mask = np.logical_and(data[b'nbJets'] >= 1, mask)
    mask = np.logical_and(data[b"passedTrigger"] == 1, mask)
    mask = np.logical_and(data[b'flavL1'] == data[b'flavL2'], mask)
    mask = np.logical_and(np.logical_or(data[b'flavL1'] == 1, data[b'chargeL1'] != data[b'chargeL2']), mask)
    #mask =np.logical_and(np.logical_or(data[b'flavL1'] == 0, abs(data[b'etaL1']) < 2.5), mask)
    # cut on mll
    lower_limit = [max(40, 87 - 0.03 * each) for each in (data[b"mVH"]/1000.)]
    higher_limit = 97 + 0.013 * data[b"mVH"]/1000.
    mask = np.logical_and(lower_limit < data[b"mLL"]/1000., mask) # essential
    mask = np.logical_and(data[b"mLL"]/1000. < higher_limit, mask)

    mask = np.logical_and(data[b"METHT"]/(1000.**0.5) < 1.15 + 8 * (10**(-3))*data[b"mVH"]/1000., mask)
    mask = np.logical_and(data[b'pTV']/1000. > ptll_cut(data[b'mVH']/1000.), mask)
    mask = np.logical_and(data[b'pTB1']/1000. > 45, mask) # essential
    mask = np.logical_and(data[b'ptL1']/1000. > 27, mask) # essential
    mask = np.logical_and(data[b'ptL2']/1000. > 20, mask) # essential
    mask = np.logical_and(np.logical_or(data[b'mBBres']/1000. >= 145, 100 >= data[b'mBBres']/1000.), mask)
    mask = np.logical_and(data[b'mBBres']/1000. < 200, mask)
    mask = np.logical_and(50 < data[b'mBBres']/1000., mask)
    return mask


def crtopcut(data):
    mask = data[b"nSigJets"] >= 2 # essential
    mask = np.logical_and(data[b'nbJets'] < 2, mask)
    mask = np.logical_and(data[b"passedTrigger"] == 1, mask)
    mask = np.logical_and(data[b'flavL1'] != data[b'flavL2'], mask)
    mask = np.logical_and(data[b'chargeL1'] != data[b'chargeL2'], mask)
    #mask =np.logical_and(np.logical_or(data[b'flavL1'] == 0, abs(data[b'etaL1']) < 2.5), mask)
    # cut on mll
    lower_limit = [max(40, 87 - 0.03 * each) for each in (data[b"mVH"]/1000.)]
    higher_limit = 97 + 0.013 * data[b"mVH"]/1000.
    mask = np.logical_and(lower_limit < data[b"mLL"]/1000., mask) # essential
    mask = np.logical_and(data[b"mLL"]/1000. < higher_limit, mask)
    mask = np.logical_and(data[b'pTV']/1000. > ptll_cut(data[b'mVH']/1000.), mask)
    mask = np.logical_and(data[b'pTB1']/1000. > 45, mask) # essential
    mask = np.logical_and(data[b'mBBres']/1000. < 145, mask)
    mask = np.logical_and(100 < data[b'mBBres']/1000., mask)
    mask = np.logical_and(data[b'ptL1']/1000. > 27, mask) # essential
    mask = np.logical_and(data[b'ptL2']/1000. > 20, mask) # essential
    return mask


def cut_mcid(mccn):
    def return_cut(data):
        mask = data[b'MCChannelNumber'] == mccn
        return mask
    return return_cut
    