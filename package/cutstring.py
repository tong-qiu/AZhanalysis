import numpy as np
import zlib

def s_resolved(mata):
    mask = mata["Regime"] == zlib.adler32(b'resolved')
    return mask
def s_merged(mata):
    mask = mata["Regime"] == zlib.adler32(b'merged')
    return mask
def s_mbbcr(mata):
    mask = mata["Description"] == zlib.adler32(b'mBBcr')
    return mask
def s_sr(mata):
    mask = mata["Description"] == zlib.adler32(b'SR')
    return mask

def s_zhf(mata):
    mask = mata["Sample"] == zlib.adler32(b'Zbb')
    mask = np.logical_or(mata["Sample"] == zlib.adler32(b'Zbc'), mask)
    mask = np.logical_or(mata["Sample"] == zlib.adler32(b'Zcc'), mask)
    return mask

def s_zlf(mata):
    mask = mata["Sample"] == zlib.adler32(b'Zbl')
    mask = np.logical_or(mata["Sample"] == zlib.adler32(b'Zcl'), mask)
    mask = np.logical_or(mata["Sample"] == zlib.adler32(b'Zl'), mask)
    return mask