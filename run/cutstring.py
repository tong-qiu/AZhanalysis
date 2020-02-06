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