import os

import uproot

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pickle
from matplotlib import rc
import sys
import zlib
import re
import copy

lib_path = os.path.abspath(os.path.join(__file__, '..', '..'))
sys.path.append(lib_path)
sys.path.append('../package')

from package.events import *
from mlcut import *
from package.stackplot import *
from curveplot import *
from cutstring import *
import multiprocessing

with open("background.pickle",'rb') as f:
    sample = pickle.load(f)