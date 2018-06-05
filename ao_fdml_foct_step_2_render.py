import numpy as np
import scipy.io as sio
from matplotlib import pyplot as plt
import sys,os
from cuttlefish import make_configuration, Hive, Series
from cuttlefish import utils #.utils import shear,find_peaks,ascend,descend
import scipy.optimize as soo
from fig2gif import GIF
import glob

try:
    tag = sys.argv[1].strip('/').strip('\\')
except Exception as e:
    print e
    tag = 'foo'
    
do_plot = True

s = Series(tag)
oversample_factor = s.hive.get('oversample_factor')[0]

s.goodness_histogram()
s.render(layer_names=['CONES'],goodness_threshold=2.0,do_plot=do_plot,left_crop=10,oversample_factor=oversample_factor)
