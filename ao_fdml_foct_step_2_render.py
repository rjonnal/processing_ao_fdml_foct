import numpy as np
import scipy.io as sio
from matplotlib import pyplot as plt
import sys,os
from cuttlefish import make_configuration, Hive, Series
from cuttlefish import utils #.utils import shear,find_peaks,ascend,descend
import scipy.optimize as soo
from fig2gif import GIF
import glob

goodness_threshold = 2.0
layer_names = ['CONES']

try:
    tag = sys.argv[1].strip('/').strip('\\')
except Exception as e:
    print e
    tag = 'foo'
    
do_plot = True

s = Series(tag)
oversample_factor = s.hive.get('oversample_factor')[0]

s.goodness_histogram(goodness_threshold)

s.render(layer_names=layer_names,goodness_threshold=goodness_threshold,do_plot=do_plot,oversample_factor=oversample_factor)
