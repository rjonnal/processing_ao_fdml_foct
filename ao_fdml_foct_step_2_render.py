import numpy as np
import scipy.io as sio
from matplotlib import pyplot as plt
import sys,os
from cuttlefish import make_configuration, Hive, Series
from cuttlefish import utils #.utils import shear,find_peaks,ascend,descend
import scipy.optimize as soo
from fig2gif import GIF
import glob


left_crop = 40
right_crop = 45

layer_names = ['CONES']

try:
    tag = sys.argv[1].strip('/').strip('\\')
except Exception as e:
    sys.exit('Please pass registered dataset name as an argument.')

try:
    goodness_threshold = float(sys.argv[2])
except:
    goodness_threshold = None
    
do_plot = True

s = Series(tag)
oversample_factor = s.hive.get('oversample_factor')[0]

if goodness_threshold is None:
    s.goodness_histogram()
else:
    s.render(layer_names=layer_names,goodness_threshold=goodness_threshold,do_plot=do_plot,oversample_factor=oversample_factor,left_crop=left_crop,right_crop=right_crop)
