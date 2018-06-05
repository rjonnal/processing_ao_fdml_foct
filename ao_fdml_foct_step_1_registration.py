import numpy as np
import scipy.io as sio
from matplotlib import pyplot as plt
import sys,os
from cuttlefish import make_configuration, Hive, Series
from cuttlefish import utils #.utils import shear,find_peaks,ascend,descend
import scipy.optimize as soo
from fig2gif import GIF
import glob
from registration_references import ref_info

tag = sys.argv[1].strip('/').strip('\\')

oversample_factor = 7
strip_width = 5
refine = False
do_plot = False
ref_idx = 100
ref_vol = 0

refstr = ['','_r'][refine]

ref_idx,ref_vol = ref_info[tag][0]

all_files = glob.glob(os.path.join(tag,'0*'))
# get current working directory to convert relative paths to absolute paths
cwd = os.getcwd()
# horizontal (fast) coordinates to extract strips for registration and rendering

n_files = len(all_files)

try:
    idx1 = int(sys.argv[2])
    idx2 = int(sys.argv[3])
except Exception as e:
    idx1 = 0
    idx2 = n_files

# layer name for rendering only
# layer for registration is always CONES
try:
    layer_name = sys.argv[4]
except Exception as e:
    layer_name = 'CONES'
    
ref = np.load('%s/%04d/projections/CONES.npy'%(tag,ref_idx))[ref_vol,:,:]
    
flist = [os.path.join(cwd,os.path.join('%s'%tag,'%04d'%k)) for k in range(idx1,idx2)]

# strip register
s = Series('%s_registered_ref_%04d_%02d_%dx_%dw%s'%(tag,ref_idx,ref_vol,oversample_factor,strip_width,refstr),ref)

s.hive.put('oversample_factor',oversample_factor)
s.hive.put('strip_width',strip_width)
s.hive.put('refine',int(refine))

if not s.is_registered()>=n_files:
    for f in flist:
        for vidx in range(2):
            s.add(f,vidx,layer_names=['CONES'],strip_width=strip_width,do_plot=do_plot,refine=refine,oversample_factor=oversample_factor)
else:
    s.goodness_histogram()
    s.render(layer_names=[layer_name],goodness_threshold=0.4,do_plot=do_plot,left_crop=10,oversample_factor=oversample_factor)
