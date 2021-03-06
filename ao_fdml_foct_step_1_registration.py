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

oversample_factor = 3
strip_width = 9
# ref_idx/ref_vol values: for DataSet1, 89/0; for 2018.05.09_10.50.28_RE_3DT_20uw-320x-PROC-5Reps-0.3_PVthresh, 50
ref_idx = 50
ref_vol = 0
do_plot = False
left_crop = 45
right_crop = 35

#For Dataset1
#left_crop = None
#right_crop = None


#oversample_factor = 7
#strip_width = 7
#ref_idx = 173
#ref_vol = 0
#do_plot = False
#left_crop = 45
#right_crop = 35

all_files = glob.glob(os.path.join(tag,'volumes/*'))
all_files.sort()


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

ref = np.load('%s/volumes/%04d/projections/CONES.npy'%(tag,ref_idx))[ref_vol,:,:]

flist = all_files[idx1:idx2]

# strip register
s = Series('%s_registered_ref_%04d_%02d_%dx_%dw'%(tag,ref_idx,ref_vol,oversample_factor,strip_width),ref)

s.hive.put('oversample_factor',oversample_factor)
s.hive.put('strip_width',strip_width)
s.hive.put('reference_index',ref_idx)
s.hive.put('reference_volume',ref_vol)
s.hive.put('left_crop',left_crop)
s.hive.put('right_crop',right_crop)


for f in flist:
    for vidx in range(2):
        s.add(f,vidx,layer_names=['CONES'],strip_width=strip_width,do_plot=do_plot,oversample_factor=oversample_factor,left_crop=left_crop,right_crop=right_crop)
