import numpy as np
import scipy.io as spi
from matplotlib import pyplot as plt
import sys,os
from cuttlefish import make_configuration, Hive, Series
from cuttlefish import utils #.utils import shear,find_peaks,ascend,descend
import scipy.ndimage as spn
from fig2gif import GIF
import glob
import cProfile
import time
from registration_references import ref_info
from cone_analysis_tools import analyze_cone

tag = sys.argv[1].strip('/').strip('\\')
all_files = glob.glob(os.path.join(tag,'frames/*'))
all_files.sort()
n_files = len(all_files)

s = Series(tag)
oversample_factor = s.hive.get('oversample_factor')[0]

keys = sorted(s.db.dictionary.keys())

try:
    idx1 = int(sys.argv[2])
    idx2 = int(sys.argv[3])
except Exception as e:
    idx1 = 0
    idx2 = n_files
 
series_hive = s.hive
layers = ['CONES']

efp = series_hive.get('average_image/CONES')
efpsy,efpsx = efp.shape
efp[np.where(efp==0)]=efp[np.where(efp)].mean()
cone_diameter = int(round(2*oversample_factor))

cones_x,cones_y = utils.find_cones(efp,cone_diameter,do_plot=True)

cones_x = np.round(cones_x).astype(np.integer)
cones_y = np.round(cones_y).astype(np.integer)

cone_rad = cone_diameter//2
region_rad = cone_rad*2

def drawbox(x,y):
    x1 = x-cone_rad
    x2 = x+cone_rad
    y1 = y-cone_rad
    y2 = y+cone_rad
    plt.plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],'r-')

cone_dictionary = {}

test_fn = os.path.split(os.path.split(s.db.get(keys[0])[0])[0])[0]
config = Hive('%s/config'%test_fn)
n_depth = config.get('n_depth')[0]

n_cones = len(cones_x)

cropped_depth = 100 # this is more than enough to get from ELM to RPE
intensity_cube = np.zeros((n_files,cropped_depth,n_cones))
time_cube = np.zeros((n_files,n_cones))

tidx = 0

running_sum = np.zeros(efp.shape)
running_count = np.ones(efp.shape)*1e-16

stimulus_time = 1.96

for key in keys[idx1:idx2]:#[89*2:]:
    try:
        foo = bar
        intensity_sheet = series_hive['frames/%s/cones'%key].get('intensity')
        phase_sheet = series_hive['frames/%s/cones'%key].get('phase')
        time_sheet = series_hive['frames/%s/cones'%key].get('time')
    except Exception as e:
        intensity_sheet = np.ones((cropped_depth,n_cones))*np.nan
        phase_sheet = np.ones((cropped_depth,n_cones))*np.nan
        time_sheet = np.ones((n_cones))*np.nan

        f,vidx = s.db.get(key)
        
        target_hive = Hive(f)
        label_dictionary = {}
        label_keys = target_hive['model/labels'].keys()
        for lk in label_keys:
            label_dictionary[lk] = target_hive['model/labels/%s'%lk][vidx]

        isos_idx = label_dictionary['ISOS']
        cost_idx = label_dictionary['COST']
        vol,t = s.map_into_reference_space(f,vidx,layers,goodness_threshold=0.4,oversample_factor=oversample_factor)
        n_slow,n_depth,n_fast = vol.shape

        model = target_hive['model/profiles'][vidx]
        this_efp = series_hive['frames/%s/projections'%key].get('CONES')

        average_interval = 20
        if (tidx+1)%average_interval:
            temp = this_efp.copy()
            temp[np.where(np.isnan(this_efp))] = 0.0
            running_sum = running_sum + temp
            running_count[np.where(temp)] = running_count[np.where(temp)] + 1.0
        else:
            running_average = running_sum/running_count
            running_sum = np.zeros(efp.shape)
            running_count = np.ones(efp.shape)*1e-16
            interval_time = np.round(np.nanmean(t)*1000)

        if False:
            XX,YY = np.meshgrid(np.arange(n_fast)-n_fast//2,np.arange(n_slow)-n_slow//2)
            sigma = 0.5
            g = np.exp(-(XX**2+YY**2)/2.0/sigma**2)
            this_efp = np.abs(np.fft.fftshift(np.fft.ifft2(np.fft.fft2(this_efp)*np.conj(np.fft.fft2(g)))))

        for cone_idx,(cx,cy) in enumerate(zip(cones_x,cones_y)):

            # first, cross_correlate neighborhood with efp to see if there's any local shift
            ry1 = cy-region_rad
            ry2 = cy+region_rad+1
            rx1 = cx-region_rad
            rx2 = cx+region_rad+1
            
            ry1 = max(ry1,0)
            ry2 = min(ry2,n_slow-1)
            rx1 = max(rx1,0)
            rx2 = min(rx2,n_fast-1)

            this_region = this_efp[ry1:ry2,rx1:rx2]
            all_region = efp[ry1:ry2,rx1:rx2]

            y1 = cy-cone_rad
            y2 = cy+cone_rad+1
            x1 = cx-cone_rad
            x2 = cx+cone_rad+1
            
            y1 = max(y1,0)
            y2 = min(y2,n_slow-1)
            x1 = max(x1,0)
            x2 = min(x2,n_fast-1)

            this_cone = this_efp[y1:y2,x1:x2]
            all_cone = efp[y1:y2,x1:x2]

            this_profile = np.ones(n_depth)*np.nan

            # we're only going to use this cone in this volume if it has adequate coverage, e.g. 75% of pixels in the region and cone
            if len(np.where(this_cone)[0])>(2*cone_rad)**2*.75 and len(np.where(this_region)[0])>(2*region_rad)**2*.75:

                try:

                    nxc = np.fft.fftshift(np.abs(np.fft.ifft2(np.fft.fft2(this_region)*np.conj(np.fft.fft2(all_region)))))
                    yerr,xerr = np.where(nxc==np.max(nxc))
                    ny,nx = np.shape(nxc)
                    yerr = yerr - ny//2
                    xerr = xerr - nx//2

                    if np.abs(xerr)<=cone_rad and np.abs(yerr)<=cone_rad:
                        newcx = cx+xerr
                        newcy = cy+yerr
                    else:
                        newcy = cy
                        newcx = cx

                    y1 = int(round(newcy-cone_rad))
                    y2 = int(round(newcy+cone_rad+1))
                    x1 = int(round(newcx-cone_rad))
                    x2 = int(round(newcx+cone_rad+1))

                    y1 = max(y1,0)
                    y2 = min(y2,n_slow-1)
                    x1 = max(x1,0)
                    x2 = min(x2,n_fast-1)


                    subvol = vol[y1:y2,:,x1:x2]
                    subt = t[y1:y2,x1:x2]
                    #test = np.mean(np.abs(subvol[:,isos_idx-5:cost_idx+5:,:]),axis=1)
                    #nany,nanx = np.where(np.isnan(test))
                    subvol = np.transpose(subvol,(1,0,2))

                    sz,sy,sx = subvol.shape
                    this_sheet = np.reshape(subvol,(sz,sy*sx))
                    this_profile = np.nanmean(np.abs(this_sheet),1)
                    this_time = np.nanmean(subt)

                    shift = int(utils.findShift(this_profile,model))

                    this_sheet_cropped = this_sheet[isos_idx-40-shift:isos_idx+60-shift,:]
                    brightest_idx = np.sum(np.abs(this_sheet),0)
                    brightest_idx[np.where(np.isnan(brightest_idx))] = -np.inf
                    brightest_idx = np.argmax(brightest_idx)
                    this_profile_cropped = np.nanmean(np.abs(this_sheet_cropped),1)
                    cropped_model = model[isos_idx-40:isos_idx+60]
                    brightest_profile_phase = np.angle(this_sheet_cropped[:,brightest_idx])

                    intensity_sheet[:,cone_idx] = this_profile_cropped
                    phase_sheet[:,cone_idx] = brightest_profile_phase
                    time_sheet[cone_idx] = this_time

                except Exception as e:
                    print e

        series_hive['frames/%s/cones'%key].put('intensity',intensity_sheet)
        series_hive['frames/%s/cones'%key].put('phase',phase_sheet)
        series_hive['frames/%s/cones'%key].put('time',time_sheet)

    intensity_cube[tidx,:,:] = intensity_sheet
    time_cube[tidx,:] = time_sheet
    tidx = tidx + 1


try:
    os.mkdir(os.path.join(tag,'cone_catalog'))
except Exception as e:
    print e

# time interval over which to integrate efps and cone profiles
interval = 0.2
starts = np.arange(0.0,9.5,interval)
ends = starts+interval
interval_keys = np.round((starts/interval)).astype(np.integer)

n_intervals = len(interval_keys)
efp_stack = np.zeros((n_intervals,efpsy,efpsx))
efp_counter_stack = np.ones((n_intervals,efpsy,efpsx))*1e-16

make_interval_efps = False
if make_interval_efps:
    for key in keys[idx1:idx2]:#[89*2:]:
        time_sheet = series_hive['frames/%s/cones'%key].get('time')
        this_efp = series_hive['frames/%s/projections'%key].get('CONES')
        try:
            t = np.nanmean(time_sheet)
            t_idx = int(np.floor(t/interval))
            efp_stack[t_idx,:,:] = efp_stack[t_idx,:,:]+this_efp
            #efp_counter_stack[np.where(this_efp[t_idx,:,:])] = efp_counter_stack[np.where(this_efp[t_idx,:,:])]+1.0
            efp_counter_stack[t_idx,:,:][np.where(this_efp)] = efp_counter_stack[t_idx,:,:][np.where(this_efp)] + 1.0
        except ValueError as ve:
            print ve
            continue

    efp_stack = efp_stack/efp_counter_stack
    efp_stack[np.where(efp_stack==0)] = np.nan
    cone_coordinates = np.zeros((len(cones_x),3))
    #plt.figure(figsize=(efpsx/100.,efpsy/100.))
    #plt.axes([0,0,1,1])
    temp = np.nanmean(efp_stack,0)
    temp = utils.nanreplace(temp,np.nanmean(temp))
    #plt.imshow(temp,cmap='gray')
    #plt.savefig('for_mehdi/en_face_projection.png',dpi=300)
    #plt.savefig(os.path.join(tag,'cone_catalog/en_face_projection.png'),dpi=300)
    #plt.autoscale(False)
    for cone_index,(cx,cy) in enumerate(zip(cones_x,cones_y)):
        cone_coordinates[cone_index,0] = cone_index
        cone_coordinates[cone_index,1] = cx
        cone_coordinates[cone_index,2] = cy
        #plt.text(cx,cy,'%04d'%cone_index,ha='center',va='center',color='g',fontsize=4)
    #plt.savefig('for_mehdi/en_face_projection_marked.png',dpi=300)
    #plt.savefig(os.path.join(tag,'cone_catalog/en_face_projection_marked.png'),dpi=300)

    mdict = {}
    mdict['en_face_projection_stack'] = efp_stack
    mdict['time_interval_starts'] = starts-2.0
    mdict['time_interval_ends'] = ends-2.0
    mdict['cone_coordinates'] = cone_coordinates
    spi.savemat('for_mehdi/en_face_projection_stack.mat',mdict)

    np.save(os.path.join(tag,'cone_catalog/interval_stack.npy'),efp_stack)
    
    #plt.show()
efp_stack = np.load(os.path.join(tag,'cone_catalog/interval_stack.npy'))
np.save(os.path.join(tag,'cone_catalog/cones_x.npy'),cones_x)
np.save(os.path.join(tag,'cone_catalog/cones_y.npy'),cones_y)

valid_cone_indices = []

print intensity_cube.shape
all_cone_mean = np.nanmean(np.nanmean(intensity_cube,2),0)
peaks = utils.find_peaks(all_cone_mean)

#plt.plot(all_cone_mean)
#plt.plot(peaks,all_cone_mean[peaks],'rs')
#plt.show()

try:
    os.mkdir(os.path.join(tag,'cone_catalog/individual_cone_plots_png'))
    os.mkdir(os.path.join(tag,'cone_catalog/individual_cone_plots_pdf'))
except Exception as e:
    print e

plt.figure(figsize=(7,6))
for cidx in range(n_cones):
    cone_sheet = intensity_cube[:,:,cidx]
    time_sheet = time_cube[:,cidx]
    test = np.sum(cone_sheet,axis=1)
    total_valid_count = np.sum(1-np.isnan(test))
    prestim_valid_count = np.sum(1-np.isnan(test[np.where(time_sheet<stimulus_time)]))
    
    #if total_valid_count<300 or prestim_valid_count<30:
    if total_valid_count<150 or prestim_valid_count<5:
        continue

    valid_cone_indices.append(cidx)

    analyze_cone(cone_sheet,time_sheet,cidx,stimulus_time,cones_x[cidx],cones_y[cidx],efp_stack,starts,ends)
    plt.pause(.000001)
    plt.savefig(os.path.join(tag,'cone_catalog/individual_cone_plots_png/%04d.png'%cidx),dpi=300)
    plt.savefig(os.path.join(tag,'cone_catalog/individual_cone_plots_pdf/%04d.pdf'%cidx))
    continue

plt.figure(figsize=(efpsx/100.,efpsy/100.))
plt.axes([0,0,1,1])
plt.imshow(efp,cmap='gray')
plt.autoscale(False)
for valid_index in valid_cone_indices:
    cx = cones_x[valid_index]
    cy = cones_y[valid_index]
    plt.text(cx,cy,'%04d'%valid_index,ha='center',va='center',color='r',fontsize=4)

plt.xticks([])
plt.yticks([])

plt.savefig(os.path.join(tag,'cone_catalog/cone_locations.png'),dpi=300)
plt.show()

