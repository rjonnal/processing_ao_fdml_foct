import numpy as np
import scipy.io as sio
from matplotlib import pyplot as plt
import sys,os
from cuttlefish import make_configuration, Hive
from cuttlefish import utils #.utils import shear,find_peaks,ascend,descend
import scipy.optimize as soo
from fig2gif import GIF
import glob

# This program converts the series of B-scan MAT files (ideally complex data) into
# numpy arrays, each of which contains 1 volume of data

# use sys.argv for running as a shell script, but use explicit definition of 
# dataset tag for running in Spyder/IPython

# important constants:
volume_left_x_start = 0
volume_left_x_end = 160
redo = False # delete everything and start over

try:
    tag = sys.argv[1].strip('/').strip('\\')
except Exception as e:
    print e
    tag = '2018.05.09_09.59.25_RE_3DT_10uw-320x-PROC-5Reps-0.3_PVthresh'

cone_projection_directory = os.path.join(tag,'cone_projections')

try:
    os.mkdir(cone_projection_directory)
except Exception as e:
    print e

hive = Hive(tag)

# coordinates for forward and backward linear regions
# DataSet1 start 30 end 100


#flist = ['./%s/AMP_L_%05d.mat'%(tag,k) for k in range(1,48000)]
flist = glob.glob('%s/mat/*.mat'%tag)
flist.sort()
n_scans = len(flist)

mat_keys = sio.loadmat(flist[0]).keys()
mat_keys.remove('__version__')
mat_keys.remove('__header__')
mat_keys.remove('__globals__')
data_key = mat_keys[0]

initial_image = sio.loadmat(flist[0])[data_key]
n_depth = initial_image.shape[0]
hive.put('config/n_depth',n_depth)

try:
    os.mkdir('./tmp')
except:
    pass

def get(n):
    return sio.loadmat(flist[n])[data_key]

# Create an en-face projection of the entire volume series
# This will be used initially to determine the locations of y-axis
# turnarounds.
# Try to load it from a file if possible. If file doesn't exist, make the projection
# and save it so we don't have to do it again.

# Projection range selected by empirically examining profile.
# It will be used to specify the range of projection relative
# to the center of mass (com).
projection_range = np.arange(-30,40).astype(np.int16) # relative to center of mass
en_face_projection_fn = os.path.join('./tmp','%s_enface_projection.npy'%tag)

scan_fraction = .1 # fraction of scans used to generate bulk en face projection

try:
    assert redo==False
    en_face_projection = np.load(en_face_projection_fn)
except Exception as e:
    en_face_projection = []
    for k in range(int(round(n_scans*scan_fraction))):
        try:
            frame = np.abs(get(k))
            prof = frame.mean(1)
            c = utils.com(prof)
            p = np.mean(frame[c+projection_range,:],0)
            en_face_projection.append(p)
            print k
        except Exception as e:
            print e
    en_face_projection = np.array(en_face_projection)
    np.save(en_face_projection_fn,en_face_projection)


scan_width = en_face_projection.shape[1]
    
#plt.imshow(en_face_projection[:1000,:],aspect='auto')
#plt.show()
#sys.exit()
# Try to identify the x-scanner turnaround by walking a vertical
# strip across the image and correlating it with its mirror image
test_width = 5
x_start = en_face_projection.shape[1]//2-test_width
x_end = en_face_projection.shape[1]//2+test_width
test_range = range(x_start,x_end)
hcorrs = []
for k in test_range:
    left = en_face_projection[:,k-test_width:k].ravel()
    right = en_face_projection[:,k+test_width:k:-1].ravel()
    hcorrs.append(utils.pearson(left,right))

    
x_turnaround = test_range[np.argmax(hcorrs)]

volume_right_x_start = x_turnaround+(x_turnaround-volume_left_x_start)
volume_right_x_end = x_turnaround+(x_turnaround-volume_left_x_end)

x_max = en_face_projection.shape[1]

# We need to make sure the start of the right side (flyback) volume
# is within the measured data; this may not happen if the range is
# set to 0,160 for the left side volume but the resonant scanner phase
# is less than 0 at that pixel. The result will be that the desired
# mirror image start will be beyond the edge of the image.
while volume_right_x_start>=x_max:
    volume_right_x_start-=1
    volume_left_x_start+=1

n_fast = volume_left_x_end-volume_left_x_start
hive.put('config/volume_left_x_start',volume_left_x_start)
hive.put('config/volume_left_x_end',volume_left_x_end)
hive.put('config/n_fast',n_fast)
hive.put('config/volume_right_x_start',volume_right_x_start)
hive.put('config/volume_right_x_end',volume_right_x_end)


    
# Now we try to identify the y-galvo turnarounds
# by walking down the projection and correlating adjacent
# strips with one of them inverted.
y_turnaround_corrs_fn = os.path.join('./tmp','%s_y_galvo_turnaround_corrs.npy'%tag)
try:
    corrs = np.load(y_turnaround_corrs_fn)
except IOError as e:
    corrs = []
    for k in range(en_face_projection.shape[0]-test_width*2):
        try:
            y1 = k
            y2 = k+test_width
            y3 = k+test_width*2
            # cut the strips out, flip one of them, and
            # arrange the strips as two column vectors for
            # quick correlation computation
            im1 = en_face_projection[y1:y2,:].ravel()
            im2 = en_face_projection[y3:y2:-1,:].ravel()

            
            corrs.append(utils.pearson(im1,im2))
            if False: # make true to visualize
                plt.cla()
                plt.plot(corrs)
                plt.axhline(0.5)
                plt.pause(.000001)
        except Exception as e:
            corrs.append(-1)
            print e

    corrs = np.array(corrs)
    np.save(y_turnaround_corrs_fn,corrs)

N = len(corrs)*4
#plt.plot(corrs)
#plt.show()

spectrum = np.fft.fftshift(np.fft.fft(corrs,n=N))[N//2:]
freq = np.fft.fftshift(np.fft.fftfreq(N))[N//2:]

# zero a safe region, < 500 px, to avoid DC
spectrum[np.where(freq<=1.0/500.0)] = 0
peak_locations = utils.find_peaks(np.abs(spectrum))
peak_heights = np.abs(spectrum[peak_locations])
peak_index = peak_locations[np.argmax(peak_heights)]
fundamental = freq[peak_index]
n_slow = int(round(1.0/fundamental))

print n_slow

# Now we need to identify the first turnaround to start the first volume.
corr_means = []
# Go through the first n_slow rows and compute the mean turnaround correlations
# for each row and its n_slow-separated counterparts through the series
for offset in range(n_slow):
    comb = range(offset,len(corrs),n_slow)
    corr_means.append(np.sum(corrs[comb])/float(len(comb)))
    
y_start = np.argmax(corr_means)+test_width
    
hive.put('config/n_slow',n_slow)
hive.put('config/y_start',y_start)
      
# Now let's make a time matrix indicating the time of acquisition for each
# pixel in the volume series. Cutouts from this matrix will be stored with volumes
# for easy timing.
time_table = np.reshape(np.arange(n_scans*scan_width)*1.0/1.6e6,(n_scans,scan_width))

hive.put('config/n_scans',n_scans)
      
# this series is too big to hold in memory, so let's treat
# each book-matched pair as a separate volume series, and put
# it in its own hive
do_plot = False

if do_plot:
    fig = plt.figure(figsize=(12,6))
for idx,y1 in enumerate(range(y_start,n_scans-n_slow,n_slow)):
    root_name = os.path.join(os.path.join('.',tag),'volumes/%04d'%idx)
    hive = Hive(root_name)
    y2 = y1+n_slow

    # we need to reverse the arrangement of
    # alternate scans because the scanner was
    # in triangle not sawtooth mode
    odd_scan = idx%2==1
    
    if not hive.has('data_time') or redo:
        data_time = np.zeros((2,n_slow,n_fast))
        for zidx,frame_index in enumerate(range(y1,y2)):
            # we need to reverse the arrangement of
            # alternate scans because the scanner was
            # in triangle not sawtooth mode
            if odd_scan:
                insertion_index = n_slow - zidx - 1
            else:
                insertion_index = zidx
            
            # frame_index references the whole series of B-scans
            # insertion_index references the location in this volume series
            # in which to insert the data
            f = time_table[frame_index,:]
            left = f[volume_left_x_start:volume_left_x_end]
            right = f[volume_right_x_start:volume_right_x_end:-1]
            data_time[0,insertion_index,:] = left
            data_time[1,insertion_index,:] = right
        hive.put('data_time',data_time)
    else:
        data_time = hive.get('data_time')    
        
    if not hive.has('processed_data') or redo:
        processed_data = np.zeros((2,n_slow,n_depth,n_fast),dtype=np.complex64)
        for zidx,frame_index in enumerate(range(y1,y2)):
            if odd_scan:
                insertion_index = n_slow - zidx - 1
            else:
                insertion_index = zidx
            
            # frame_index references the whole series of B-scans
            # insertion_index references the location in this volume series
            # in which to insert the data
            f = get(frame_index)
            left = f[:,volume_left_x_start:volume_left_x_end]
            right = f[:,volume_right_x_start:volume_right_x_end:-1]
            processed_data[0,insertion_index,:,:] = left
            processed_data[1,insertion_index,:,:] = right
        hive.put('processed_data',processed_data)
    else:
        processed_data = hive.get('processed_data')

    if do_plot:
        # log this process
        log_dir = '%s_projection_log'%tag
        try:
            os.makedirs(log_dir)
        except Exception as e:
            print e
            pass
        
    if not (hive.has('projections/ISOS') and hive.has('projections/COST')) or redo:
        # we have to make double data structures for each
        # model, projection, and label set, since we're dealing with
        # two separate volumes
        profiles = []
        label_dictionary = {}
        projection_dictionary = {}
        layers = ['ELM','ISOS','COST','RPE','SISOS']
        for l in layers:
            label_dictionary[l] = []
            projection_dictionary[l] = np.zeros((2,n_slow,n_fast))

        flat_offsets_pair = np.zeros((2,n_slow,n_fast))

        model_directory = os.path.join(root_name,'model')
        try:
            os.mkdir(model_directory)
        except Exception as e:
            print e
        
        for vidx in range(2):
            avol = np.abs(processed_data[vidx,:,:,:])
            sy,sz,sx = avol.shape
            flat_offsets = utils.get_flat_offsets(avol)
            flat_offsets_pair[vidx,:,:] = flat_offsets
            flattened_avol = np.zeros((sy,sz+flat_offsets.max(),sx))
            for y in range(sy):
                for x in range(sx):
                    flattened_avol[y,flat_offsets[y,x]:flat_offsets[y,x]+sz,x] = avol[y,:,x]

            # we should just make the model and try to label it here--which allows
            # us to skip ao_fdml_foct_step_1_make...
            model_profile = flattened_avol.mean(2).mean(0)
            profiles.append(model_profile)
            peaks = utils.find_peaks(model_profile)
            peak_heights = model_profile[peaks]
            # Let's try to assign is/os, cost, and rpe to the first three peaks
            # that are above 1/3 of the profile maximum
            try:
                isos,cost,rpe = peaks[np.where(peak_heights>model_profile.max()/3.0)]
            except Exception as e:
                sys.exit(e)

            # now let's guess the SISOS location--between ISOS and COST but
            # slightly closer to the ISOS:
            sisos = (isos+cost)//2-2
            # now let's try to locate the ELM based on approximate location:
            elm_min = isos-30
            elm_max = isos-10
            try:
                inner_peaks = peaks[np.where(np.logical_and(peaks<elm_max,peaks>elm_min))]
                inner_peak_heights = model_profile[inner_peaks]
                elm = inner_peaks[np.argmax(inner_peak_heights)]
            except Exception as e:
                print e, 'attempt 1 failed'
                inner_peaks = peaks[np.where(peaks<elm_max)]
                inner_peak_heights = model_profile[inner_peaks]
                try:
                    elm = inner_peaks[np.argmax(inner_peak_heights)]
                except Exception as e:
                    print e, 'attempt 2 failed; elm = -1'
                    elm = -1

            label_dictionary['ELM'].append(elm)
            label_dictionary['ISOS'].append(isos)
            label_dictionary['SISOS'].append(sisos)
            label_dictionary['COST'].append(cost)
            label_dictionary['RPE'].append(rpe)
            
            plt.plot(model_profile)
            # how many pixels (in depth) should the projections include:
            projection_half_depth = 3
            for k in label_dictionary.keys():
                depth = label_dictionary[k][-1]
                height = model_profile[depth]
                plt.text(depth,height,k,ha='center',va='bottom')
                if k in ['ELM','ISOS','SISOS','COST']:
                    projection = flattened_avol[:,depth-projection_half_depth:depth+projection_half_depth+1,:].mean(1)
                    projection_dictionary[k][vidx,:,:] = projection
                # need to handle RPE separately because projecting from both sides
                # includes too much COST
                elif k in ['RPE']:
                    projection = flattened_avol[:,depth:depth+projection_half_depth+1,:].mean(1)
                    projection_dictionary[k][vidx,:,:] = projection
                    
            plot_filename = os.path.join(model_directory,'volume_%03d_model.png'%vidx)
            plt.savefig(plot_filename)
            plt.close()

        for k in projection_dictionary.keys():
            proj = projection_dictionary[k]
            hive.put('projections/%s'%k,proj)
        
        for lkey in label_dictionary.keys():
            hive.put('model/labels/%s'%lkey,np.array(label_dictionary[lkey]))
        plengths = []
        pcount = len(profiles)
        max_length = -np.inf
        for p in profiles:
            if len(p)>max_length:
                max_length = len(p)

        profile_array = np.zeros((pcount,max_length))
        for pidx,p in enumerate(profiles):
            profile_array[pidx,:len(p)] = p
        hive.put('model/profiles',profiles)

        hive.put('flat_offsets',flat_offsets_pair)

        cones = (projection_dictionary['ISOS']+projection_dictionary['COST'])/2.0
        hive.put('projections/CONES',cones)
        
        for vidx in range(2):
            cone_proj = cones[vidx,:,:]
            plt.cla()
            plt.imshow(cone_proj,cmap='gray',clim=np.percentile(cone_proj,(1,99.9)))
            plt.title('volume %04d %02d'%(idx,vidx))
            plt.savefig(os.path.join(cone_projection_directory,'%04d_%02d.png'%(idx,vidx)))
            plt.close()
