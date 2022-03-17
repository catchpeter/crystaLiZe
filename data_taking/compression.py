import numpy as np
import matplotlib.pyplot as pl
import gzip
import h5py
import time

start_t = time.time()

save_mode = "npy" # options are "npy" or "h5py"

#data_dir = "C:/Users/ryanm/Documents/Research/Work zone/20220111_test/20220111_test/"
#data_dir = "G:/.shortcut-targets-by-id/11qeqHWCbcKfFYFQgvytKem8rulQCTpj8/crystalize/data/data-202110/20211012/20211012_1656_Po_Co_OCVtop_0.0g_0.0c_1.19bar_3mv_25us_circ_20min/"
data_dir = "G:/.shortcut-targets-by-id/11qeqHWCbcKfFYFQgvytKem8rulQCTpj8/crystalize/data/data-202201/20220131/testForAlign/"
save_dir = ""

n_boards = 3
n_sipms = [16,8,8]
n_all_ch = int(np.sum(n_sipms))

wsize = 3000+8 # 8 = size of header 
block_size = 5000 # This will also be number of events saved per file
n_blocks = 1

load_dtype = "int16"

# Loop over blocks of events
for bk in range(n_blocks):

    # First, check board alignment 
    # Get list of event numbers in header
    evNum = np.array([])
    n_events = []
    for bd in range(n_boards):
        ch0_data = np.fromfile(data_dir+"waveforms_"+str(bd)+"_0.dat", dtype=load_dtype, offset=block_size*wsize*bk, count=wsize*block_size)
        n_events.append(int(ch0_data.size/wsize))
        evNum = np.concatenate((evNum, ch0_data[2::wsize]))

    max_n_events = int(np.max(n_events))


    # Check for misaligned events
    tosser = []
    lastToCheck = int(np.max(evNum) )
    for i in range(1,lastToCheck):
        if np.count_nonzero(evNum == i) != 3 && np.count_nonzero(evNum == i) > 0: tosser.append(i)


    # Load data, loop over all boards and channels
    all_data = np.zeros((n_all_ch, max_n_events, wsize-8))
    ch_ind = 0
    for bd in range(n_boards):
        for ch in range(n_sipms[bd]):
            ch_data = np.fromfile(data_dir + "waveforms_"+str(bd)+"_"+str(ch)+".dat", dtype=load_dtype, offset=block_size*wsize*bk, count=wsize*block_size)
            ch_data = np.reshape(ch_data, (n_events[bd], wsize))
            
            # Loop over events to check alignment, then to place in all data array
            toss_counts = 0
            for ev in range(n_events[bd]):
                if ch_data[ev, 2] in tosser:
                    toss_counts += 1
                else:
                    all_data[ch_ind, ev-toss_counts, :] = ch_data[ev, 8:wsize] 
                    # need to scale by spe size!!!
            
            ch_ind += 1


    # Reshape into pods
    pod_size = 5 # samples. 5 samples = 10 ns
    n_pods = int( (wsize-8)/pod_size )
    all_data_pods = np.reshape(all_data, (n_all_ch, max_n_events, n_pods, pod_size) )


    # Take sum, do baseline subtraction, do cut on area of pods
    sum_data_pods = np.sum(all_data_pods, axis=0)
    baselines = np.mean(sum_data_pods[:,0:10,:], axis=(1,2))
    sum_data_pods_area = np.sum( np.subtract(sum_data_pods, baselines[:,None,None]) , axis=2)
    area_threshold = 100 # one day this will be phe 
    toSaveOrNotToSave = sum_data_pods_area > area_threshold 


    # Create array to save
    stuffToSave = np.zeros_like(all_data_pods)
    for ch in range(n_all_ch):
        stuffToSave[ch, toSaveOrNotToSave, :] = all_data_pods[ch, toSaveOrNotToSave, :]

    stuffToSave = np.reshape(stuffToSave, (n_all_ch, max_n_events*(wsize-8)))
    print("Percentage suppressed: "+str(np.count_nonzero(stuffToSave)/stuffToSave.size) )


    # Save that mf
    if save_mode == "npy":
        with open(save_dir+"compressed_"+str(bk)+".npy", "wb") as f:
            np.savez_compressed(f, stuffToSave.flatten() )
    elif save_mode == "h5py":
        with h5py.File(save_dir+"compressed_"+str(bk)+".h5", "w") as f:
            f.create_dataset("dataset", data=stuffToSave, compression="gzip" )
    else:
        with open(save_dir+"compressed_"+str(bk)+".npy", "wb") as f:
            np.savez_compressed(f, stuffToSave )

    


end_t = time.time()
print("Finished zero baseline reduction", end_t-start_t, "sec")



