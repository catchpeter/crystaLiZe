import numpy as np
import matplotlib.pyplot as pl
import gzip
import h5py
import shutil
import time

start_t = time.time()

save_mode = "npy" # options are "npy" or "h5py"

#data_dir = "C:/Users/ryanm/Documents/Research/Work zone/20220111_test/20220111_test/"
#data_dir = "G:/.shortcut-targets-by-id/11qeqHWCbcKfFYFQgvytKem8rulQCTpj8/crystalize/data/data-202110/20211012/20211012_1656_Po_Co_OCVtop_0.0g_0.0c_1.19bar_3mv_25us_circ_20min/"
data_dir = "G:/.shortcut-targets-by-id/11qeqHWCbcKfFYFQgvytKem8rulQCTpj8/crystalize/data/data-202201/20220131/testForAlign/"
save_dir = ""

n_boards = 3
n_sipms = [16,8,8]

wsize = 3000+8 # 8 = size of header 

load_dtype = "int16"


# First, check board alignment 
# Get list of event numbers in header
evNum = np.array([])
n_events = []
for bd in range(n_boards):
    ch0_data = np.fromfile(data_dir+"waveforms_"+str(bd)+"_0.dat", dtype="int16")
    n_events.append(int(ch0_data.size/wsize))
    evNum = np.concatenate((evNum, ch0_data[2::wsize]))



# Check for misaligned events
tosser = []
lastToCheck = int(np.max(evNum) )
for i in range(1,lastToCheck):
    if np.count_nonzero(evNum == i) != 3: tosser.append(i)



# Do compression and remove any misaligned events

# Loop over boards
for bd in range(n_boards):
    print("Board",bd)

    # Loop over channels
    for si in range(n_sipms[bd]):

        # Load data
        ch_data = np.fromfile(data_dir + "waveforms_"+str(bd)+"_"+str(si)+".dat", dtype=load_dtype)
        ch_data = np.reshape(ch_data, (n_events[bd], wsize))
        only_header = ch_data[:,0:7]
        only_data = ch_data[:,8:wsize]
        stuffToSave = np.zeros_like(only_data)

        # Loop over events
        toss_counts = 0
        for ev in range(n_events[bd]):
            baseline = np.mean(only_data[ev,0:50])
            sigma = np.std(only_data[ev,0:50])

            if only_header[ev,2] in tosser:
                toss_counts += 1
                continue
            else:
                toSaveOrNotToSave = only_data[ev,:] > baseline + 2*sigma
                stuffToSave[ev-toss_counts, toSaveOrNotToSave] = only_data[ev, toSaveOrNotToSave] - baseline

        # Save that mf
        if save_mode == "npy":
            with open(save_dir+"compressed_"+str(bd)+"_"+str(si)+".npy", "wb") as f:
                np.savez_compressed(f, stuffToSave.flatten() )
        elif save_mode == "h5py":
            with h5py.File(save_dir+"compressed_"+str(bd)+"_"+str(si)+".h5", "w") as f:
                f.create_dataset("dataset", data=stuffToSave.flatten(), compression="gzip" )
        else:
            with open(save_dir+"compressed_"+str(bd)+"_"+str(si)+".npy", "wb") as f:
                np.savez_compressed(f, stuffToSave.flatten() )

                   
end_t = time.time()


print("Finished zero baseline reduction", end_t-start_t, "sec")


#test = np.load(save_dir+"compressed_0_0.npy")
#print(np.sum(test["arr_0"])/3000)
