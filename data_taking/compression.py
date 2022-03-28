import numpy as np
import matplotlib.pyplot as pl
import gzip
import h5py
import time
import os



start_t = time.time()

save_mode = "npy" # options are "npy" or "h5py"

#data_dir = "C:/Users/ryanm/Documents/Research/Work zone/20220111_test/20220111_test/"
#data_dir = "G:/.shortcut-targets-by-id/11qeqHWCbcKfFYFQgvytKem8rulQCTpj8/crystalize/data/data-202110/20211012/20211012_1656_Po_Co_OCVtop_0.0g_0.0c_1.19bar_3mv_25us_circ_20min/"
#data_dir = "G:/.shortcut-targets-by-id/11qeqHWCbcKfFYFQgvytKem8rulQCTpj8/crystalize/data/data-202201/20220131/testForAlign/"
data_dir = "/home/xaber/Data/data-202203/20220324/202203241628_1.4bar_2600C2400G500A_54B/"
save_dir = data_dir+"compressed_data/"

os.mkdir(save_dir)

n_boards = 3
n_sipms = [16,8,8]
n_all_ch = int(np.sum(n_sipms))

wsize = 12500+8 #12500+8 #12500+8 #12500+8 #3000+8 # 8 = size of header 
block_size = 1500 # This will also be number of events saved per file
n_blocks = 1000

delay = 24 #  48 #48

load_dtype = "int16"


ch0_data = np.fromfile(data_dir+"waveforms_"+str(0)+"_0.dat",dtype=load_dtype)
tot_ev = int(ch0_data.size/wsize)
tot_fi = int(np.ceil(tot_ev/block_size))
print("Total events: "+str(tot_ev) )
print("Number of compressed files = "+str(tot_fi))
time.sleep(5)


# Loop over blocks of events
for bk in range(tot_fi+1):

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
    firstToCheck = int(np.min(evNum) ) 
    for i in range(firstToCheck,lastToCheck+1):
        if np.count_nonzero(evNum == i) != 3 and np.count_nonzero(evNum == i) > 0: tosser.append(i)

    print("Number of misaligned events: "+str(len(tosser)))

    # Load data, loop over all boards and channels
    all_data = np.zeros((n_all_ch, max_n_events, wsize-8))
    test_ev = np.zeros((n_boards, max_n_events))
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
                    if bd == 0:
                        all_data[ch_ind, ev-toss_counts, :] = ch_data[ev, 8:wsize] 
                        test_ev[bd,ev-toss_counts] = ch_data[ev,2]
                    elif bd == 1:
                        all_data[ch_ind, ev-toss_counts, delay:] = ch_data[ev, 8:wsize-delay]
                        test_ev[bd,ev-toss_counts] = ch_data[ev,2]
                    elif bd == 2:
                        all_data[ch_ind, ev-toss_counts, 2*delay:] = ch_data[ev, 8:wsize-2*delay]
                        test_ev[bd,ev-toss_counts] = ch_data[ev,2]
                    # need to scale by spe size!!!
            
            ch_ind += 1






    for ch in range(n_all_ch):
        #print(np.mean(all_data[ch,0,8+2*delay:100+2*delay]))
        all_data[ch,:,:] -= np.mean(all_data[ch,0,8+2*delay:100+2*delay])
        if ch < 4:
            all_data[ch,:,:] = np.zeros_like(all_data[0,:,:])


    # Reshape into pods
    pod_size = 5 # samples. 5 samples = 10 ns
    n_pods = int( (wsize-8)/pod_size )
    all_data_pods = np.reshape(all_data, (n_all_ch, max_n_events, n_pods, pod_size) )


    # Take sum, do baseline subtraction, do cut on area of pods
    sum_data_pods = np.sum(all_data_pods, axis=0)
    #baselines = np.mean(sum_data_pods[:,0:10,:], axis=(1,2))
    sum_data_pods_area = np.sum(sum_data_pods, axis=2)  #np.sum( np.subtract(sum_data_pods, baselines[:,None,None]) , axis=2)
    area_threshold = 50 # one day this will be phe 
    toSaveOrNotToSave = sum_data_pods_area > area_threshold 


    # Create array to save
    stuffToSave = np.zeros_like(all_data_pods)
    for ch in range(n_all_ch):
        stuffToSave[ch, toSaveOrNotToSave, :] = all_data_pods[ch, toSaveOrNotToSave, :]

    stuffToSave = np.reshape(stuffToSave, (n_all_ch, max_n_events*(wsize-8)))
    
    #for i in range(max_n_events):
        #print(test_ev[:,i])
        #pl.figure()
        #pl.plot(np.sum(stuffToSave[0:15, i*(wsize-8):(i+1)*(wsize-8)], axis=0),"r")
        #pl.plot(np.sum(stuffToSave[16:22, i*(wsize-8):(i+1)*(wsize-8)], axis=0),"b")
        #pl.plot(np.sum(stuffToSave[23:31, i*(wsize-8):(i+1)*(wsize-8)], axis=0),"m")
        
        #pl.show()
    
    
    
    print("Percentage saved: "+str(np.count_nonzero(stuffToSave)/stuffToSave.size) )


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



