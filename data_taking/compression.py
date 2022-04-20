import numpy as np
import matplotlib.pyplot as pl
import gzip
import h5py
import time
import os



start_t = time.time()

save_mode = "npy" #"npy" # options are "npy", "h5py", "none"

#data_dir = "C:/Users/ryanm/Documents/Research/Work zone/20220111_test/20220111_test/"
#data_dir = "G:/.shortcut-targets-by-id/11qeqHWCbcKfFYFQgvytKem8rulQCTpj8/crystalize/data/data-202110/20211012/20211012_1656_Po_Co_OCVtop_0.0g_0.0c_1.19bar_3mv_25us_circ_20min/"
#data_dir = "G:/.shortcut-targets-by-id/11qeqHWCbcKfFYFQgvytKem8rulQCTpj8/crystalize/data/data-202201/20220131/testForAlign/"

with open("/home/xaber/Analysis/solid_xenon_tpc/pulse_finding_scripts/path.txt", 'r') as path:
    data_dir = path.read()
    data_dir = data_dir[:-1]

#data_dir = "/home/xaber/Data/data-202204/20220419/202204191736_1.2bar_3500C_3200G_1000A_54B_15us_0.5Vpp_2coin_3mVtrig/"

save_dir = data_dir+"compressed_data/"

try:
    os.mkdir(save_dir)
except:
    print("Directory already exists")

n_boards = 3
n_sipms = [16,8,8]
n_all_ch = int(np.sum(n_sipms))

wsize = 7500+8 #7500+8 #12500+8 #12500+8 #12500+8 #12500+8 #3000+8 # 8 = size of header 
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
    all_data_front = np.zeros((n_all_ch, max_n_events, wsize-8))
    all_data_back = np.zeros((n_all_ch, max_n_events, wsize-8))
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
                        all_data_front[ch_ind, ev-toss_counts, :] = ch_data[ev, 8:wsize] - np.mean(ch_data[ev, 8:8+150])
                        all_data_back[ch_ind, ev-toss_counts, :] = ch_data[ev, 8:wsize] - np.mean(ch_data[ev, -150:-1])
                        test_ev[bd,ev-toss_counts] = ch_data[ev,2]
                    elif bd == 1:
                        all_data_front[ch_ind, ev-toss_counts, delay:] = ch_data[ev, 8:(wsize-delay)] - np.mean(ch_data[ev, 8:8+150])
                        all_data_back[ch_ind, ev-toss_counts, delay:] = ch_data[ev, 8:(wsize-delay)] - np.mean(ch_data[ev, -150:-1])
                        test_ev[bd,ev-toss_counts] = ch_data[ev,2]
                    elif bd == 2:
                        all_data_front[ch_ind, ev-toss_counts, 2*delay:] = ch_data[ev, 8:(wsize-2*delay)] - np.mean(ch_data[ev, 8:8+150])
                        all_data_back[ch_ind, ev-toss_counts, 2*delay:] = ch_data[ev, 8:(wsize-2*delay)] - np.mean(ch_data[ev, -150:-1])
                        test_ev[bd,ev-toss_counts] = ch_data[ev,2]
                    # need to scale by spe size!!!
            
            ch_ind += 1






    #for ch in range(n_all_ch):
        #print(np.mean(all_data_front[ch,0,8+2*delay:100+2*delay]))
        #all_data_front[ch,:,:] -= np.mean(all_data_front[ch,0,200:300])
        #if ch < 4:
        #    all_data_front[ch,:,:] = np.zeros_like(all_data_front[0,:,:])
        #    all_data_back[ch,:,:] = np.zeros_like(all_data_back[0,:,:])


    # Reshape into pods
    pod_size = 10 #5 # samples. 5 samples = 10 ns
    n_pods = int( (wsize-8)/pod_size )
    all_data_front_pods = np.reshape(all_data_front, (n_all_ch, max_n_events, n_pods, pod_size) )
    all_data_back_pods = np.reshape(all_data_back, (n_all_ch, max_n_events, n_pods, pod_size) )

    # Sum channels
    sum_data_front_pods = np.sum(all_data_front_pods, axis=0)
    sum_data_back_pods = np.sum(all_data_back_pods, axis=0)

    # Sum within pods
    sum_data_front_pods_area = np.sum(sum_data_front_pods, axis=2) 
    sum_data_back_pods_area = np.sum(sum_data_back_pods, axis=2) 

    # Do cuts on areas
    area_threshold_front = 1000 #200 # one day this will be phe 
    area_threshold_back = 1000 #475 #area_threshold_front
    toSaveOrNotToSave = (np.abs(sum_data_front_pods_area) > area_threshold_front)*(np.abs(sum_data_back_pods_area) > area_threshold_back)
   
    nBefore = 25
    for i in range(nBefore):
        diffSave = np.diff(toSaveOrNotToSave, axis=1)
        #print(np.sum(diffSave[diffSave==1]))
        if i < 8: toSaveOrNotToSave[:,:-1][diffSave != 0] = 1
        toSaveOrNotToSave[:,1:][diffSave != 0] = 1
   

    # Create array to save
    stuffToSave = np.zeros_like(all_data_front_pods)
    for ch in range(n_all_ch):
        stuffToSave[ch, toSaveOrNotToSave, :] = all_data_front_pods[ch, toSaveOrNotToSave, :]

    stuffToSave = np.reshape(stuffToSave, (n_all_ch, max_n_events*(wsize-8)))
    
    """
    for i in range(max_n_events):
        #print(test_ev[:,i])
        pl.figure()
        #for ch in range(24,32):
        #    pl.plot(all_data_front[ch, i, :],"g")
        pl.plot(np.sum(all_data_front[:, i, :],axis=0),"b")
        pl.plot(np.sum(stuffToSave[:, i*(wsize-8):(i+1)*(wsize-8)], axis=0),"r")
        #pl.legend(["Raw","Compressed"])
        #pl.plot(np.sum(stuffToSave[16:22, i*(wsize-8):(i+1)*(wsize-8)], axis=0),"b")
        #pl.plot(np.sum(stuffToSave[23:31, i*(wsize-8):(i+1)*(wsize-8)], axis=0),"m")
        
        pl.show()
    """
    
    
    
    print("Percentage saved: "+str(np.count_nonzero(stuffToSave)/stuffToSave.size) )


    # Save that mf
    if save_mode == "npy":
        with open(save_dir+"compressed_"+str(bk)+".npy", "wb") as f:
            np.savez_compressed(f, stuffToSave.flatten() )
    elif save_mode == "h5py":
        with h5py.File(save_dir+"compressed_"+str(bk)+".h5", "w") as f:
            f.create_dataset("dataset", data=stuffToSave, compression="gzip" )
    elif save_mode == "none":
        continue

    


end_t = time.time()
print("Finished zero baseline reduction", end_t-start_t, "sec")



