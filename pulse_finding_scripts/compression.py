import numpy as np
import time
import os

from read_settings import get_event_window

def compression(
    data_dir,
    threshold=300, 
    save_mode="npy", 
    save_everything=False, 
    ret_block='all', 
    tellblocks=False):
    """Checks board alignment, does zero-baseline suppression, compresses data
    """
    
    if data_dir[-1] == "\n": data_dir = data_dir[:-1]
    if data_dir[-1] != "/": data_dir += "/" # in case you forgot...
    if save_everything: threshold = -9999999 # if you want to save everything
    
    if save_mode not in ('none','None',None):
        # Make save directory
        save_dir = data_dir+"compressed_data/"
        try:
            os.mkdir(save_dir)
        except:
            print("Directory already exists")
    
    # Channel variables
    n_boards = 3
    n_sipms = [16,8,8]
    n_all_ch = int(np.sum(n_sipms))
    
    # Get event window 
    event_window = get_event_window(data_dir)
    if event_window < 0: 
        raise ValueError("Invalide event window")
    
    wsize = int(500 * event_window) + 8 # Calculate size of waveform + header
    block_size = int(1500*15/event_window) # Number of events loaded, then saved per compressed file
    delay = 24 # Hardcoded delay between boards. DO NOT CHANGE
    load_dtype = "int16" # int16 saves on memory, fewer compressed files
    
    # Get how many total events there are
    ch0_data = np.fromfile(data_dir+"waveforms_"+str(0)+"_0.dat",dtype=load_dtype)
    tot_ev = int(ch0_data.size/wsize)
    tot_fi = int(np.ceil(tot_ev/block_size))
    if tellblocks:
        return tot_fi
    else:
        print("Total events: "+str(tot_ev) )
        print("Number of compressed files = "+str(tot_fi))
        time.sleep(1)
    
    headers = np.zeros((tot_ev+10,8),dtype=int)
    
    if ret_block in ("all", "All", "ALL"):
        tot_fi_to_loop = range(tot_fi)
    else:
        tot_fi_to_loop = [ret_block]
        print(f"Giving block {ret_block} of {tot_fi}")
    # Loop over blocks of events
    #tot_fi = 2 # custom number of files
    for bk in tot_fi_to_loop:
        # First, check board alignment 
        # Get list of event numbers in header
        evNum = np.array([])
        n_events = []
        for bd in range(n_boards):
            # Offset is in bytes!!!!!!!! Each sample is 2 bytes (int16)
            ch0_data = np.fromfile(
                data_dir+"waveforms_"+str(bd)+"_0.dat", 
                dtype=load_dtype, 
                offset=block_size*wsize*bk*2, 
                count=wsize*block_size)
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
        all_data_front = np.empty((n_all_ch, max_n_events, wsize-8))
        all_data_back = np.empty((n_all_ch, max_n_events, wsize-8))
        test_ev = np.empty((n_boards, max_n_events))
        ch_ind = 0
        for bd in range(n_boards):
            for ch in range(n_sipms[bd]):
                ch_data = np.fromfile(
                    data_dir + "waveforms_"+str(bd)+"_"+str(ch)+".dat", 
                    dtype=load_dtype, 
                    offset=block_size*wsize*bk*2, 
                    count=wsize*block_size)
                try:
                    ch_data = np.reshape(ch_data, (n_events[bd], wsize))
                except:
                    print("Error in reshaping array. If final block, ignore this.")
                    continue
                
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
                            if ch == 0:
                                headers[ev-toss_counts+bk*block_size,:] = ch_data[ev,0:8]
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
        area_threshold_front = threshold 
        area_threshold_back = threshold 
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
        
        
        
        # Save that mf
        if save_mode == "npy":
            np.savez_compressed(f'{save_dir}compressed_{bk}.npy', stuffToSave.flatten())
            np.savez_compressed(f'{save_dir}headers.npy', headers.flatten())
        elif save_mode in ("none","None", None):
            return stuffToSave

def main():
    with open("path.txt", 'r') as path:
        data_dir = path.read()    
    compression(data_dir)

if __name__ == "__main__":
    main()