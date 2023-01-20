import numpy as np
import time
import os
try:
    import h5py
except ImportError:
    pass
import matplotlib.pyplot as pl

from read_settings import get_event_window
from ch_evt_filter_compress import baseline_suppress, filter_channel_event




def compression(
    data_dir,
    threshold=5, 
    save_mode="npy", 
    save_everything=False, 
    ret_block='all',
    channel='all', 
    tellblocks=False,
    verbose=True,
    debug=False):
    """
    Checks board alignment, does zero-baseline suppression, compresses data
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
            if verbose:
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
        if verbose:
            print("Total events: "+str(tot_ev) )
            print("Number of compressed files = "+str(tot_fi))
        time.sleep(1)
    
    headers = np.zeros((tot_ev+10,8),dtype=int)
    
    if ret_block in ("all", "All", "ALL"):
        tot_fi_to_loop = range(tot_fi)
    else:
        tot_fi_to_loop = [ret_block]
        if verbose:
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
        tosser = [] # tosser? I barely knew her!
        lastToCheck = int(np.max(evNum) )
        firstToCheck = int(np.min(evNum) ) 
        for i in range(firstToCheck,lastToCheck+1):
            if np.count_nonzero(evNum == i) != 3 and np.count_nonzero(evNum == i) > 0: tosser.append(i)
        if verbose:
            print("Number of misaligned events: "+str(len(tosser)))
        
        
        # Initializing arrays to store data
        # Note: np.zeros is preferred over np.empty because we want default to be zero
        all_bool_front = np.zeros((n_all_ch, max_n_events, wsize-8))
        all_bool_back = np.zeros((n_all_ch, max_n_events, wsize-8))
        #all_data_front = np.zeros((n_all_ch, max_n_events, wsize-8))
        #all_data_back = np.zeros((n_all_ch, max_n_events, wsize-8))
        raw_data_front = np.zeros((n_all_ch,max_n_events,wsize-8))
        filtered_data_front = np.zeros((n_all_ch,max_n_events,wsize-8))
        data_to_save = np.zeros((n_all_ch,max_n_events,wsize-8))

        use_aarons_baseline_suppress = True

        # Load data one channel at a time
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
                
                # Loop over events to check alignment, correct for delay between boards, then do baseline quashing    
                toss_counts = 0
                for ev in range(n_events[bd]):
                    if ch_data[ev, 2] in tosser:
                        toss_counts += 1
                    else:
                        # Get which board to determine time delay
                        if bd == 0: tot_delay = 0
                        elif bd == 1: tot_delay = delay
                        elif bd == 2: tot_delay = 2*delay
                            
                        # Subtract baseline in two ways
                        ch_data_BSF = ch_data[ev, 8:(wsize-tot_delay)] - np.mean(ch_data[ev, 8:8+150])
                        ch_data_BSB = ch_data[ev, 8:(wsize-tot_delay)] - np.mean(ch_data[ev, -150:-1])
                        raw_data_front[ch_ind,ev-toss_counts,tot_delay:] = ch_data_BSF

                        # Filter data
                        ch_data_BSF_filtered = filter_channel_event(ch_data_BSF)
                        ch_data_BSB_filtered = filter_channel_event(ch_data_BSB)
                        filtered_data_front[ch_ind,ev-toss_counts,tot_delay:] = ch_data_BSF_filtered

                        # Get bool for baseline suppression
                        if use_aarons_baseline_suppress:
                            all_bool_front[ch_ind,ev-toss_counts,tot_delay:] = baseline_suppress(ch_data_BSF_filtered, pls_thresh=threshold)
                            all_bool_back[ch_ind,ev-toss_counts,tot_delay:] = baseline_suppress(ch_data_BSB_filtered, pls_thresh=threshold)
                        else:
                            all_bool_front[ch_ind,ev-toss_counts,tot_delay:] = np.absolute(ch_data_BSF_filtered) > threshold
                            all_bool_back[ch_ind,ev-toss_counts,tot_delay:] = np.absolute(ch_data_BSB_filtered) > threshold

                        # Save header info for the first channel
                        if ch_ind == 0: headers[ev-toss_counts+bk*block_size,:] = ch_data[ev,0:8]
                        
                    # end of event loop

                ch_ind += 1
                # end of channel loop
        

        # Decide on what to save 
        bool_front = np.sum(all_bool_front,axis=0)
        bool_front[bool_front > 0] = np.ones_like(bool_front[bool_front > 0]) 
        bool_back = np.sum(all_bool_back,axis=0)
        bool_back[bool_back > 0] = np.ones_like(bool_back[bool_back > 0])
        
        to_save_or_not_to_save = np.logical_and(bool_front,bool_back)


        # Do some extra work if not using Aaron's function
        if not use_aarons_baseline_suppress:

            # Buffers before and after stuff passing threshold
            # Not particuarly efficient, but this only runs once per block
            # Assumption: buffL <= buffR (which is fine, before a pulse is baseline, after a pulse could be something else)
            buffL=100
            buffR=250
            for i in range(buffR):
                diff_bool = np.diff(to_save_or_not_to_save, axis=1)
                if i < buffL: to_save_or_not_to_save[:,:-1][diff_bool != 0] = 1 # left buffers
                to_save_or_not_to_save[:,1:][diff_bool != 0] = 1 # right buffers

            # Condense. Even more inefficient
            condense_thresh = 100
            for ev in range(n_events[0]):
                save_ind = np.nonzero(to_save_or_not_to_save[ev,:] == 1)[0]
                diff_save_ind = np.diff(save_ind)
                for i in range(diff_save_ind.size):
                    if diff_save_ind[i] < condense_thresh and diff_save_ind[i] > 1:
                        to_save_or_not_to_save[ev,save_ind[i]:int(min(save_ind[i]+condense_thresh,wsize-8))] = 1





        data_to_save[:,to_save_or_not_to_save] = raw_data_front[:,to_save_or_not_to_save]


        # Some plotting for debugging
        if debug:
            t = 2*np.arange(wsize-8)
            for ev in range(n_events[0]):
                pl.figure()
                pl.plot(t, raw_data_front[0,ev,:], "black")
                pl.plot(t, filtered_data_front[0,ev,:] - 200, "red")
                pl.plot(t, data_to_save[0,ev,:], "blue", alpha=0.5)
                pl.show()
        

        
        # Save that mf
        if save_mode == "npy":
            1 == 2
            #np.savez_compressed(f'{save_dir}compressed_{bk}', all_data_front.flatten())
            #np.savez_compressed(f'{save_dir}headers', headers.flatten())
        elif save_mode in ("none","None", None):
            print("Disabled saving, moving to next block")
            #return all_data_front #stuffToSave

def main():
    with open("path.txt", 'r') as path:
        data_dir = path.read()    
    compression(data_dir)

if __name__ == "__main__":
    main()