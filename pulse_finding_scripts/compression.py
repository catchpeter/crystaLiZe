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
    filtered=True, 
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
        if filtered:
            save_dir = data_dir+"compressed_filtered_data_dev/"
        else:
            save_dir = data_dir+"compressed_data_dev/"
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
        tosser = [] 
        lastToCheck = int(np.max(evNum) )
        firstToCheck = int(np.min(evNum) ) 
        for i in range(firstToCheck,lastToCheck+1):
            if np.count_nonzero(evNum == i) != 3 and np.count_nonzero(evNum == i) > 0: tosser.append(i)
        if verbose:
            print("Number of misaligned events: "+str(len(tosser)))
        
        
        # Initializing arrays to store data
        # Note: np.zeros is preferred over np.empty because we want default to be zero
        all_bool_front = np.zeros((n_all_ch, max_n_events, wsize-8), dtype=bool)
        all_bool_back = np.zeros((n_all_ch, max_n_events, wsize-8), dtype=bool)
        raw_data_front = np.zeros((n_all_ch,max_n_events,wsize-8))
        filtered_data_front = np.zeros((n_all_ch,max_n_events,wsize-8))
        data_to_save = np.zeros((n_all_ch,max_n_events,wsize-8))


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
                        bool_front = baseline_suppress(ch_data_BSF_filtered, pls_thresh=threshold, buffL=200, buffR=300) #, buffL=0, buffR=100, condense_thresh=0 )
                        bool_back = baseline_suppress(ch_data_BSB_filtered, pls_thresh=threshold, buffL=200, buffR=300)
                        all_bool_front[ch_ind,ev-toss_counts,tot_delay:] = bool_front
                        all_bool_back[ch_ind,ev-toss_counts,tot_delay:] = bool_back
                        

                        if debug: 
                            t = 2*np.arange(wsize-8-tot_delay)
                            pl.figure()
                            pl.plot(t, ch_data_BSF, "black")
                            pl.plot(t, ch_data_BSF_filtered-100, "red")
                            pl.plot(t, threshold*np.ones_like(t) - 100, "cyan" )
                            pl.plot(t, 100*bool_front, "green")
                            pl.xlabel("Time [ns]")
                            pl.ylabel("ADCC")
                            pl.title("Channel 0, Event "+str(ev))
                            pl.grid("major","major")
                            pl.show()


                        # Save header info for the first channel
                        if ch_ind == 0: headers[ev-toss_counts+bk*block_size,:] = ch_data[ev,0:8]
                        
                    # end of event loop

                ch_ind += 1
                # end of channel loop
        

        # Decide on what to save 
        to_save_or_not_to_save = np.logical_and(all_bool_front, all_bool_back)
        if filtered:
            data_to_save[to_save_or_not_to_save] = filtered_data_front[to_save_or_not_to_save]
        else:
            data_to_save[to_save_or_not_to_save] = raw_data_front[to_save_or_not_to_save]


        # Some plotting for debugging
        if debug:
            t = 2*np.arange(wsize-8)
            for ev in range(n_events[0]):
                pl.figure()
                pl.plot(t, raw_data_front[0,ev,:], "black")
                pl.plot(t, filtered_data_front[0,ev,:] - 100, "red")
                pl.plot(t, threshold*np.ones_like(t) - 100, "cyan" )
                pl.plot(t, data_to_save[0,ev,:], "blue", alpha=0.5)
                pl.plot(t, 100*to_save_or_not_to_save[0,ev,:], "green")
                #pl.plot(t, 200*all_bool_front[0,ev,:], "green")
                #pl.legend(("Raw","Filtered (shifted down)","Saved raw"))
                pl.xlabel("Time [ns]")
                pl.ylabel("ADCC")
                pl.title("Channel 0, Event "+str(ev))
                pl.grid("major","major")
                pl.show()
        

        
        # Save that mf
        if save_mode == "npy":
            if filtered:
                np.savez_compressed(f'{save_dir}compressed_filtered_dev{bk}', data_to_save.flatten())
            else:
                np.savez_compressed(f'{save_dir}compressed_{bk}', data_to_save.flatten())
            
        elif save_mode in ("none","None", None):
            print("Disabled saving, moving to next block")
            #return all_data_front #stuffToSave

    # Headers are saved at the very end!
    if save_mode not in ("none","None", None): 
        np.savez_compressed(f'{save_dir}headers', headers.flatten())


def main():
    with open("path.txt", 'r') as path:
        data_dir = path.read()    
    compression(data_dir)

if __name__ == "__main__":
    main()