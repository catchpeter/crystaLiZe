import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mpl
import time
import sys
import glob
from natsort import natsorted
#from collections import deque

#from scipy.signal import spectogram

import PulseFinderScipy as pf
import PulseQuantities as pq
import PulseClassification as pc
import PulseFinderVerySimple as vs
from read_settings import get_event_window, get_vscale, get_sipm_bias

from ch_evt_filter_compress import baseline_suppress, filter_channel_event

from compression import compression

from scipy import signal
from scipy.signal import find_peaks

#import scipy.signal.spectrogram as spectrogram
#from scipy.signal import spectrogram
#from c_process import S2filter
#from ch_evt_filter_compress import filter_channel_event





def find_SPEs(data_dir, handscan=False, max_pulses=5, filtered=True, correct_swap=False):
    # ====================================================================================================================================
    # Plotting parameters

    # set plotting style
    mpl.rcParams['font.size']=10
    mpl.rcParams['legend.fontsize']='small'
    mpl.rcParams['figure.autolayout']=True
    mpl.rcParams['figure.figsize']=[8.0,6.0]

    # use for coloring pulses
    pulse_class_colors = np.array(['blue', 'green', 'red', 'magenta', 'darkorange'])
    pulse_class_labels = np.array(['Other', 'S1-like LXe', 'S1-like gas', 'S2-like', 'Merged S1/S2'])
    pc_legend_handles=[]
    for class_ind in range(len(pulse_class_labels)):
        pc_legend_handles.append(mpl.patches.Patch(color=pulse_class_colors[class_ind], label=str(class_ind)+": "+pulse_class_labels[class_ind]))

    inn="" # used to control hand scan


    # ====================================================================================================================================
    # define DAQ and SiPM parameters
    
    vscale = get_vscale(data_dir) # ADCC to mV conversion
    
    # Timing parameters
    event_window = get_event_window(data_dir)
    if event_window < 0: 
        print("Invalid event window")
        return
    wsize = int(500 * event_window)  # samples per waveform # 12500 for 25 us
    tscale = (8.0/4096.0)     # = 0.002 µs/sample, time scale
    post_trigger = 0.5 # Was 0.2 for data before 11/22/19
    trigger_time_us = event_window*(1-post_trigger)
    trigger_time = int(trigger_time_us/tscale)

    # SiPM numbering
    n_sipms = 32    
    n_channels = n_sipms + 1 # include sum
    n_top = int((n_channels-1)/2)
    top_channels=np.arange(0,16)  #np.array(range(n_top),int)
    bottom_channels=np.arange(16,32) #np.array(range(n_top,2*n_top),int)
    
    block_size = int(1500*15/event_window) # number of events per compressed file

    l_window = 120
    r_window = 280

    #if avg_wf and spe_height[k,i,0] > 0.16 and spe_height[k,i,0] < 0.67: # 54
    #if avg_wf and spe_height[k,i,0] > 0.16 and spe_height[k,i,0] < 0.6: # 53
    #if avg_wf and spe_height[k,i,0] > 0.16 and spe_height[k,i,0] < 0.5: # 52
    #if avg_wf and spe_height[k,i,0] > 0.16 and spe_height[k,i,0] < 0.42: # 51
    #if avg_wf and spe_height[k,i,0] > 0.16 and spe_height[k,i,0] < 0.34: # 50

    sipm_bias = get_sipm_bias(data_dir)
    if sipm_bias == 54: spe_upper_height = 0.67
    if sipm_bias == 53: spe_upper_height = 0.6
    if sipm_bias == 52: spe_upper_height = 0.5
    if sipm_bias == 51: spe_upper_height = 0.42
    if sipm_bias == 50: spe_upper_height = 0.34
    if sipm_bias == 49: spe_upper_height = 0.25
    else: spe_upper_height = 99999

    print(sipm_bias,spe_upper_height)

    # ====================================================================================================================================
    # Configure header data and event time

    # Get list of compressed files and calculate number of events
    if filtered:
        compressed_file_list = natsorted(glob.glob(data_dir+"/compressed_filtered*/compressed_*.npz") )
    else:
        compressed_file_list = natsorted(glob.glob(data_dir+"/compressed_data/compressed_*.npz") )
    if len(compressed_file_list) < 1:
        print("No compressed files found in "+data_dir)
        return

    n_events = (len(compressed_file_list)+5)*block_size # with some extra room 

    # ====================================================================================================================================
    # Initialize rq's to save
    # np.zeros is preferred over np.empty bc/ we want zero to be default value

    
    spe_area = np.zeros((n_sipms,n_events,max_pulses))
    spe_height = np.zeros((n_sipms,n_events,max_pulses))
    spe_height_time = np.zeros((n_sipms,n_events,max_pulses))
    spe_rms = np.zeros((n_sipms,n_events,max_pulses))
    spe_nPulses = np.zeros((n_sipms,n_events),dtype=int)
    spe_sum_wf = np.zeros((n_sipms,n_events))

    spe_sum = np.zeros((n_sipms,1000))
    spe_nSaved = np.zeros(n_sipms)

    # ====================================================================================================================================
    # Begin loop over compressed files

    j = 0 # index for number of compressed files
    counter = 0 # index for total events

    for compressed_file in compressed_file_list:
        #if j > 2: break
        # load data
        try:
            with np.load(compressed_file) as data:
                ch_data_adcc = data["arr_0"]
        except:
            print("Error in loading "+compressed_file)
            continue
        
        n_tot_samp_per_ch = int( (ch_data_adcc.size)/n_sipms )
        n_events_b = int((ch_data_adcc.size)/(n_sipms*wsize)) # n events per compressed file (same as block_size)
    
        #dothebetterbaselinesubtraction = asdf

        """
        # Better baseline subtraction
        # Super slow
        ch_data_adcc = np.reshape(ch_data_adcc, (n_sipms,n_events_b,wsize))
        for s in range(n_sipms):
            print(j,str(s+1)+"/32")
            for t in range(n_events_b):
                suppressed = np.nonzero( (ch_data_adcc[s,t,:] == 0) )[0]
                for u in range(suppressed.size - 1):
                    if suppressed[u+1] - suppressed[u] > 100:
                        baseline = np.mean(ch_data_adcc[s,t,suppressed[u]+1:suppressed[u]+1+75])
                        ch_data_adcc[s,t,suppressed[u]+1:suppressed[u+1]-1] -= baseline

                        spe_area.append(tscale*1000*vscale*np.sum( ch_data_adcc[s,t,suppressed[u]+1+75:suppressed[u+1]-1]) )
                        spe_height.append(vscale*max(ch_data_adcc[s,t,suppressed[u]+1+75:suppressed[u+1]-1]) )
                        spe_rms.append( get_rms(vscale*ch_data_adcc[s,t,suppressed[u]+1+75:suppressed[u+1]-1]) )
                        spe_width.append(tscale*1000*(suppressed[u+1]-1 - suppressed[u]+1+75)  )


                    #if ch_data_adcc[u] == 0 and ch_data_adcc[u+1] != 0:
        
                        
        ch_data_adcc = np.reshape(ch_data_adcc, int(n_sipms*n_events_b*wsize))
        """




        # Convert from ADCC to phd/sample and get summed waveform
        ch_data_adcc = np.concatenate((ch_data_adcc, np.zeros(n_tot_samp_per_ch) ))
        ch_data_mV = vscale*np.reshape(ch_data_adcc, (n_channels,n_events_b,wsize))
        ch_data_mV[-1,:,:] = np.sum(ch_data_mV, axis=0)

        if correct_swap:
            print("\nCorrecting for swapped channels\n")
            temp_31 = ch_data_mV[7,:,:]
            temp_7 = ch_data_mV[31,:,:]
            ch_data_mV[7,:,:] = temp_7
            ch_data_mV[31,:,:] = temp_31



        
        # create a time axis in units of µs:
        x = np.arange(0, wsize, 1)
        t = tscale*x
        n_events = int(ch_data_mV[0].shape[0])
        if n_events == 0: break
            
        # ====================================================================================================================================
        # Loop over events
        
        avg_wf = True

        for i in range(j*block_size, j*block_size+n_events):

            # Loop over sipms
            for k in range(n_sipms):

                all_max_i = np.zeros(max_pulses, dtype=int)

                wf = np.copy(ch_data_mV[k,i-j*block_size,:])

                for r in range(max_pulses):
                    max_i, code = analyze_biggest_pulse(wf, l_window, r_window)
                    if code == 0:
                        if r > 0: 
                            if np.count_nonzero( np.absolute(max_i - all_max_i) < l_window + r_window) == 0:
                                all_max_i[r] = max_i
                        elif r == 0:
                            all_max_i[r] = max_i
                        wf[max_i-l_window:max_i+r_window] = 0
                    elif code == 1:
                        wf[:max_i+r_window] = 0
                    elif code == 2:
                        wf[max_i-l_window:] = 0


                spe_sum_wf[k,i] = np.sum(ch_data_mV[k,i-j*block_size,:])
                spe_nPulses[k,i] = np.count_nonzero(all_max_i > 0) # haha
                all_max_i[:spe_nPulses[k,i]] = all_max_i[all_max_i > 0]
                for ps in range(spe_nPulses[k,i]):

                    spe_height_time[k,i,ps] = all_max_i[ps]
                    spe_height[k,i,ps] = ch_data_mV[k,i-j*block_size,all_max_i[ps]]
                    spe_area[k,i,ps] = tscale*1000*np.sum(ch_data_mV[k,i-j*block_size,all_max_i[ps]-l_window:all_max_i[ps]+r_window])
                    spe_rms[k,i,ps] = get_rms(ch_data_mV[k,i-j*block_size,all_max_i[ps]-l_window:all_max_i[ps]+r_window])








          

            # ==========================================================================================================================
            # Plotting code, resist the urge to use this for handscanning blobs

            # Code to allow skipping to another event index for plotting
            plot_event_ind = i
            try:
                plot_event_ind = int(inn)
                if plot_event_ind < i:
                    inn = ''
                    plot_event_ind = i
                    print("Can't go backwards! Continuing to next event.")
            except ValueError:
                plot_event_ind = i

            # Condition to skip the individual plotting, hand scan condition
            if np.size(handscan) == 1: 
                plotyn = handscan
            else:
                plotyn = handscan[i]
           
            if inn == 's': sys.exit()

            if not inn == 'q' and plotyn and plot_event_ind == i:

 
                fig = pl.figure()
                ax = pl.gca()
                for ch in range(0,31):
                    pl.plot(x*tscale, ch_data_mV[ch,i-j*block_size,:] )
                    #for ps in range(spe_nPulses[0,i]):
                    #    #pl.axvspan(tscale*(spe_height_time[ch,i,ps]-l_window), tscale*(spe_height_time[ch,i,ps]+r_window),alpha=0.2,color="blue",zorder=0  )

  
  
                pl.xlabel(r"Time [$\mu$s]")
                pl.ylabel("mV")
                #pl.ylim(-0.02,0.1)
                pl.title("Event {}".format(i))  
                #pl.legend()
                pl.grid(which="both",axis="both",linestyle="--")
                #pl.xlim(0,event_window)
                #pl.ylim(-2*max(ch_data_mV[-1,i-j*block_size,:]), 5*max(ch_data_mV[-1,i-j*block_size,:]))
                pl.show()
                inn = input("Press enter to continue, q to stop plotting, evt # to skip to # (forward only)")
                #fig.clf()

            counter += 1
                
        # end of loop over events
        


        #n_events = i
        t_end = time.time()
        print("total number of events processed:", n_events)
        #print("Time used: {}".format(t_end-t_start))
        j += 1

    # end of loop over compressed files


    # ==========================================================================================================================
    # Final step: save rq's

    #create a dictionary with all RQs
    list_rq = {}
    list_rq['spe_area'] = spe_area
    list_rq['spe_height'] = spe_height
    list_rq['spe_height_time'] = spe_height_time
    list_rq['spe_rms'] = spe_rms
    list_rq['spe_sum'] = spe_sum 
    list_rq['spe_nSaved'] = spe_nSaved
    list_rq['spe_sum_wf'] = spe_sum_wf
   


    if filtered:
        save_name = "/rq_SPE_filtered_new.npy"
    else:
        save_name = "/rq_SPE_test.npy"
    rq = open(data_dir + save_name,'wb')
    np.savez(rq, **list_rq)
    rq.close()

  


def analyze_biggest_pulse(wf, l_window, r_window):

    wf_s = wf.size
    if wf_s < 2: return -99999, 4

    max_i = -99999
    code = 0

    max_i = np.argmax(wf)
    if max_i > l_window+5 and max_i < wf_s-r_window-5:
        area = np.sum(wf[max_i-l_window:max_i+r_window])
        rms = get_rms(wf[max_i-l_window:max_i+r_window])
    if max_i < l_window+5: code += 1
    if max_i > wf_s-r_window-5: code += 2

    return max_i, code




def get_rms(wf):

    return np.sqrt(np.sum(np.power(wf,2) )/wf.size )




def main():
    #with open(sys.path[0]+"/path.txt", 'r') as path:
    #    data_dir = path.read()
   

    print("Sleeping before rq-ing")
    #time.sleep(4*60*60)
    #data_dir_list = glob.glob("/media/xaber/G-Drive2/crystalize_data/data-202307/20230710/*49SiPM*SPE*60min/")
    data_dir_list = glob.glob("/media/xaber/G-Drive2/crystalize_data/data-202307/20230711/*50SiPM*SPE*60min/")

    data_dir_list = glob.glob("/media/xaber/G-Drive2/crystalize_data/data-202307/20230720/*SPE*/")

    print("\n")
    for i in data_dir_list: print(i)
    print("\n")
    #time.sleep(9.5*60*60)

    for data_dir in data_dir_list:
        print(data_dir)
        #compression(data_dir, save_everything = True, debug=False)
        find_SPEs(data_dir)

if __name__ == "__main__":
    main()