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
from read_settings import get_event_window, get_vscale

from ch_evt_filter_compress import baseline_suppress, filter_channel_event

from scipy import signal
#import scipy.signal.spectrogram as spectrogram
#from scipy.signal import spectrogram
#from c_process import S2filter
#from ch_evt_filter_compress import filter_channel_event


def find_single_electrons(data_dir, handscan=False, max_pulses=4, filtered=True, simpleS2=True, save_avg_wfm=False, phase="liquid"):
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

    # To do: create txt files of SPE sizes to read in

    # SPE sizes, as of April 2022
    #spe_sizes_0 = np.array([85.406,86.876,84.763,83.986,85.470,85.032,85,968,85.452,84.126,84.825,84.340,85.217,84.285,85.226,83.753,84.609])
    #spe_sizes_1 = np.array([79.897,78.625,81.818,81.189,74.952,77.289,79.880,76.970])
    #spe_sizes_2 = np.array([79.597,79.023,80.213,81.023,78.173,79.883,79.069,75.496])

    # SPE sizes, as of May 24, 2022
    #spe_sizes_0 = np.array([83.325,85.437,84.449,83.025,84.827,84.184,85.503,84.656,85.029,84.984,84.961,84.562,84.014,85.917,83.846,86.926])
    #spe_sizes_1 = np.array([79.617,79.891,81.203,80.859,75.585,77.097,79.112,78.256])
    #spe_sizes_2 = np.array([79.232,78.950,78.912,80.191,79.115,77.171,73.503,78.465])

    # SPE sizes in SOLID, June 15, 2022
    #spe_sizes_0 = np.array([89.524,87.773,85.071,86.106,86.781,86.102,86.417,86.840,86.919,86.675,85.908,86.630,85.881,87.946,87.629,88.216])
    #spe_sizes_1 = np.array([84.134,83.419,83.957,84.037,79.348,81.228,82.678,81.617])
    #spe_sizes_2 = np.array([82.342,82.528,82.477,82.523,84.209,81.481,78.693,81.669])

    if phase == "liquid":
        # SPE sizes in LIQUID, Sept, 28, 2022
        spe_sizes_0 = np.array([90.010,88.944,88.831,88.209,90.296,91.249,93.823,92.238,87.540,89.733,86.509,83.149,90.761,91.263,91.641,93.016])
        spe_sizes_1 = np.array([92.661,93.194,92.746,94.623,89.254,94.524,93.302,93.410])
        spe_sizes_2 = np.array([92.955,93.619,93.944,93.199,95.470,96.461,92.317,93.509])

    elif phase == "solid":
        # SPE sizes in SOLID, Oct 2022
        spe_sizes_0 = np.array([90.836,93.329,90.721,90.831,93.071,91.682,93.485,95.265,88.747,91.275,88.771,89.520,93.875,94.136,94.966,94.632])
        spe_sizes_1 = np.array([94.666,95.533,93.915,99.042,97.783,94.895,97.134,97.501])
        spe_sizes_2 = np.array([94.553,95.514,96.554,96.465,96.711,96.920,95.460,95.705])

    spe_sizes = np.concatenate((spe_sizes_0,spe_sizes_1,spe_sizes_2))


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

    

    se_area = np.zeros((n_events,20))
    se_width = np.zeros((n_events,20))
    se_coincidence = np.zeros((n_events,20))
    se_max_height = np.zeros(n_events)
    se_rms = np.zeros(n_events)
    se_area_ch = np.zeros((n_events, n_channels-1))
    se_tba = np.zeros(n_events)
    se_top_x = np.zeros(n_events)
    se_top_y = np.zeros(n_events)

    s1_area = np.zeros(n_events)
    s2_area = np.zeros(n_events)

    se_a2_height = np.zeros((n_events,20))

    drift_time = np.zeros(n_events)

    right_area = np.zeros(n_events)


    # ====================================================================================================================================
    # Begin loop over compressed files

    j = 0 # index for number of compressed files
    counter = 0 # index for total events

    for compressed_file in compressed_file_list:
        #if j > 10: break
        # load data
        try:
            with np.load(compressed_file) as data:
                ch_data_adcc = data["arr_0"]
        except:
            print("Error in loading "+compressed_file)
            continue
        
        n_tot_samp_per_ch = int( (ch_data_adcc.size)/n_sipms )
        n_events_b = int((ch_data_adcc.size)/(n_sipms*wsize)) # n events per compressed file (same as block_size)
    
        # Convert from ADCC to phd/sample and get summed waveform
        ch_data_adcc = np.concatenate((ch_data_adcc, np.zeros(n_tot_samp_per_ch) ))
        ch_data_mV = vscale*np.reshape(ch_data_adcc, (n_channels,n_events_b,wsize))
        ch_data_phdPerSample = np.zeros_like(ch_data_mV) 
        for ch in range(n_sipms): # using np.divide to get rid of this for loop is more trouble than it's worth
            ch_data_phdPerSample[ch,:,:] = ch_data_mV[ch,:,:]*tscale*(1000)/spe_sizes[ch]
        ch_data_phdPerSample[-1,:,:] = np.sum(ch_data_phdPerSample, axis=0)   
        #ch_data_phdPerSample = np.concatenate((ch_data_phdPerSample, np.sum(ch_data_phdPerSample, axis=0)))
            
        # create a time axis in units of µs:
        x = np.arange(0, wsize, 1)
        t = tscale*x
        n_events = int(ch_data_phdPerSample[0].shape[0])
        if n_events == 0: break
            
       
        # ====================================================================================================================================
        # Loop over events

        for i in range(j*block_size, j*block_size+n_events):

            plot_debug = True

            """
            start_flag = False
            ps_found = 0
            se_start = 0
            se_end = 0
            for s in range(wsize):

                if ps_found > 19: continue

                if not start_flag and ch_data_phdPerSample[-1,i-j*block_size,s] != 0:
                    start_flag = True
                    se_start = s


                elif start_flag and ch_data_phdPerSample[-1,i-j*block_size,s] == 0:    
                    start_flag = False
                    se_end = s
                    se_area[i,ps_found] = np.sum(ch_data_phdPerSample[-1,i-j*block_size,se_start:se_end])
                    se_width[i,ps_found] = se_end-se_start
                    se_coincidence[i,ps_found] = np.count_nonzero( np.sum(ch_data_phdPerSample[:-1,i-j*block_size,se_start:se_end], axis=1) != 0  )
                    print(se_coincidence[i,ps_found])

                    ps_found += 1

            """




            """

            se_start  = 0
            se_end = 0

            #for bl in range(4):
            #    bl_s = int(wsize/4)
            #    ch_data_phdPerSample[-1,i-j*block_size,i*bl_s:(i+1)*bl_s] -= np.median(ch_data_phdPerSample[-1,i-j*block_size,i*bl_s:(i+1)*bl_s])

            all_peaks = np.array([],dtype=int)
            all_ind = np.array([],dtype=int)
            for ch in range(32):
                peaks, properties = signal.find_peaks(ch_data_phdPerSample[ch,i-j*block_size,:],height=0.015,distance=200,width=10)
                ind = ch*np.ones_like(peaks)

                all_peaks = np.concatenate((all_peaks,peaks),dtype=int)
                all_ind = np.concatenate((all_ind,ind),dtype=int)

            sorted_i = np.argsort(all_peaks)
            sorted_peaks = all_peaks[sorted_i]
            sorted_ind = all_ind[sorted_i]

            current_se = []
            se_start = 0
            se_end = 0
            nFound = 0
            for p in range(0,sorted_peaks.size-3):
                if p in current_se: continue
                peaks_in_range = (sorted_peaks < sorted_peaks[p] + 900)&(sorted_peaks >= sorted_peaks[p])
                peaks_before = (sorted_peaks > sorted_peaks[p] - 900)&(sorted_peaks < sorted_peaks[p])
                peaks_after = (sorted_peaks < sorted_peaks[p] + 900 + 900)&(sorted_peaks > sorted_peaks[p] + 900)
                if np.count_nonzero(peaks_in_range) > 2 and np.count_nonzero(peaks_before) == 0 and np.count_nonzero(peaks_after) == 0 and np.ptp( sorted_peaks[peaks_in_range] ) > 30:
                    se_start = sorted_peaks[p] - 50
                    se_end = sorted_peaks[peaks_in_range][-1] + 200
                    se_width[i,nFound] = se_end - se_start
                    current_se = sorted_peaks[peaks_in_range]

                    contributing_ch = sorted_ind[peaks_in_range]
                    for cch in contributing_ch:
                        se_area[i,nFound] += np.sum(ch_data_phdPerSample[cch,i-j*block_size,se_start:se_end])
                    nFound += 1
                    if nFound > 3: break

                    #break


            
            """



            #max_loc = np.argmax(ch_data_phdPerSample[-1,i-j*block_size,:])
            #if max_loc < 1000 or wsize - max_loc < 1000: continue


            #ps_baseline = np.mean(ch_data_phdPerSample[-1,i-j*block_size,max_loc-600:max_loc-100])

            a1 = s1_filter(ch_data_phdPerSample[-1,i-j*block_size,:],120)

            a2 = s2_filter(ch_data_phdPerSample[-1,i-j*block_size,:],a1,800)
            
            se_start = 0
            se_end = 0
            max_val_a2 = max(a2)
            max_loc = np.argmax(a2)
            thresh = 1e-8
            peaks, properties = signal.find_peaks(a2,height=1.5,distance=800,width=10)
            se_start = np.zeros(20,dtype=int)
            se_end = np.zeros(20,dtype=int)
            for p in range(min(20,len(peaks))):
                try:
                    se_start[p] = peaks[p] - min(peaks[p] - (ch_data_phdPerSample[-1,i-j*block_size,:peaks[p]] == 0).nonzero()[0] )
                    se_end[p] = 2*peaks[p] +  min(  (ch_data_phdPerSample[-1,i-j*block_size,peaks[p]:] == 0).nonzero()[0] - peaks[p] )
                except:
                    continue
                
                se_area[i,p] = np.sum(ch_data_phdPerSample[-1,i-j*block_size,se_start[p]:se_end[p]])
                se_width[i,p] = se_end[p] - se_start[p]
                se_a2_height[i,p] = a2[peaks[p]]

 


            """
            # New pulse finder
            start_times = []
            end_times = []
            areas = []
            lh_cut = wsize
            for g in range(max_pulses):
                temp_start, temp_end = vs.PulseFinderVerySimple(ch_data_phdPerSample[-1,i-j*block_size,:lh_cut], verbose=False)
                if temp_start != temp_end:
                    if g==0: right_area[i] = np.sum(ch_data_phdPerSample[-1,i-j*block_size,temp_end:])
                    start_times.append(temp_start)
                    end_times.append(temp_end)
                    areas.append(np.sum(ch_data_phdPerSample[-1,i-j*block_size,temp_start:temp_end]))
                    lh_cut = temp_start
                else:
                    continue

            """

            

            


            """
            if len(start_times) > 1 and areas[1] > 5:
                temp_start = end_times[1] + 300
                temp_end = start_times[0] - 100
                start_window = temp_start
                if temp_end - temp_start < 1000: continue
            else:
                continue
            
            elif len(start_times) == 1 and areas[0] < 5000:
                temp_start = end_times[0] + 300
                temp_end = wsize - 200
                start_window = temp_start
                if temp_end - temp_start < 1000: continue
            
            

            data_window = ch_data_phdPerSample[-1,i-j*block_size,temp_start:temp_end]
            data_window_ch = ch_data_phdPerSample[:-1,i-j*block_size,temp_start:temp_end]

            # Find points above threshold
            # If there exists a cluster, tag it 
            se_thresh = 2*0.015
            above_thresh_i = (data_window > se_thresh).nonzero()[0]
            se_start = 0
            se_end = 0
            values, edges = np.histogram(above_thresh_i,bins=int(wsize*2/200),range=(0,wsize)  )
            for b in range(1,values.size-4):
                if values[b] > 0 and values[b+1] > 0 and values[b+2] > 0 and values[b+3] > 0:
                    se_start = int(edges[b]) - 200
                    for c in range(values.size-4-b):
                        if values[b+3+c] == 0:
                            se_end = int(edges[b+3+c]) + 200
                            break
            
            if se_end - se_start <= 0: continue # or se_end - se_start > 2000: continue

            se_area[i] = np.sum(data_window[se_start:se_end])
            #print(se_area[i])
            se_width[i] = se_end - se_start
            se_max_height[i] = max(data_window[se_start:se_end])
            se_rms[i] = np.sqrt(np.sum(np.power(data_window[se_start:se_end],2))/(data_window[se_start:se_end].size))

            se_area_ch[i,:] = pq.GetPulseAreaChannel(se_start, se_end, data_window_ch[:,se_start:se_end])
            se_area_top = np.sum(se_area_ch[i,:16])
            se_area_bottom = np.sum(se_area_ch[i,16:])
            se_tba[i] = (se_area_top-se_area_bottom)/(se_area_top+se_area_bottom)
            se_bot_x, se_bot_y, se_top_x[i], se_top_y[i] = pq.GetCentroids(se_area_ch[i,:])

            """





            """
            data_window = ch_data_phdPerSample[-1,i-j*block_size,se_start:se_end]
            data_window_ch = ch_data_phdPerSample[:-1,i-j*block_size,se_start:se_end]
            se_area[i] = np.sum(data_window)
            #print(se_area[i])
            se_width[i] = se_end - se_start
            se_max_height[i] = max(data_window)
            se_rms[i] = np.sqrt(np.sum(np.power(data_window,2))/(data_window.size))

            se_area_ch[i,:] = pq.GetPulseAreaChannel(se_start, se_end, data_window_ch)
            se_area_top = np.sum(se_area_ch[i,:16])
            se_area_bottom = np.sum(se_area_ch[i,16:])
            se_tba[i] = (se_area_top-se_area_bottom)/(se_area_top+se_area_bottom)
            se_bot_x, se_bot_y, se_top_x[i], se_top_y[i] = pq.GetCentroids(se_area_ch[i,:])
            """
            


            """
            if len(start_times) > 0:
                area2 = np.sum(ch_data_phdPerSample[-1,i-j*block_size,start_times[0]:end_times[0] ])

                if area2 < 5000:
                    start_window = end_times[0] + 1500
                    if wsize - start_window < 1000: continue # 3000

                    data_window = ch_data_phdPerSample[-1,i-j*block_size,start_window:]
                    data_window_ch = ch_data_phdPerSample[:-1,i-j*block_size,start_window:]

                    max_loc = np.argmax(data_window)
                    max_val = max(data_window)

                    try:
                        left_zeros = (data_window[:max_loc] == 0).nonzero()[0]
                        se_start = max_loc - min(max_loc - left_zeros) # + 100
                        right_zeros = (data_window[max_loc:] == 0).nonzero()[0]
                        se_end = 2*max_loc + min(right_zeros - max_loc) #- 200
                    except:
                        continue

                    if se_start != se_end:

                        se_area[i] = np.sum(data_window[se_start:se_end])
                        #print(se_area[i])
                        se_width[i] = se_end - se_start
                        se_max_height[i] = max_val
                        se_rms[i] = np.sqrt(np.sum(np.power(data_window,2))/(data_window.size))

                        se_area_ch[i,:] = pq.GetPulseAreaChannel(se_start, se_end, data_window_ch)
                        se_area_top = np.sum(se_area_ch[i,:16])
                        se_area_bottom = np.sum(se_area_ch[i,16:])
                        se_tba[i] = (se_area_top-se_area_bottom)/(se_area_top+se_area_bottom)
                        se_bot_x, se_bot_y, se_top_x[i], se_top_y[i] = pq.GetCentroids(se_area_ch[i,:])

                        #s1_area[i] = np.sum(ch_data_phdPerSample[-1,i-j*block_size,start_times[1]:end_times[1] ])
                        #s2_area[i] = np.sum(ch_data_phdPerSample[-1,i-j*block_size,start_times[0]:end_times[0] ])
                        #drift_time[i] = tscale*(start_times[0]-start_times[1])
                
                    else: continue
                else: continue
            else: continue
            """


            """
            if len(start_times) > 1 and False:
                start_window = end_times[1] + 200 #200
                end_window = start_times[0] - 200   #200
                if end_window - start_window < 600:
                    continue
                else:
                    data_window = ch_data_phdPerSample[-1,i-j*block_size,start_window:end_window]
                    data_window_ch = ch_data_phdPerSample[:-1,i-j*block_size,start_window:end_window]

                    max_loc = np.argmax(data_window)
                    max_val = max(data_window)

                    try:
                        left_zeros = (data_window[:max_loc] == 0).nonzero()[0]
                        se_start = max_loc - min(max_loc - left_zeros) # + 100
                        right_zeros = (data_window[max_loc:] == 0).nonzero()[0]
                        se_end = 2*max_loc + min(right_zeros - max_loc) #- 200
                    except:
                        continue
                    
                    if se_start != se_end:

                        se_area[i] = np.sum(data_window[se_start:se_end])
                        se_width[i] = se_end - se_start
                        se_max_height[i] = max_val
                        se_rms[i] = np.sqrt(np.sum(np.power(data_window,2))/(data_window.size))

                        se_area_ch[i,:] = pq.GetPulseAreaChannel(se_start, se_end, data_window_ch)
                        se_area_top = np.sum(se_area_ch[i,:16])
                        se_area_bottom = np.sum(se_area_ch[i,16:])
                        se_tba[i] = (se_area_top-se_area_bottom)/(se_area_top+se_area_bottom)
                        se_bot_x, se_bot_y, se_top_x[i], se_top_y[i] = pq.GetCentroids(se_area_ch[i,:])

                        s1_area[i] = np.sum(ch_data_phdPerSample[-1,i-j*block_size,start_times[1]:end_times[1] ])
                        s2_area[i] = np.sum(ch_data_phdPerSample[-1,i-j*block_size,start_times[0]:end_times[0] ])
                        drift_time[i] = tscale*(start_times[0]-start_times[1])

                #continue
            else:
                continue
            """



          

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
            #plotyn = se_area[i] > 15 and se_area[i] < 25
            #print(se_area[i,:])
            plotyn = True #True #np.any( a2 > 7  )
            #plotyn = np.any( (se_area[i,:] > 0)&(se_area[i,:] < 2)&(se_width[i,:] > 1.3/tscale)&(se_coincidence[i,:] == 3)  )
            if not inn == 'q' and plotyn and plot_event_ind == i:
                #print((se_area[i,:] > 0)&(se_area[i,:] < 2)&(se_width[i,:] > 1.3/tscale)&(se_coincidence[i,:] == 3))

                fspect, tspect, Sxx = signal.spectrogram(ch_data_phdPerSample[-1,i-j*block_size,:], fs=500e6)

                fig = pl.figure()
                ax = pl.gca()
                #pl.pcolormesh(tspect, fspect, np.log10(Sxx), shading="gouraud",)
                pl.plot(x*tscale, ch_data_phdPerSample[-1,i-j*block_size,:],color='black',lw=1.2, label = "Summed All" )
                pl.plot(x*tscale, a2,color='red')
                #pl.plot(x*tscale, a1,color='green')
                #pl.plot(all_peaks*tscale, ch_data_phdPerSample[-1,i-j*block_size,all_peaks] ,"ro")
                pl.ylim(-1,5)
                #for ch in range(32):
                #    pl.plot(x*tscale, ch_data_phdPerSample[ch,i-j*block_size,:])
                #pl.plot(all_peaks*tscale, np.zeros_like(all_peaks) ,"ro")
                #pl.plot(x*tscale,  wf_filtered, color="red", lw=0.7)
                #pl.plot(x*tscale, data_conv,"blue", label="S2 Filtered")
                #pl.axvspan(tscale*(start_window+se_start),tscale*(start_window+se_end),alpha=0.3,color="blue",zorder=0)
                for p in range(20):
                    pl.axvspan(tscale*(se_start[p]),tscale*(se_end[p]),alpha=0.2,color="blue",zorder=0)
                #ax.text(tscale*(start_window+se_end), 0.03 + se_max_height[i]/ax.get_ylim()[1], '{:.2f} phd'.format(se_area[i]), fontsize=9, color="blue")
                pl.xlabel(r"Time [$\mu$s]")
                pl.ylabel("phd/sample")
                #pl.title("Event {}".format(i))  
                #pl.legend()
                pl.grid(which="both",axis="both",linestyle="--")
                #pl.xlim(0,event_window)
                #pl.ylim(-2*max(ch_data_phdPerSample[-1,i-j*block_size,:]), 5*max(ch_data_phdPerSample[-1,i-j*block_size,:]))
                pl.show()
                inn = input("Press enter to continue, q to stop plotting, evt # to skip to # (forward only)")
                #fig.clf()

            counter += 1
                
        # end of loop over events


        n_events = i
        t_end = time.time()
        print("total number of events processed:", n_events)
        #print("Time used: {}".format(t_end-t_start))
        print(np.count_nonzero(se_coincidence))
        j += 1

    # end of loop over compressed files


    # ==========================================================================================================================
    # Final step: save rq's

    #create a dictionary with all RQs
    list_rq = {}
    print(np.count_nonzero(se_area))
    list_rq['se_area'] = se_area
    list_rq['se_width'] = se_width
    list_rq['se_max_height'] = se_max_height
    list_rq['se_rms'] = se_rms
    list_rq['se_area_ch'] = se_area_ch
    list_rq['se_tba'] = se_tba
    list_rq['se_top_x'] = se_top_x
    list_rq['se_top_y'] = se_top_y
    list_rq['se_coincidence'] = se_coincidence
    list_rq['se_a2_height'] = se_a2_height
    
    list_rq['s1_area'] = s1_area
    list_rq['s2_area'] = s2_area
    list_rq['drift_time'] = drift_time


    #list_rq[''] =    #add more rq

    #remove zeros in the end of each RQ array. 
    #for rq in list_rq.keys():
    #    if rq != 'n_events':
    #        list_rq[rq] = list_rq[rq][:n_events]

    if filtered:
        save_name = "/rq_SE_filtered_test2.npy"
    else:
        save_name = "/rq_SE_test.npy"
    rq = open(data_dir + save_name,'wb')
    np.savez(rq, **list_rq)
    rq.close()

  

def high_pass_filter(wf,alpha):
    wf_filtered = np.zeros_like(wf)
    wf_filtered[0] = wf[0]
    for i in range(1,wf.size):
        wf_filtered[i] = alpha*(wf[i] - wf[i-1] + wf_filtered[i-1])
    
    return wf_filtered



def s1_filter(wf,n1):

    n1_12 = int(n1/2)
    a1 = np.zeros_like(wf)

    a1[n1_12] = np.sum(wf[0:2*n1_12])
    for i in range(n1_12+1,wf.size-n1_12):
        a1[i] = a1[i-1] + wf[i+n1_12-1] - wf[i-n1_12-1]
 
    return a1


def s2_filter(wf,a1,n2):

    n2_12 = int(n2/2)
    a1 = np.round(a1,8)
    a2 = np.zeros_like(wf)

    window_a1 = np.asarray( sorted(a1[0:2*n2_12]) )
    max_a1 = window_a1[-1]
    a2[n2_12] = np.sum(wf[0:2*n2_12]) - max_a1
    old_max_a1 = max_a1

    # Main loop over samples
    for i in range(n2_12+1,wf.size-n2_12):
        
        # Remove sample not in current window
        to_remove = a1[i-n2_12-1]
        ind_to_remove = (window_a1 == to_remove).nonzero()[0]
        if ind_to_remove.size > 0: window_a1[ind_to_remove[0]] = -999999999

        # Determining max A1 in window
        window_a1 = window_a1[window_a1 >= a1[i+n2_12-1]]
        window_a1 = np.concatenate(([a1[i+n2_12-1]],window_a1))
        max_a1 = window_a1[-1]

        # Calculating A2
        a2[i] = a2[i-1] + old_max_a1 + wf[i+n2_12-1] - wf[i-n2_12-1] - max_a1

        # Reset old max for next iteration
        old_max_a1 = max_a1
     

    return a2


def test_filter(wf,n12):

    aTest = np.zeros_like(wf)

    for i in range(n12,wf.size-n12):
        aTest[i] = np.count_nonzero( wf[i-n12:i+n12] > 0.01)

    return aTest




def main():
    #with open(sys.path[0]+"/path.txt", 'r') as path:
    #    data_dir = path.read()
    #data_dir = "/media/xaber/f5d91b31-9a7d-3278-ac5b-4f9ae16edd60/crystalize_data/data-202211/20221120/20221120-1302_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_54SiPM_1.44bar_-94.64ICVbot_2fold_beta_120min/"
    #data_dir = "/media/xaber/f5d91b31-9a7d-3278-ac5b-4f9ae16edd60/crystalize_data/data-202303/20230307/20230307-1744_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_54SiPM_1.67bar_-99.78ICVbot_2fold_blank_afterFlow_60min/"
    data_dir = "/media/xaber/extradrive2/crystalize_data/data-202303/20230320/20230320-0254_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_54SiPM_1.58bar_77.51ICVbot_2fold_degradedNew_60min/"
    find_single_electrons(data_dir)

if __name__ == "__main__":
    main()