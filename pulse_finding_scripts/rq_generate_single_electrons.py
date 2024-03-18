import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mpl
import time
import sys
import glob
from natsort import natsorted
#from collections import deque

#from scipy.signal import spectogram

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


def find_single_electrons(data_dir, handscan=False, max_pulses=2, filtered=True, simpleS2=True, save_avg_wfm=False, phase="liquid"):
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


    #spe_sizes = np.concatenate((spe_sizes_0,spe_sizes_1,spe_sizes_2))

    spe_dir = "/home/xaber/crystalize/Analysis/spe_calibration/202403/50V_3-6-2024.txt"
    spe_sizes = np.loadtxt(spe_dir, dtype='float')
    print(spe_sizes)

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

    
    """
    se_area = np.zeros((n_events,20))
    se_width = np.zeros((n_events,20))
    se_coincidence = np.zeros((n_events,20))
    se_max_height = np.zeros((n_events,20))
    se_rms = np.zeros(n_events)
    se_area_ch = np.zeros((n_events, n_channels-1))
    se_tba = np.zeros(n_events)
    se_top_x = np.zeros(n_events)
    se_top_y = np.zeros(n_events)

    s1_area = np.zeros(n_events)
    s2_area = np.zeros(n_events)

    se_a2_height = np.zeros((n_events,20))
    se_aft10 = np.zeros((n_events,20))
    se_aft90 = np.zeros((n_events,20))

    drift_time = np.zeros(n_events)

    right_area = np.zeros(n_events)
    """

    gate_start = np.zeros(n_events, dtype=int)
    gate_end = np.zeros(n_events, dtype=int)
    gate_max_loc = np.zeros(n_events)
    gate_area = np.zeros(n_events)
    gate_height = np.zeros(n_events)
    gate_rms = np.zeros(n_events)

    cathode_start = np.zeros(n_events, dtype=int)
    cathode_end = np.zeros(n_events, dtype=int)
    cathode_area = np.zeros(n_events)
    cathode_height = np.zeros(n_events)

    n_pulses = np.zeros(n_events, dtype=int)
    p_start = np.zeros((n_events,4),dtype=int)
    p_end = np.zeros((n_events,4), dtype=int)

    s1_area = np.zeros(n_events)



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
    

        """
        # Better baseline subtraction
        # Super slow
        ch_data_adcc = np.reshape(ch_data_adcc, (n_sipms,n_events_b,wsize))
        for s in range(n_sipms):
            for t in range(n_events_b):
                suppressed = np.nonzero( (ch_data_adcc[s,t,:] == 0) )[0]
                for u in range(suppressed.size - 1):
                    if suppressed[u+1] - suppressed[u] > 1:
                        baseline = np.mean(ch_data_adcc[s,t,suppressed[u]+1:suppressed[u]+1+75])
                        ch_data_adcc[s,t,suppressed[u]+1:suppressed[u+1]-1] -= baseline


                    #if ch_data_adcc[u] == 0 and ch_data_adcc[u+1] != 0:
        """
                        
        ch_data_adcc = np.reshape(ch_data_adcc, int(n_sipms*n_events_b*wsize))
        




        # Convert from ADCC to phd/sample and get summed waveform
        ch_data_adcc = np.concatenate((ch_data_adcc, np.zeros(n_tot_samp_per_ch) ))
        ch_data_mV = vscale*np.reshape(ch_data_adcc, (n_channels,n_events_b,wsize))
        ch_data_phdPerSample = np.zeros_like(ch_data_mV) 
        for ch in range(n_sipms): # using np.divide to get rid of this for loop is more trouble than it's worth
            if spe_sizes[ch] == 0: continue
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


            # Pulse finder
            start_times, end_times = vs.find_pulses_in_event(ch_data_phdPerSample[-1,i-j*block_size,:], max_pulses=max_pulses, verbose=False)

            # Sort pulses by start times
            startinds = np.argsort(start_times)
            n_pulses[i] = min(max_pulses,len(start_times))
            mp = 0
            for m in startinds:
                if m >= max_pulses:
                    continue
                if start_times[m] < 0.25/tscale: continue
                p_start[i,mp] = start_times[m]
                p_end[i,mp] = end_times[m]
                mp += 1

            # If drift time is too short, continue
            if (p_start[i,1] - p_start[i,0])*tscale < 3.8: continue
                
            # First pulse should be an S1
            temp_area = pq.GetPulseArea(p_start[i,0], p_end[i,0], ch_data_phdPerSample[-1,i-j*block_size,:] )
            (p_afs_2l, p_afs_1, p_afs_10, p_afs_25, p_afs_50, p_afs_75, p_afs_90, p_afs_99) = pq.GetAreaFraction(p_start[i,0], p_end[i,0], ch_data_phdPerSample[-1,i-j*block_size,:] )
            if (p_afs_50-p_afs_2l)*tscale > 0.125 or temp_area < 10: continue
            s1_area[i] = temp_area
            s1_start_afs2 = p_afs_2l

            # Make first window, find max
            
            temp_start = int(s1_start_afs2 + 0.8/tscale)
            temp_end = int(temp_start + (1.9)/tscale)
            gate_max_loc[i] = np.argmax(ch_data_phdPerSample[-1,i-j*block_size,temp_start:temp_end]) + temp_start
            if temp_end - gate_max_loc[i] < 0.7/tscale or gate_max_loc[i] - temp_start < 0.7/tscale: continue 

            """
            gate_max_loc[i] = np.argmax(ch_data_phdPerSample[-1,i-j*block_size,temp_start:temp_end]) + temp_start
            #if gate_max_loc[i] - temp_start > 0.8: continue
            if temp_end - gate_max_loc[i] > 0.75/tscale:
                gate_start[i] = temp_start
                gate_end[i] = int(gate_max_loc[i] + 0.75/tscale)
            """
            gate_start[i] = int(gate_max_loc[i] - 0.7/tscale)
            gate_end[i] = int(gate_max_loc[i] + 0.7/tscale)

            if gate_end[i] == gate_start[i] or gate_start[i] == 0 or gate_end[i] == 0: continue

            # "Precise" baseline subtract from right before the S1
            base_precise_end = p_start[i,0] - 50
            base_precise_start = p_start[i,0] - 50 - 100
            bases = np.sum( ch_data_phdPerSample[:,i-j*block_size,base_precise_start:base_precise_end], axis=1)/(base_precise_end-base_precise_start)
            for ch in range(32):
                ch_data_phdPerSample[ch,i-j*block_size,:] -= bases[ch]
            ch_data_phdPerSample[-1,i-j*block_size,:] = np.sum(ch_data_phdPerSample[:-1,i-j*block_size,:], axis=0)

            # RQ's
            gate_area[i] = pq.GetPulseArea(gate_start[i], gate_end[i], ch_data_phdPerSample[-1,i-j*block_size,:] )
            gate_height[i] = max(ch_data_phdPerSample[-1,i-j*block_size,gate_start[i]:gate_end[i]])
            gate_rms[i] = np.sqrt(np.sum(np.power(ch_data_phdPerSample[-1,i-j*block_size,gate_start[i]:gate_end[i]],2))/(gate_end[i]-gate_start[i]))
            #(p_afs_2l[i], p_afs_1, p_afs_10, p_afs_25, p_afs_50, p_afs_75, p_afs_90, p_afs_99) = pq.GetAreaFraction(gate_start[i], gate_end[i], ch_data_phdPerSample[-1,i-j*block_size,:] )

            """
            # Find an S1
            temp_area = np.zeros(max_pulses)
            s1_afs = -1
            for pp in range(max_pulses):
                temp_area[pp] = pq.GetPulseArea(p_start[i,pp], p_end[i,pp], ch_data_phdPerSample[-1,i-j*block_size,:] )
                (p_afs_2l, p_afs_1, p_afs_10, p_afs_25, p_afs_50, p_afs_75, p_afs_90, p_afs_99) = pq.GetAreaFraction(p_start[i,pp], p_end[i,pp], ch_data_phdPerSample[-1,i-j*block_size,:] )
                if (p_afs_50-p_afs_2l)*tscale < 0.125 and temp_area[pp] > 10:
                    s1_afs = p_afs_2l
            if s1_afs == -1: continue

            if (p_start[i,1] - p_start[i,0])*tscale < 3.8 or temp_area[0] < 10: continue

            # If found S1, look for SE's in the S1 echos
            try:
                gate_start[i] = int(s1_afs + 0.8/tscale)
                gate_end[i] = int(s1_afs + (0.8+2.8)/tscale)
                gate_area[i] = pq.GetPulseArea(gate_start[i], gate_end[i], ch_data_phdPerSample[-1,i-j*block_size,:] )
                gate_height[i] = max(ch_data_phdPerSample[-1,i-j*block_size,gate_start[i]:gate_end[i]])
                gate_rms[i] = np.sqrt(np.sum(np.power(ch_data_phdPerSample[-1,i-j*block_size,gate_start[i]:gate_end[i]],2))/(gate_end[i]-gate_start[i]))

                cathode_start[i] = int(s1_afs + 5/tscale)
                cathode_end[i] = int(s1_afs + (5+2.8)/tscale)
                cathode_area[i] = pq.GetPulseArea(cathode_start[i], cathode_end[i], ch_data_phdPerSample[-1,i-j*block_size,:] )
                cathode_height[i] = max(ch_data_phdPerSample[-1,i-j*block_size,cathode_start[i]:cathode_end[i]])

                s1_area[i] = temp_area[0]

            except:
                print(f"error ev {i}")
                continue

            
            plot_debug = True


            delayed = True


            

            #max_loc = np.argmax(ch_data_phdPerSample[-1,i-j*block_size,:])
            #if max_loc < 1000 or wsize - max_loc < 1000: continue


            #ps_baseline = np.mean(ch_data_phdPerSample[-1,i-j*block_size,max_loc-600:max_loc-100])
            if delayed:

                
                a1 = s1_filter(ch_data_phdPerSample[-1,i-j*block_size,:],120)

                
                a2 = s2_filter(ch_data_phdPerSample[-1,i-j*block_size,:],a1,700)
                
                
                se_start = 0
                se_end = 0
                max_val_a2 = max(a2)
                max_loc = np.argmax(a2)
                thresh = 1e-8
                peaks, properties = signal.find_peaks(a2,height=0.5,distance=800,width=10)
                se_start = np.zeros(10,dtype=int)
                se_end = np.zeros(10,dtype=int)
                for p in range(min(10,len(peaks))):
                    try:
                        se_start[p] = peaks[p] - min(peaks[p] - (ch_data_phdPerSample[-1,i-j*block_size,:peaks[p]] == 0).nonzero()[0] )
                        se_end[p] = 2*peaks[p] +  min(  (ch_data_phdPerSample[-1,i-j*block_size,peaks[p]:] == 0).nonzero()[0] - peaks[p] )
                    except:
                        continue
                    
                    se_area[i,p] = np.sum(ch_data_phdPerSample[-1,i-j*block_size,se_start[p]:se_end[p]])
                    afs_2l,afs_1,afs_10,afs_25,afs_50,afs_75,afs_90,afs_99 = pq.GetAreaFraction(se_start[p],se_end[p],ch_data_phdPerSample[-1,i-j*block_size,:])
                    se_aft10[i,p] = afs_10
                    se_aft90[i,p] = afs_90 
                    se_width[i,p] = se_end[p] - se_start[p]
                    se_a2_height[i,p] = a2[peaks[p]]

            else:

                # New pulse finder
                start_times = []
                end_times = []
                areas = []
                lh_cut = wsize
                for g in range(2):
                    temp_start, temp_end = vs.PulseFinderVerySimple(ch_data_phdPerSample[-1,i-j*block_size,:lh_cut], verbose=False)
                    if temp_start != temp_end:
                        if g==0: right_area[i] = np.sum(ch_data_phdPerSample[-1,i-j*block_size,temp_end:])
                        start_times.append(temp_start)
                        end_times.append(temp_end)
                        areas.append(np.sum(ch_data_phdPerSample[-1,i-j*block_size,temp_start:temp_end]))
                        lh_cut = temp_start
                    else:
                        continue
                
                nFound = len(start_times)
                if nFound == 0: continue

                afs_2l,afs_1,afs_10,afs_25,afs_50,afs_75,afs_90,afs_99 = pq.GetAreaFraction(start_times[0],end_times[0],ch_data_phdPerSample[-1,i-j*block_size,:])
                this_right_area = np.sum(  ch_data_phdPerSample[-1,i-j*block_size,end_times[0]:] )


                b_s1 = (tscale*(afs_50 - afs_2l) < 0.15)&(this_right_area < 200)&(areas[0] < 5000)

                max_wf_loc = np.argmax(ch_data_phdPerSample[-1,i-j*block_size,:]) 

                #b_s1s2 = (nFound == 2) and (areas[0] > 1000) and (areas[1] > 10) and (end_times[1] - start_times[0] > 1000  )
                #b_s1 = (nFound == 1) and (areas[0] > 100) and (areas[0] < 5000) and (end_times[0] < wsize - 2000)

                #if b_s1s2:
                #    wf = ch_data_phdPerSample[-1,i-j*block_size,end_times[1]:start_times[0]]
                #    offset = end_times[1]

                se_start = np.zeros(10,dtype=int)
                se_end = np.zeros(10,dtype=int)
                if b_s1:

                    se_start[0] = max_wf_loc + 2/tscale
                    se_end[0] = se_start[0] + 1.2/tscale

                    se_start[1] = max_wf_loc + 6/tscale
                    se_end[1] = se_start[1] + 1.2/tscale   

                    
                    se_start[0] = 11.5/tscale 
                    se_end[0] = 11.5/tscale + 800

                    se_start[1] = 15/tscale 
                    se_end[1] = 15/tscale + 800
                    
                    offset = 0
                    #wf = ch_data_phdPerSample[-1,i-j*block_size,end_times[0]:]
                    #offset = end_times[0]
                #else:
                    #continue


                


                #se_start, se_end, a1, a2 = se_pf(wf)
                #se_start, se_end = se_pf_dumb(wf)
                #se_start += offset
                #se_end += offset

                #for p in range(2):
                #    se_area[i,p] = np.sum( ch_data_phdPerSample[-1,i-j*block_size, se_start[p]:se_end[p]] )
                #    se_max_height[i,p] = max(ch_data_phdPerSample[-1,i-j*block_size, se_start[p]:se_end[p]])
                    #se_area[i,p] = np.sum(wf[se_start[p]:se_end[p]])
                    #se_width[i,p] = tscale*(se_end[p]-se_start[p])



            


            
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
            plotyn = False #(gate_area[i] > 20 and gate_area[i] < 40) #True #True #np.any( a2 > 7  )
            #plotyn = np.any( (se_area[i,:] > 0)&(se_area[i,:] < 2)&(se_width[i,:] > 1.3/tscale)&(se_coincidence[i,:] == 3)  )
            if not inn == 'q' and plotyn and plot_event_ind == i:
                #print((se_area[i,:] > 0)&(se_area[i,:] < 2)&(se_width[i,:] > 1.3/tscale)&(se_coincidence[i,:] == 3))

                fspect, tspect, Sxx = signal.spectrogram(ch_data_phdPerSample[-1,i-j*block_size,:], fs=500e6)
                #print(np.round(se_area[i,:],1))
                fig = pl.figure()
                ax = pl.gca()
                #pl.pcolormesh(tspect, fspect, np.log10(Sxx), shading="gouraud",)
                pl.plot(x*tscale, ch_data_phdPerSample[-1,i-j*block_size,:],color='black',lw=1.2, label = "Summed All" )
                for ch in range(32):
                    pl.plot(x*tscale, ch_data_phdPerSample[ch,i-j*block_size,:] )

                pl.axvspan(tscale*(gate_start[i]),tscale*(gate_end[i]),alpha=0.2,color="blue",zorder=0)
                pl.axvspan(tscale*(base_precise_start),tscale*(base_precise_end),alpha=0.2,color="orange",zorder=0)
                #pl.axvspan(tscale*(cathode_start[i]),tscale*(cathode_end[i]),alpha=0.2,color="blue",zorder=0)
                pl.annotate(f"{round(gate_area[i],2)} phd", xy=(tscale*(gate_start[i])+1.5, 0.03 + gate_height[i]/ax.get_ylim()[1]), fontsize=15)
                #pl.plot(x*tscale, a2,color='red')
                #pl.plot(x*tscale, a1,color='green')
                #pl.plot(all_peaks*tscale, ch_data_phdPerSample[-1,i-j*block_size,all_peaks] ,"ro")
                #pl.ylim(-1,5)
                #for ch in range(32):
                #    pl.plot(x*tscale, ch_data_phdPerSample[ch,i-j*block_size,:])
                #pl.plot(all_peaks*tscale, np.zeros_like(all_peaks) ,"ro")
                #pl.plot(x*tscale,  wf_filtered, color="red", lw=0.7)
                #pl.plot(x*tscale, data_conv,"blue", label="S2 Filtered")
                #pl.axvspan(tscale*(se_start[p]),tscale*(se_end[p]),alpha=0.3,color="blue",zorder=0)
                #for p in range(10):
                #    pl.axvspan(tscale*(offset + se_start[p]),tscale*(offset + se_end[p]),alpha=0.2,color="blue",zorder=0)
                #    pl.axvspan(tscale*( se_start[p]),tscale*(se_end[p]),alpha=0.2,color="blue",zorder=0)
                #ax.text(tscale*(start_window+se_end), 0.03 + se_max_height[i]/ax.get_ylim()[1], '{:.2f} phd'.format(se_area[i]), fontsize=9, color="blue")
                pl.xlabel(r"Time [$\mu$s]")
                pl.ylabel("phd/sample")
                #pl.ylim(-0.02,0.1)
                pl.title("Event {}".format(i))  
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
        #print(np.count_nonzero(se_coincidence))
        j += 1

    # end of loop over compressed files


    # ==========================================================================================================================
    # Final step: save rq's

    #create a dictionary with all RQs
    list_rq = {}
    #print(np.count_nonzero(se_area))

    list_rq['gate_start'] = gate_start
    list_rq['gate_end'] = gate_end
    list_rq['gate_area'] = gate_area
    list_rq['gate_height'] = gate_height
    list_rq['gate_rms'] = gate_rms

    list_rq['cathode_start'] = cathode_start
    list_rq['cathode_end'] = cathode_end
    list_rq['cathode_area'] = cathode_area
    list_rq['cathode_height'] = cathode_height

    list_rq['s1_area'] = s1_area
    

    """
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
    list_rq['se_aft10'] = se_aft10
    list_rq['se_aft90'] = se_aft90
    """


    #list_rq[''] =    #add more rq

    #remove zeros in the end of each RQ array. 
    #for rq in list_rq.keys():
    #    if rq != 'n_events':
    #        list_rq[rq] = list_rq[rq][:n_events]

    if filtered:
        save_name = "/rq_SE_v1.npy"
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
    """
    Boxcar filter for S1s

    Inputs:
      wf: Waveform to filter
      n1: Filter width
    """

    n1_12 = int(n1/2)
    a1 = np.zeros_like(wf)

    a1[n1_12] = np.sum(wf[0:2*n1_12])
    for i in range(n1_12+1,wf.size-n1_12):
        a1[i] = a1[i-1] + wf[i+n1_12-1] - wf[i-n1_12-1]
 
    return a1


def s2_filter(wf,a1,n2):
    """
    Boxcar filter with S1 filter subtracted.
    This uses a sliding window maximum finder algorithm on a1

    Inputs:
      wf: Waveform to filter
      a1: S1 filtered waveform
      n2: Filter width
    """

    n2_12 = int(n2/2)
    a1 = np.round(a1,8) # saw some floating point errors
    a2 = np.zeros_like(wf)

    # Initial A2 sample before loop
    window_a1 = np.asarray( sorted(a1[0:2*n2_12]) )
    max_a1 = window_a1[-1]
    a2[n2_12] = np.sum(wf[0:2*n2_12]) - max_a1
    old_max_a1 = max_a1

    # Loop over samples
    for i in range(n2_12+1,wf.size-n2_12):
        
        # Remove A1 sample not in current window
        to_remove = a1[i-n2_12-1]
        ind_to_remove = (window_a1 == to_remove).nonzero()[0]
        if ind_to_remove.size > 0: window_a1[ind_to_remove[0]] = -99999

        


        """
        # Add new A1 sample in window, determine max A1
        to_change = a1[i+n2_12-1]
        ind_to_change = (window_a1 < to_change).nonzero()[0][-1]
        window_a1[:ind_to_change+1] = -99999
        window_a1[ind_to_change] = to_change
        
        window_a1 = window_a1[window_a1 != -99999]
        
        max_a1 = window_a1[-1]
        """
 

        """
        # Add new A1 sample in window, determine max A1
        to_change = a1[i+n2_12-1]
        ind_to_change = (window_a1 < to_change).nonzero()[0]
        if ind_to_change.size > 0:
            window_a1[ind_to_change] = -999999999
            window_a1[ind_to_change[-1]] = to_change
        max_a1 = window_a1[-1]
        """

        
        # Slower version
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



def se_pf(wf):

    a1 = s1_filter(wf,120)
    a2 = s2_filter(wf,a1,700)

    se_start = 0
    se_end = 0
    peaks, properties = signal.find_peaks(a2,height=0.5,distance=800,width=10)
    se_start = np.zeros(10,dtype=int)
    se_end = np.zeros(10,dtype=int)
    for p in range(min(10,len(peaks))):
        try:
            se_start[p] = peaks[p] - min(peaks[p] - (wf[:peaks[p]] == 0).nonzero()[0] )
            se_end[p] = 2*peaks[p] +  min(  (wf[peaks[p]:] == 0).nonzero()[0] - peaks[p] )
        except:
            continue

    
    return se_start, se_end, a1, a2



def se_pf_dumb(wf):

    se_start = 0
    se_end = 0
    peaks, properties = signal.find_peaks(wf,height=0.05,distance=800,width=10)
    se_start = np.zeros(10,dtype=int)
    se_end = np.zeros(10,dtype=int)
    for p in range(min(10,len(peaks))):
        try:
            se_start[p] = peaks[p] - min(peaks[p] - (wf[:peaks[p]] == 0).nonzero()[0] )
            se_end[p] = 2*peaks[p] +  min(  (wf[peaks[p]:] == 0).nonzero()[0] - peaks[p] )
        except:
            continue


    return se_start, se_end




def main():
    #with open(sys.path[0]+"/path.txt", 'r') as path:
    #    data_dir = path.read()
   
    #data_dir_list = glob.glob("/media/xaber/outSSD2/crystalize_data/data-202403/20240309/*")
    data_dir_list = ["/media/xaber/extradrive1/crystalize_data/data-202403/20240315/20240315-091729/"]

    data_dir_list = glob.glob("/media/xaber/extradrive1/crystalize_data/data-202403/20240315/*/")
    # next
    #"/media/xaber/extradrive1/crystalize_data/data-202403/20240315/20240315-091729/"

    #data_dir_list = ["/media/xaber/outSSD2/crystalize_data/data-202403/20240309/20240309-184356/"]

    for data_dir in data_dir_list:
        print(data_dir)
        try:
            find_single_electrons(data_dir, phase="liquid")
        except:
            print("error")
            continue

    #data_dir = "/media/xaber/outSSD2/crystalize_data/data-202403/20240306/20240306-181727/"
    #find_single_electrons(data_dir, phase="liquid")

if __name__ == "__main__":
    main()