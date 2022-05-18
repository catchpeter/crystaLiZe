
import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mpl
import time
import sys
import glob

import PulseFinderScipy as pf
import PulseQuantities as pq
import PulseClassification as pc

#data_dir = "G:/.shortcut-targets-by-id/11qeqHWCbcKfFYFQgvytKem8rulQCTpj8/crystalize/data/data-202103/031121/Po_2.8g_3.0c_0.78bar_circ_30min_1312/"
#data_dir = "/home/xaber/caen/wavedump-3.8.2/data/041921/Po_2.8g_3.0c_0.72bar_circ_20min_0928/"
#data_dir = "G:/My Drive/crystalize/data/data-202104/041421/Po_2.8g_3.0c_1.1bar_circ_60min_1747/"

def make_rq(data_dir, handscan = False):   
    # set plotting style
    mpl.rcParams['font.size']=10
    mpl.rcParams['legend.fontsize']='small'
    mpl.rcParams['figure.autolayout']=True
    mpl.rcParams['figure.figsize']=[8.0,6.0]

    # ==================================================================
    # define DAQ and other parameters
    
    # Get vscale
    if data_dir.find("0.5DR") != -1:
        vscale = (500.0/16384.0)
    elif data_dir.find("2DR") != -1:
        vscale = (2000.0/16384.0) # = 0.122 mV/ADCC, vertical scale
    else:
        vscale = (2000.0/16384.0) # default to 2V

    # Get window size
    if data_dir.find("3us") != -1:
        event_window = 3
    elif data_dir.find("6us") != -1:
        event_window = 6
    elif data_dir.find("15us") != -1:
        event_window = 15
    elif data_dir.find("25us") != -1:
        event_window = 25
    else:
        print("Need to input window size")
        return
    

    wsize = int(500 * event_window)  # samples per waveform # 12500 for 25 us
    tscale = (8.0/4096.0)     # = 0.002 µs/sample, time scale

    save_avg_wfm = False # get the average waveform passing some cut and save to file

    post_trigger = 0.5 # Was 0.2 for data before 11/22/19
    trigger_time_us = event_window*(1-post_trigger)
    trigger_time = int(trigger_time_us/tscale)


    n_sipms = 32
    n_channels = n_sipms + 1 # include sum
    
    block_size = 1500


    # define top, bottom channels
    n_top = int((n_channels-1)/2)
    top_channels=np.array(range(n_top),int)
    bottom_channels=np.array(range(n_top,2*n_top),int)

    # SPE sizes, as of spring 2022
    spe_sizes_0 = np.array([85.406,86.876,84.763,83.986,85.470,85.032,85,968,85.452,84.126,84.825,84.340,85.217,84.285,85.226,83.753,84.609])
    spe_sizes_1 = np.array([79.897,78.625,81.818,81.189,74.952,77.289,79.880,76.970])
    spe_sizes_2 = np.array([79.597,79.023,80.213,81.023,78.173,79.883,79.069,75.496])
    spe_sizes = np.concatenate((spe_sizes_0,spe_sizes_1,spe_sizes_2))


    # ==================================================================

    
    compressed_file_list = sorted(glob.glob(data_dir+"compressed_data/compressed_*.npy") )
    if len(compressed_file_list) < 1:
        print("No compressed files found in "+data_dir+"compressed_data/")
        return

    
    n_events = (len(compressed_file_list)+5)*1500 # some extra room 

    
    
    
    
    # RQ's to save
    
    # RQs to add:
    # Pulse level: channel areas (fracs; max fracs), TBA, rise time? (just difference of AFTs...)
    # Event level: drift time; S1, S2 area
    # Pulse class (S1, S2, other)
    # max number of pulses per event
    max_pulses = 4
    p_start = np.zeros(( n_events, max_pulses), dtype=int)
    p_end   = np.zeros(( n_events, max_pulses), dtype=int)
    p_found = np.zeros(( n_events, max_pulses), dtype=int)

    #center of mass
    center_top_x = np.zeros(( n_events, max_pulses))
    center_top_y = np.zeros(( n_events, max_pulses))
    center_bot_x = np.zeros(( n_events, max_pulses))
    center_bot_y = np.zeros(( n_events, max_pulses))

    p_area = np.zeros(( n_events, max_pulses))
    p_max_height = np.zeros(( n_events, max_pulses))
    p_min_height = np.zeros(( n_events, max_pulses))
    p_width = np.zeros(( n_events, max_pulses))

    p_afs_2l = np.zeros((n_events, max_pulses) )
    p_afs_2r = np.zeros((n_events, max_pulses) )
    p_afs_1 = np.zeros((n_events, max_pulses) )
    p_afs_25 = np.zeros((n_events, max_pulses) )
    p_afs_50 = np.zeros((n_events, max_pulses) )
    p_afs_75 = np.zeros((n_events, max_pulses) )
    p_afs_99 = np.zeros((n_events, max_pulses) )
                
    p_hfs_10l = np.zeros((n_events, max_pulses) )
    p_hfs_50l = np.zeros((n_events, max_pulses) )
    p_hfs_10r = np.zeros((n_events, max_pulses) )
    p_hfs_50r = np.zeros((n_events, max_pulses) )

    p_mean_time = np.zeros((n_events, max_pulses) )
    p_rms_time = np.zeros((n_events, max_pulses) )

    # Channel level (per event, per pulse, per channel)
    p_start_ch = np.zeros((n_events, max_pulses, n_channels-1), dtype=int)
    p_end_ch = np.zeros((n_events, max_pulses, n_channels-1), dtype=int )
    p_area_ch = np.zeros((n_events, max_pulses, n_channels-1) )
    p_area_ch_frac = np.zeros((n_events, max_pulses, n_channels-1) )

    p_area_top = np.zeros((n_events, max_pulses))
    p_area_bottom = np.zeros((n_events, max_pulses))
    p_tba = np.zeros((n_events, max_pulses))

    p_class = np.zeros((n_events, max_pulses), dtype=int)

    # Event-level variables
    n_pulses = np.zeros(n_events, dtype=int)

    n_s1 = np.zeros(n_events, dtype=int)
    n_s2 = np.zeros(n_events, dtype=int)
    sum_s1_area = np.zeros(n_events)
    sum_s2_area = np.zeros(n_events)
    drift_Time = np.zeros(n_events)
    drift_Time_AF = np.zeros(n_events)
    drift_Time_AF50 = np.zeros(n_events)
    drift_Time_AF150 = np.zeros(n_events)
    drift_Time_AS = np.zeros(n_events) # for multi-scatter drift time, defined by the first S2. 
    s1_before_s2 = np.zeros(n_events, dtype=bool)

    waveform_area = np.zeros(n_events)  #integrate the whole waveform
    n_wfms_summed = 0
    avg_wfm = np.zeros(wsize)

    # Temporary, for testing low area, multiple-S1 events
    dt = np.zeros(n_events)
    small_weird_areas = np.zeros(n_events)
    big_weird_areas = np.zeros(n_events)
    


    
    n_golden = 0
    inn=""

    inn="" # used to control hand scan

    empty_evt_ind = np.zeros(n_events)


    j = 0
    for compressed_file in compressed_file_list:
    
        # load compressed data
        try:
            with np.load(compressed_file) as data:
                ch_data = data["arr_0"]
        except:
            print("Error in loading "+compressed_file)
            continue
        
        n_tot_samp_per_ch = int( (ch_data.size)/n_sipms )
        n_events_b = int((ch_data.size)/(n_sipms*wsize))
    
        ch_data = np.concatenate((ch_data.astype(int), np.zeros(n_tot_samp_per_ch, dtype="int")), dtype="int")
        ch_data = np.reshape(ch_data, (n_channels,n_events_b,wsize))

        v_matrix_all_ch = ch_data*vscale
    
    
        v_bls_matrix_all_ch = np.zeros_like(v_matrix_all_ch)
        for ch in range(n_channels-1):
            v_bls_matrix_all_ch[ch,:,:] = v_matrix_all_ch[ch,:,:] - np.mean(v_matrix_all_ch[ch,:,0:100])
            v_bls_matrix_all_ch[ch,:,:] = v_matrix_all_ch[ch,:,:]*tscale*(1000)/spe_sizes[ch] # scaled by tscale and spe size
    
        v_bls_matrix_all_ch[-1,:,:] = np.sum(v_bls_matrix_all_ch[:,:,:], axis=0)
    
    
    
        # create a time axis in units of µs:
        x = np.arange(0, wsize, 1)
        t = tscale*x
        #t_matrix = np.repeat(t[np.newaxis,:], V.size/wsize, 0)
        # Note: if max_evts != -1, we won't load in all events in the dataset
        n_events = int(v_matrix_all_ch[0].shape[0])
        if n_events == 0: break
            
        # perform baseline subtraction for full event (separate from pulse-level baseline subtraction):
        baseline_start = int(0./tscale)
        baseline_end = np.min((int(wsize*0.2), 1000))

        # baseline subtracted (bls) waveforms saved in this matrix:
        #v_bls_matrix_all_ch = np.zeros( np.shape(v_matrix_all_ch), dtype=array_dtype) # dims are (chan #, evt #, sample #)

        #t_end_wfm_fill = time.time()
        #print("Time to fill all waveform arrays: ", t_end_wfm_fill - t_end_load)

        #print("Events to process: ",n_events)
        for i in range(0, n_events):
            
            sum_baseline = np.mean( v_matrix_all_ch[-1][i,baseline_start:baseline_end] ) #avg ~us, avoiding trigger
            baselines = [ np.mean( ch_j[i,baseline_start:baseline_end] ) for ch_j in v_matrix_all_ch ]
            
            sum_data = v_matrix_all_ch[-1][i,:] - sum_baseline
            ch_data = [ch_j[i,:]-baseline_j for (ch_j,baseline_j) in zip(v_matrix_all_ch,baselines)]
            
            #v_bls_matrix_all_ch[:,i,:] = ch_data
            
        #v_bls_matrix_all_ch = v_matrix_all_ch
        


        # ==================================================================
        # ==================================================================
        # now setup for pulse finding on the baseline-subtracted sum waveform

        
    #check mark

        #print("Running pulse finder on {:d} events...".format(n_events))

        # use for coloring pulses
        pulse_class_colors = np.array(['blue', 'green', 'red', 'magenta', 'darkorange'])
        pulse_class_labels = np.array(['Other', 'S1-like LXe', 'S1-like gas', 'S2-like', 'Merged S1/S2'])
        pc_legend_handles=[]
        for class_ind in range(len(pulse_class_labels)):
            pc_legend_handles.append(mpl.patches.Patch(color=pulse_class_colors[class_ind], label=str(class_ind)+": "+pulse_class_labels[class_ind]))

        for i in range(j*block_size, j*block_size+n_events):
            if (i)%2000==0: print("Event #",i)
            
            # Find pulse locations; other quantities for pf tuning/debugging
            start_times, end_times, peaks, data_conv, properties = pf.findPulses( v_bls_matrix_all_ch[-1,i-j*block_size,:], max_pulses , SPEMode=False)

            # Sort pulses by start times, not areas
            startinds = np.argsort(start_times)
            n_pulses[i] = min(max_pulses,len(start_times))
            if (n_pulses[i] < 1):
                #print("No pulses found for event {0}; skipping".format(i))
                empty_evt_ind[i] = i
                #continue
            for m in startinds:
                if m >= max_pulses:
                    continue
                p_start[i,m] = start_times[m]
                p_end[i,m] = end_times[m]

            # Individual channel pulse locations, in case you want this info
            # Can't just ":" the the first index in data, findPulses doesn't like it, so have to loop 
            #for j in range(n_channels-1):
            #    start_times_ch, end_times_ch, peaks_ch, data_conv_ch, properties_ch = pf.findPulses( v_bls_matrix_all_ch[j,i,:], max_pulses )
                # Sorting by start times from the sum of channels, not each individual channel
            #    for k in startinds:
            #        if k >= len(start_times_ch):
            #            continue
            #        p_start_ch[i,k,j] = start_times_ch[k]
            #        p_end_ch[i,k,j] = end_times_ch[k]
                

            # More precisely estimate baselines immediately before each pulse
            baselines_precise = pq.GetBaselines(p_start[i,:n_pulses[i]], p_end[i,:n_pulses[i]], v_bls_matrix_all_ch[:,i-j*block_size,:])

            # Calculate interesting quantities, only for pulses that were found
            for pp in range(n_pulses[i]):
                # subtract out more precise estimate of baseline for better RQ estimates
                baselines_pulse = baselines_precise[pp] # array of baselines per channel, for this pulse
                v_pulse_bls = np.array([ch_j - baseline_j for (ch_j, baseline_j) in zip(v_bls_matrix_all_ch[:,i-j*block_size,:], baselines_pulse)])

                # Version w/o pulse-level baseline subtraction
                #v_pulse_bls = v_bls_matrix_all_ch[:,i-j*block_size,:]

                # copied from above, for reference
                #sum_data = v_matrix_all_ch[-1][i, :] - sum_baseline
                #ch_data = [ch_j[i, :] - baseline_j for (ch_j, baseline_j) in zip(v_matrix_all_ch, baselines)]

                # Area, max & min heights, width, pulse mean & rms
                p_area[i,pp] = pq.GetPulseArea(p_start[i,pp], p_end[i,pp], v_pulse_bls[-1] )
                p_max_height[i,pp] = pq.GetPulseMaxHeight(p_start[i,pp], p_end[i,pp], v_pulse_bls[-1] )
                p_min_height[i,pp] = pq.GetPulseMinHeight(p_start[i,pp], p_end[i,pp], v_pulse_bls[-1] )
                p_width[i,pp] = p_end[i,pp] - p_start[i,pp]
                #(p_mean_time[i,pp], p_rms_time[i,pp]) = pq.GetPulseMeanAndRMS(p_start[i,pp], p_end[i,pp], v_bls_matrix_all_ch[-1,i,:])

                # Area and height fractions      
                (p_afs_2l[i,pp], p_afs_1[i,pp], p_afs_25[i,pp], p_afs_50[i,pp], p_afs_75[i,pp], p_afs_99[i,pp]) = pq.GetAreaFraction(p_start[i,pp], p_end[i,pp], v_pulse_bls[-1] )
                try: (p_hfs_10l[i,pp], p_hfs_50l[i,pp], p_hfs_10r[i,pp], p_hfs_50r[i,pp]) = pq.GetHeightFractionSamples(p_start[i,pp], p_end[i,pp], v_pulse_bls[-1] )
                except: 1==1
            
                # Areas for individual channels and top bottom
                p_area_ch[i,pp,:] = pq.GetPulseAreaChannel(p_start[i,pp], p_end[i,pp], v_pulse_bls )
                p_area_ch_frac[i,pp,:] = p_area_ch[i,pp,:]/p_area[i,pp]
                p_area_top[i,pp] = sum(p_area_ch[i,pp,top_channels])
                p_area_bottom[i,pp] = sum(p_area_ch[i,pp,bottom_channels])
                p_tba[i, pp] = (p_area_top[i, pp] - p_area_bottom[i, pp]) / (p_area_top[i, pp] + p_area_bottom[i, pp])
                

                # Centroids

                board_offset = 0 # Gap between SiPM quadrants (0 = assuming flush)
                l = 0.59 # SiPM width/length (not exactly a square...see specs)
                d1 = 0.75 + board_offset # distance from center of board to quadrant center 
                d2 = 0.025 # distance from quadrant center to near SiPM edge
                d3 = d2 + l # distance from quadrant center to far SiPM edge
                d4 = 0.32 # distance from quadrant center to SiPM center
                r_tpc = (1.175/2)*2.54 # TPC (inner) radius

                w1 = (d1-d4) # Weight for near SiPMs
                w2 = (d1+d4) # Weight for far SiPMs

                fudge = 1


                # X
                # +1: 2,3,14,15
                # +3: 1,4,13,16
                # -1: 5,8,9,12
                # -3: 6,7,10,11
                
                # Y
                # +1: 7,8,3,4
                # +3: 6,5,2,1
                # -1: 10,9,14,13
                # -3: 11,12,15,16
                
                b0 = 15 - 1
                
                center_bot_x[i,pp] += w1*(p_area_ch[i,pp,b0+2]+p_area_ch[i,pp,b0+3]+p_area_ch[i,pp,b0+14]+p_area_ch[i,pp,b0+15])
                center_bot_x[i,pp] += w2*(p_area_ch[i,pp,b0+1]+p_area_ch[i,pp,b0+4]+p_area_ch[i,pp,b0+13]+p_area_ch[i,pp,b0+16])
                center_bot_x[i,pp] += -w1*(p_area_ch[i,pp,b0+5]+p_area_ch[i,pp,b0+8]+p_area_ch[i,pp,b0+9]+p_area_ch[i,pp,b0+12])
                center_bot_x[i,pp] += -w2*(p_area_ch[i,pp,b0+6]+p_area_ch[i,pp,b0+7]+p_area_ch[i,pp,b0+10]+p_area_ch[i,pp,b0+11])
                center_bot_x[i,pp] *= (fudge/p_area_bottom[i,pp])
                
                center_bot_y[i,pp] += w1*(p_area_ch[i,pp,b0+7]+p_area_ch[i,pp,b0+8]+p_area_ch[i,pp,b0+3]+p_area_ch[i,pp,b0+4])
                center_bot_y[i,pp] += w2*(p_area_ch[i,pp,b0+6]+p_area_ch[i,pp,b0+5]+p_area_ch[i,pp,b0+2]+p_area_ch[i,pp,b0+1])
                center_bot_y[i,pp] += -w1*(p_area_ch[i,pp,b0+10]+p_area_ch[i,pp,b0+9]+p_area_ch[i,pp,b0+14]+p_area_ch[i,pp,b0+13])
                center_bot_y[i,pp] += -w2*(p_area_ch[i,pp,b0+11]+p_area_ch[i,pp,b0+12]+p_area_ch[i,pp,b0+15]+p_area_ch[i,pp,b0+16])
                center_bot_y[i,pp] *= (fudge/p_area_bottom[i,pp])
                
                
                t0 = -1
                
                # X
                # +1: 5,8,9,12
                # +3: 6,7,10,11
                # -1: 2,3,14,15
                # -3: 1,4,13,16
                
                # Y
                # +1: 4,3,8,7
                # +3: 1,2,5,6
                # -1: 13,14,9,10
                # -3: 11,12,15,16
                
                center_top_x[i,pp] += w1*(p_area_ch[i,pp,t0+2]+p_area_ch[i,pp,t0+3]+p_area_ch[i,pp,t0+14]+p_area_ch[i,pp,t0+15])
                center_top_x[i,pp] += w2*(p_area_ch[i,pp,t0+1]+p_area_ch[i,pp,t0+4]+p_area_ch[i,pp,t0+13]+p_area_ch[i,pp,t0+16])
                center_top_x[i,pp] += -w1*(p_area_ch[i,pp,t0+5]+p_area_ch[i,pp,t0+8]+p_area_ch[i,pp,t0+9]+p_area_ch[i,pp,t0+12])
                center_top_x[i,pp] += -w2*(p_area_ch[i,pp,t0+6]+p_area_ch[i,pp,t0+7]+p_area_ch[i,pp,t0+10]+p_area_ch[i,pp,t0+11])
                center_top_x[i,pp] *= (-fudge/p_area_top[i,pp])
                
                center_top_y[i,pp] += w1*(p_area_ch[i,pp,t0+7]+p_area_ch[i,pp,t0+8]+p_area_ch[i,pp,t0+3]+p_area_ch[i,pp,t0+4])
                center_top_y[i,pp] += w2*(p_area_ch[i,pp,t0+6]+p_area_ch[i,pp,t0+5]+p_area_ch[i,pp,t0+2]+p_area_ch[i,pp,t0+1])
                center_top_y[i,pp] += -w1*(p_area_ch[i,pp,t0+10]+p_area_ch[i,pp,t0+9]+p_area_ch[i,pp,t0+14]+p_area_ch[i,pp,t0+13])
                center_top_y[i,pp] += -w2*(p_area_ch[i,pp,t0+11]+p_area_ch[i,pp,t0+12]+p_area_ch[i,pp,t0+15]+p_area_ch[i,pp,t0+16])
                center_top_y[i,pp] *= (fudge/p_area_top[i,pp])
                

            waveform_area[i] = np.sum(v_pulse_bls[-1])

                
                
                
                
                
                

                
            # Pulse classifier, work in progress
            p_class[i,:] = pc.ClassifyPulses(p_tba[i, :], (p_afs_50[i, :]-p_afs_2l[i, :])*tscale, n_pulses[i], p_area[i,:])

            # Event level analysis. Look at events with both S1 and S2.
            index_s1 = (p_class[i,:] == 1) + (p_class[i,:] == 2) # S1's
            index_s2 = (p_class[i,:] == 3) + (p_class[i,:] == 4) # S2's
            n_s1[i] = np.sum(index_s1)
            n_s2[i] = np.sum(index_s2)
            
            if n_s1[i] > 0:
                sum_s1_area[i] = np.sum(p_area[i, index_s1])
            if n_s2[i] > 0:
                sum_s2_area[i] = np.sum(p_area[i, index_s2])
            if n_s1[i] == 1:
                if n_s2[i] == 1:
                    drift_Time[i] = tscale*(p_start[i, np.argmax(index_s2)] - p_start[i, np.argmax(index_s1)])
                    drift_Time_AF[i] = tscale*(p_afs_1[i, np.argmax(index_s2)] - p_afs_1[i, np.argmax(index_s1)])
                    drift_Time_AF50[i] = tscale*(p_afs_50[i, np.argmax(index_s2)] - p_afs_50[i, np.argmax(index_s1)])
                    drift_Time_AF150[i] = tscale*(p_afs_50[i, np.argmax(index_s2)] - p_afs_1[i, np.argmax(index_s1)])
                    drift_Time_AS[i] = tscale*(p_start[i, np.argmax(index_s2)] - p_start[i, np.argmax(index_s1)])
                if n_s2[i] > 1:
                    s1_before_s2[i] = np.argmax(index_s1) < np.argmax(index_s2) 
                    drift_Time_AS[i] = tscale*(p_start[i, np.argmax(index_s2)] - p_start[i, np.argmax(index_s1)]) #For multi-scatter events. 
            
            if drift_Time[i]>0:
                n_golden += 1


            # =============================================================
            # draw the waveform and the pulse bounds found

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

            afs50_2 = (p_afs_50[i,:]-p_afs_2l[i,:])*tscale
            
            # Condition to skip the individual plotting, hand scan condition
            if np.size(handscan) == 1: 
                plotyn = handscan
            else:
                plotyn = handscan[i]
            #plotyn = drift_Time[i]<2 and drift_Time[i]>0 and np.any((p_tba[i,:]>-0.75)*(p_tba[i,:]<-0.25)*(p_area[i,:]<3000)*(p_area[i,:]>1400))#np.any((p_tba[i,:]>-0.91)*(p_tba[i,:]<-0.82)*(p_area[i,:]<2800)*(p_area[i,:]>1000))# True#np.any(p_class[i,:]==4)#False#np.any(p_area[i,:]>1000) and 
            #plotyn = drift_Time[i]>2.5 and (center_bot_y[i,0]**2+center_bot_x[i,0]**2) <0.1
            #plotyn = True #np.any((p_class[i,:] == 3) + (p_class[i,:] == 4))#np.any((p_tba[i,:]>-0.75)*(p_tba[i,:]<-0.25)*(p_area[i,:]<3000)*(p_area[i,:]>1000))
            #plotyn = np.any((p_tba[i,:]>-1)*(p_tba[i,:]<-0.25)*(p_area[i,:]<30000)*(p_area[i,:]>3000))#np.any((np.log10(p_area[i,:])>3.2)*(np.log10(p_area[i,:])<3.4) )#False#
            # Pulse area condition
            # afs50_2 = (p_afs_50[i,:]-p_afs_2l[i,:])*tscale
            # lower_limit = -0.13*(np.log10(p_area[i,:])-3.2)**2-1.25
            # higher_limit = -0.2*(np.log10(p_area[i,:])-4)**2-0.4
            # afs50_2 = (p_afs_50[i,:]-p_afs_2l[i,:])*tscale
            # temp_condition = (np.log10(afs50_2)>-0.75)*(np.log10(afs50_2)<-0.6)*(np.log10(p_area[i,:])>3.2)*(np.log10(p_area[i,:])<4.4)
            # plotyn = np.any(temp_condition)
            plotyn = False
            
            areaRange = np.sum((p_area[i,:] < 50)*(p_area[i,:] > 5))
            if areaRange > 0:
                dt[i] = abs(p_start[i,1] - p_start[i,0]) # For weird double s1 data
                weird_areas =[p_area[i,0], p_area[i,1] ]
                small_weird_areas[i] = min(weird_areas)
                big_weird_areas[i] = max(weird_areas)

            # Condition to include a wfm in the average
            add_wfm = np.any((p_area[i,:]>5000)*(p_tba[i,:]<-0.75))*(n_s1[i]==1)*(n_s2[i]==0)
            if add_wfm and save_avg_wfm:
                plotyn = add_wfm # in avg wfm mode, plot the events which will go into the average
                avg_wfm += v_bls_matrix_all_ch[-1,i-j*block_size,:]
                n_wfms_summed += 1

            # Both S1 and S2 condition
            s1s2 = (n_s1[i] == 1)*(n_s2[i] == 1)

            if inn == 's': sys.exit()
            
            if not inn == 'q' and plotyn: # and drift_Time[i] > 0 and s1s2 and p_area[i,0]<500. and p_area[i,0]>20.: # and np.any((p_area[i,:] > 350)*(p_area[i,:] < 600)): # and drift_Time[i] > 0: # plotyn: #plot_event_ind == i and plotyn:


                fig = pl.figure()
                ax = pl.gca()

                #pl.plot(x*tscale, v_bls_matrix_all_ch[23,i-j*block_size,:]/(tscale*(1000)/spe_sizes[23]) , "black")
                pl.plot(x*tscale, np.sum(v_bls_matrix_all_ch[0:15,i-j*block_size,:],axis=0), color="red",lw=0.7)
                pl.plot(x*tscale, np.sum(v_bls_matrix_all_ch[16:31,i-j*block_size,:],axis=0), color="green",lw=0.7)
                pl.plot(x*tscale, v_bls_matrix_all_ch[-1,i-j*block_size,:],color='black',lw=0.7 )
                #for ch in range(32):
                #    pl.plot(x*tscale, v_bls_matrix_all_ch[ch,i-j*block_size,:])
                #pl.ylabel("mV")
                #pl.plot(x*tscale, v_bls_matrix_all_ch[-1,i-j*block_size,:]/np.max(v_bls_matrix_all_ch[-1,i-j*block_size,:]),'black' )
                #pl.plot(x*tscale, np.sum(v_bls_matrix_all_ch[0:15,i-j*block_size,:],axis=0)/np.max(np.sum(v_bls_matrix_all_ch[0:15,i-j*block_size,:],axis=0)), "red")
                #pl.plot(x*tscale, np.sum(v_bls_matrix_all_ch[16:31,i-j*block_size,:],axis=0)/np.max(np.sum(v_bls_matrix_all_ch[16:31,i-j*block_size,:],axis=0)), "green")
                pl.xlabel(r"Time [$\mu$s]")
                pl.ylabel("phd/sample")
                pl.title("Event {}".format(i))
                for ps in range(n_pulses[i]):
                    #pl.axvspan(tscale*start_times[ps],tscale*end_times[ps],alpha=0.20,color="b")
                    pl.axvspan(tscale*start_times[ps],tscale*end_times[ps],alpha=0.3,color=pulse_class_colors[p_class[i,ps]],zorder=0)
                    pl.axvline(tscale*p_afs_1[i,ps],color=pulse_class_colors[p_class[i,ps]],zorder=0,linestyle='--')
                    
                    ax.text((end_times[ps]) * tscale, (0.94-ps*0.2) * ax.get_ylim()[1], '{:.1f} phd'.format(p_area[i, ps]),
                            fontsize=9, color=pulse_class_colors[p_class[i, ps]])
                    ax.text((end_times[ps]) * tscale, (0.9-ps*0.2) * ax.get_ylim()[1], 'TBA={:.1f}'.format(p_tba[i, ps]),
                        fontsize=9, color=pulse_class_colors[p_class[i, ps]])
                    ax.text((end_times[ps]) * tscale, (0.86-ps*0.2) * ax.get_ylim()[1], 'Rise={:.1f} us'.format(afs50_2[ps]),
                        fontsize=9, color=pulse_class_colors[p_class[i, ps]])
                    #ax.text((end_times[ps]) * tscale, (0.82-ps*0.2) * ax.get_ylim()[1], 'Check={}'.format(temp_condition[ps]),
                        #fontsize=9, color=pulse_class_colors[p_class[i, ps]])
                    
                pl.legend(["All","Summed Top","Summed Bottom"])
                
                #for ch in range(32):
                    #if ch < 16:
                    #pl.plot( x*tscale, v_bls_matrix_all_ch[ch,i-j*block_size,:], "r")
                    #if ch > 15 and ch < 23
                #pl.plot( x*tscale, v_bls_matrix_all_ch[20,i-j*block_size,:], "m" )
                #pl.plot( x*tscale, v_bls_matrix_all_ch[30,i-j*block_size,:], "black" )
                pl.grid(which="both",axis="both",linestyle="--")
                pl.xlim(0,event_window)
                #pl.legend(("Summed waveform"))
                pl.show()

                """
                fig = pl.figure(1,figsize=(10, 7))
                pl.rc('xtick', labelsize=10)
                pl.rc('ytick', labelsize=10)
                
                ax = pl.subplot2grid((2,2),(0,0))
                #pl.title("Top array, event "+str(i))
                #pl.grid(b=True,which='major',color='lightgray',linestyle='--')
                #ch_labels = [int(i) for i in range(n_channels)]
                ch_colors = [pl.cm.tab10(ii) for ii in range(n_channels)]
                for pulse in range(len(start_times)): # fill found pulse regions for top
                    ax.axvspan(start_times[pulse] * tscale, end_times[pulse] * tscale, alpha=0.25,
                               color=pulse_class_colors[p_class[i, pulse]])
                
                
                for i_chan in range(n_channels-1):
                    if i_chan == (n_channels-1)/2:
                        ax = pl.subplot2grid((2,2),(0,1))
                        pl.title("Bottom array, event "+str(i))
                        pl.grid(b=True,which='major',color='lightgray',linestyle='--')
                        for pulse in range(len(start_times)):  # fill found pulse regions for bottom
                            ax.axvspan(start_times[pulse] * tscale, end_times[pulse] * tscale, alpha=0.25,
                                       color=pulse_class_colors[p_class[i, pulse]])
                    
                    pl.plot(t,v_bls_matrix_all_ch[i_chan,i-j*block_size,:],color=ch_colors[i_chan],label=ch_labels[i_chan])
                    #pl.plot( x, v_bls_matrix_all_ch[i_chan,i,:],color=ch_colors[i_chan],label=ch_labels[i_chan] )
                    pl.xlim([trigger_time_us-8,trigger_time_us+8])
                    #pl.xlim([wsize/2-4000,wsize/2+4000])
                    pl.ylim([-5, 3000/chA_spe_size])
                    pl.xlabel('Time (us)')
                    #pl.xlabel('Samples')
                    pl.ylabel('phd/sample')
                    pl.legend()
                
                
                #ax = pl.subplot2grid((2,2),(1,0),colspan=2)
                #pl.plot(t,v_bls_matrix_all_ch[-1,i,:],'blue')
                pl.plot( x*tscale, v_bls_matrix_all_ch[-1,i-j*block_size,:],'blue' )
                #pl.xlim([0,wsize])
                pl.xlim([0,event_window])
                pl.ylim( [-1, 1.01*np.max(v_bls_matrix_all_ch[-1,i-j*block_size,:])])
                pl.xlabel('Time (us)')
                #pl.xlabel('Samples')
                pl.ylabel('phd/sample')
                pl.title("Sum, event "+ str(i))
                pl.grid(b=True,which='major',color='lightgray',linestyle='--')
                pl.legend(handles=pc_legend_handles)

                for pulse in range(len(start_times)):
                    ax.axvspan(start_times[pulse] * tscale, end_times[pulse] * tscale, alpha=0.25, color=pulse_class_colors[p_class[i, pulse]])
                    ax.text((end_times[pulse]) * tscale, 0.9 * ax.get_ylim()[1], '{:.1f} phd'.format(p_area[i, pulse]),
                            fontsize=9, color=pulse_class_colors[p_class[i, pulse]])
                    ax.text((end_times[pulse]) * tscale, 0.8 * ax.get_ylim()[1], 'TBA={:.1f}'.format(p_tba[i, pulse]),
                        fontsize=9, color=pulse_class_colors[p_class[i, pulse]])
                
                #ax.axhline( 0.276, 0, wsize, linestyle='--', lw=1, color='orange')

                # Debugging of pulse finder
                debug_pf = False #True
                if debug_pf and n_pulses[i]>0:
                    #pl.plot(t_matrix[i-j*block_size, :], data_conv, 'red')
                    #pl.plot(t_matrix[i-j*block_size, :], np.tile(0., np.size(data_conv)), 'gray')
                    pl.vlines(x=peaks*tscale, ymin=data_conv[peaks] - properties["prominences"],
                               ymax=data_conv[peaks], color="C1")
                    pl.hlines(y=properties["width_heights"], xmin=properties["left_ips"]*tscale,
                               xmax=properties["right_ips"]*tscale, color="C1")
                    #print("pulse heights: ", data_conv[peaks] )
                    #print("prominences:", properties["prominences"])

                #pl.draw()
                pl.show(block=0)
                """
                inn = input("Press enter to continue, q to stop plotting, evt # to skip to # (forward only)")
                #fig.clf()
                
        # end of pulse finding and plotting event loop

        if save_avg_wfm:
            avg_wfm /= n_wfms_summed
            np.savetxt(data_dir+'average_waveform.txt',avg_wfm)
            print("Average waveform saved")

        n_events = i
        t_end = time.time()
        print("total number of events processed:", n_events)
        #print("Time used: {}".format(t_end-t_start))

        print("empty events: {0}".format(np.sum(empty_evt_ind>0)))
        # pl.hist(empty_evt_ind[empty_evt_ind>0], bins=1000)
        # pl.xlabel('Empty event index')
        # pl.show()

        j += 1

    #create a dictionary with all RQs
    list_rq = {}
    list_rq['center_top_x'] = center_top_x
    list_rq['center_top_y'] = center_top_y
    list_rq['center_bot_x'] = center_bot_x
    list_rq['center_bot_y'] = center_bot_y
    list_rq['n_s1'] = n_s1
    list_rq['n_s2'] = n_s2
    list_rq['s1_before_s2'] = s1_before_s2
    list_rq['n_pulses'] = n_pulses
    list_rq['n_events'] = n_events
    list_rq['p_area'] = p_area
    list_rq['p_class'] = p_class
    list_rq['drift_Time'] = drift_Time
    list_rq['drift_Time_AF'] = drift_Time_AF
    list_rq['drift_Time_AF50'] = drift_Time_AF50
    list_rq['drift_Time_AF150'] = drift_Time_AF150
    list_rq['drift_Time_AS'] = drift_Time_AS
    list_rq['p_max_height'] = p_max_height
    list_rq['p_min_height'] = p_min_height
    list_rq['p_width'] = p_width
    list_rq['p_afs_1'] = p_afs_1
    list_rq['p_afs_99'] = p_afs_99
    list_rq['p_afs_2l'] = p_afs_2l
    list_rq['p_afs_50'] = p_afs_50
    list_rq['p_area_ch'] = p_area_ch
    list_rq['p_area_ch_frac'] = p_area_ch_frac
    list_rq['p_area_top'] = p_area_top
    list_rq['p_area_bottom'] = p_area_bottom
    list_rq['p_tba'] = p_tba
    list_rq['p_start'] = p_start
    list_rq['p_end'] = p_end
    list_rq['sum_s1_area'] = sum_s1_area
    list_rq['sum_s2_area'] = sum_s2_area
    list_rq['waveform_area'] = waveform_area
    #list_rq[''] =    #add more rq

    #remove zeros in the end of each RQ array. 
    for rq in list_rq.keys():
        if rq != 'n_events':
            list_rq[rq] = list_rq[rq][:n_events]

    rq = open(data_dir + "rq.npz",'wb')
    #rq = open("/home/xaber/Data/20220303/202203031433 /rq.npz")
    
    np.savez(rq, **list_rq)
    rq.close()

  


def main():
    with open("path.txt", 'r') as path:
        #data_dir = "/home/xaber/Data/data-202203/20220323/202203232045_1.4bar_2600C2400G0A_54B_topCo_15us/"
        data_dir = path.read()
        #data_dir = data_dir[:-1]
    
    make_rq(data_dir)

if __name__ == "__main__":
    main()

