import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mpl
import time
import sys
import glob
from natsort import natsorted

import PulseFinderScipy as pf
import PulseQuantities as pq
import PulseClassification as pc
import PulseFinderVerySimple as vs
from read_settings import get_event_window, get_vscale
#from ch_evt_filter_compress import filter_channel_event


def make_rq(data_dir, handscan=False, max_pulses=4, filtered=True, simpleS2=True, save_avg_wfm=False, phase="liquid"):
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

    # Load headers and calculate event time
    h_file = np.load(data_dir+"/compressed_filtered_data/headers.npz")
    h_array = h_file["arr_0"]
    h_n_events = int(np.floor(h_array.size/8))
    h_array = np.reshape(h_array,(h_n_events,8))

    # Calculate event time
    # Precision up to 0.5 ms. To-do: get precision to 16 ns
    second_16 = h_array[:,5]
    second_16[second_16 < 0] = second_16[second_16 < 0] + 2**15
    second_16_next = np.zeros(h_n_events,dtype=int)
    for i in range(1,h_n_events):
        if second_16[i] - second_16[i-1] < 0:
            second_16_next[i:] += 1
    ev_time_s = 16*(second_16 + second_16_next*2**15)*(10**-9 * 2**15)

    if h_n_events < n_events:
        # Match header size to n_events which has some extra zeros
        ev_time_s = np.concatenate( (ev_time_s,np.zeros(n_events-h_n_events,dtype=int) ) )
    elif h_n_events > n_events:
        # This shouldn't happen
        print("Warning: Header file contains more events than compressed events!")
        ev_time_s = np.concatenate( (ev_time_s,np.zeros(h_n_events-n_events,dtype=int) ) )
    

    # ====================================================================================================================================
    # Initialize rq's to save
    # np.zeros is preferred over np.empty bc/ we want zero to be default value

    p_start = np.zeros(( n_events, max_pulses), dtype=int)
    p_end   = np.zeros(( n_events, max_pulses), dtype=int)
    p_found = np.zeros(( n_events, max_pulses), dtype=int)
    p_area = np.zeros(( n_events, max_pulses))
    p_max_height = np.zeros(( n_events, max_pulses))
    p_min_height = np.zeros(( n_events, max_pulses))
    p_width = np.zeros(( n_events, max_pulses))
    p_mean_time = np.zeros((n_events, max_pulses) )
    p_rms_time = np.zeros((n_events, max_pulses) )
    p_area_top = np.zeros((n_events, max_pulses))
    p_area_bottom = np.zeros((n_events, max_pulses))
    p_tba = np.zeros((n_events, max_pulses))
    p_class = np.zeros((n_events, max_pulses), dtype=int)

    # centroid 
    center_top_x = np.zeros(( n_events, max_pulses))
    center_top_y = np.zeros(( n_events, max_pulses))
    center_bot_x = np.zeros(( n_events, max_pulses))
    center_bot_y = np.zeros(( n_events, max_pulses))

    # AFT's 
    p_afs_2l = np.zeros((n_events, max_pulses) )
    p_afs_2r = np.zeros((n_events, max_pulses) )
    p_afs_1 = np.zeros((n_events, max_pulses) )
    p_afs_10 = np.zeros((n_events, max_pulses) )
    p_afs_25 = np.zeros((n_events, max_pulses) )
    p_afs_50 = np.zeros((n_events, max_pulses) )
    p_afs_75 = np.zeros((n_events, max_pulses) )
    p_afs_90 = np.zeros((n_events, max_pulses) )
    p_afs_99 = np.zeros((n_events, max_pulses) )

    # HFT's     
    p_hfs_10l = np.zeros((n_events, max_pulses) )
    p_hfs_50l = np.zeros((n_events, max_pulses) )
    p_hfs_90l = np.zeros((n_events, max_pulses) )
    p_hfs_10r = np.zeros((n_events, max_pulses) )
    p_hfs_50r = np.zeros((n_events, max_pulses) )
    p_hfs_90r = np.zeros((n_events, max_pulses) )

    # Channel level (per event, per pulse, per channel)
    p_start_ch = np.zeros((n_events, max_pulses, n_channels-1), dtype=int)
    p_end_ch = np.zeros((n_events, max_pulses, n_channels-1), dtype=int )
    p_area_ch = np.zeros((n_events, max_pulses, n_channels-1) )
    p_area_ch_frac = np.zeros((n_events, max_pulses, n_channels-1) )
    p_max_height_ch = np.zeros((n_events, max_pulses, n_channels-1) )
    #p_coincidence = np.zeros((n_events, max_pulses) ) # technically pulse level
    
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
    drift_Time_max = np.zeros(n_events)  # drift time defined by the interval between the biggest S1 and biggest S2 within an event window
    index_max_s1 = np.zeros(n_events, dtype = "int")-1
    index_max_s2 = np.zeros(n_events, dtype = "int")-1    # set default to -1
    drift_Time_AS = np.zeros(n_events) # for multi-scatter drift time, defined by the first S2. 
    s1_before_s2 = np.zeros(n_events, dtype=bool)
    right_area = np.zeros(n_events)

    # Avg waveform quantities
    waveform_area = np.zeros(n_events)  #integrate the whole waveform
    rq_ev_time_s = np.zeros(n_events,dtype="int")
    n_wfms_summed = 0
    avg_wfm = np.zeros(wsize)

    empty_evt_ind = np.zeros(n_events) # empty event quantity


    # ====================================================================================================================================
    # Begin loop over compressed files

    j = 0 # index for number of compressed files
    counter = 0 # index for total events

    for compressed_file in compressed_file_list:
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

            # Event time
            rq_ev_time_s[i] = ev_time_s[counter]
            
            # New pulse finder
            start_times = []
            end_times = []
            lh_cut = wsize
            for g in range(max_pulses):
                if lh_cut < 1: continue
                temp_start, temp_end = vs.PulseFinderVerySimple(ch_data_phdPerSample[-1,i-j*block_size,:lh_cut], verbose=False)
                if temp_start != temp_end:
                    if g==0: right_area[i] = np.sum(ch_data_phdPerSample[-1,i-j*block_size,temp_end:])
                    start_times.append(temp_start)
                    end_times.append(temp_end)
                    lh_cut = temp_start
                else:
                    continue

            # Data is already filtered
            data_conv = ch_data_phdPerSample[-1,i-j*block_size,:]


            # Sort pulses by start times, not areas as given by pf.findPulses
            startinds = np.argsort(start_times)
            n_pulses[i] = min(max_pulses,len(start_times))
            if (n_pulses[i] < 1):
                empty_evt_ind[i] = i
            mp = 0
            for m in startinds:
                if m >= max_pulses:
                    continue
                if start_times[m] < 0.25/tscale: continue
                p_start[i,mp] = start_times[m]
                p_end[i,mp] = end_times[m]
                mp += 1
                

            # More precisely estimate baselines immediately before each pulse
            baselines_precise = pq.GetBaselines(p_start[i,:n_pulses[i]], p_end[i,:n_pulses[i]], ch_data_phdPerSample[:,i-j*block_size,:])


            # ====================================================================================================================================
            # Loop over pulses, calculate some pulse level rq's

            for pp in range(n_pulses[i]):
                # subtract out more precise estimate of baseline for better RQ estimates
                baselines_pulse = baselines_precise[pp] # array of baselines per channel, for this pulse
                ch_data_sum_pulse_bls = np.array([ch_j - baseline_j for (ch_j, baseline_j) in zip(ch_data_phdPerSample[:,i-j*block_size,:], baselines_pulse)])

                # Area, max & min heights, width, pulse mean & rms
                p_area[i,pp] = pq.GetPulseArea(p_start[i,pp], p_end[i,pp], ch_data_sum_pulse_bls[-1] )
                p_max_height[i,pp] = pq.GetPulseMaxHeight(p_start[i,pp], p_end[i,pp], ch_data_sum_pulse_bls[-1] )
                p_min_height[i,pp] = pq.GetPulseMinHeight(p_start[i,pp], p_end[i,pp], ch_data_sum_pulse_bls[-1] )
                p_width[i,pp] = p_end[i,pp] - p_start[i,pp]

                # Area and height fractions      
                (p_afs_2l[i,pp], p_afs_1[i,pp], p_afs_10[i,pp],p_afs_25[i,pp], p_afs_50[i,pp], p_afs_75[i,pp], p_afs_90[i,pp], p_afs_99[i,pp]) = pq.GetAreaFraction(p_start[i,pp], p_end[i,pp], ch_data_sum_pulse_bls[-1] )
                (p_hfs_10l[i,pp], p_hfs_50l[i,pp], p_hfs_90l[i,pp], p_hfs_10r[i,pp], p_hfs_50r[i,pp], p_hfs_90r[i,pp]) = pq.GetHeightFractionSamples(p_start[i,pp], p_end[i,pp], ch_data_sum_pulse_bls[-1] )
                
                # Areas for individual channels and top bottom
                p_area_ch[i,pp,:] = pq.GetPulseAreaChannel(p_start[i,pp], p_end[i,pp], ch_data_sum_pulse_bls[:-1] )
                p_area_ch_frac[i,pp,:] = p_area_ch[i,pp,:]/p_area[i,pp]
                p_area_top[i,pp] = np.sum(p_area_ch[i,pp,top_channels])
                p_area_bottom[i,pp] = np.sum(p_area_ch[i,pp,bottom_channels])
                p_tba[i, pp] = (p_area_top[i, pp] - p_area_bottom[i, pp]) / (p_area_top[i, pp] + p_area_bottom[i, pp])

                p_max_height_ch[i,pp,:] = pq.GetPulseMaxHeightChannel(p_start[i,pp], p_end[i,pp], ch_data_sum_pulse_bls)
                #p_coincidence[i,pp] = 
                
                # Centroids
                center_bot_x[i,pp], center_bot_y[i,pp], center_top_x[i,pp], center_top_y[i,pp] = pq.GetCentroids(p_area_ch[i,pp])
                
            
            # ====================================================================================================================================
            # Event level analysis. This needs some serious work

            if n_pulses[i] != 0:
                waveform_area[i] = np.sum(ch_data_sum_pulse_bls[-1])
            else:
                waveform_area[i] = np.sum(ch_data_phdPerSample[-1,i-j*block_size,:])
            p_class[i,:] = pc.ClassifyPulses(p_tba[i, :], (p_afs_50[i, :]-p_afs_2l[i, :])*tscale, n_pulses[i], p_area[i,:]) # classifier

            # Look at events with both S1 and S2.
            index_s1 = (p_class[i,:] == 1) + (p_class[i,:] == 2) # S1's
            index_s2 = (p_class[i,:] == 3) + (p_class[i,:] == 4) # S2's
            n_s1[i] = np.sum(index_s1)
            n_s2[i] = np.sum(index_s2)
            if n_s1[i] > 0:
                sum_s1_area[i] = np.sum(p_area[i, index_s1])
            if n_s2[i] > 0:
                sum_s2_area[i] = np.sum(p_area[i, index_s2])
            if n_s1[i] > 0 and (n_s1[i] + n_s2[i]) > 1 and ((p_class[i,0] == 1) + (p_class[i,0] == 2)):
                if p_area[i,index_s1].size != 0 and p_area[i,index_s2].size != 0:
                    for ps in range(max_pulses):
                        if p_area[i,ps] == np.max(p_area[i,index_s1]):
                            index_max_s1[i] = ps
                        if p_area[i,ps] == np.max(p_area[i,index_s2]):
                            index_max_s2[i] = ps
                if p_area[i, index_max_s1[i]] > 100 and p_area[i, index_max_s2[i]] > 100:
                    drift_Time_max[i] = tscale*(p_afs_1[i, index_max_s2[i]]-p_afs_1[i, index_max_s1[i]])
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

            afs50_2 = (p_afs_50[i,:]-p_afs_2l[i,:])*tscale # wtf 
            
            # Condition to skip the individual plotting, hand scan condition
            if np.size(handscan) == 1: 
                plotyn = handscan
            else:
                plotyn = handscan[i]
           
            # Condition to include a wfm in the average
            add_wfm = np.any((p_area[i,:]>5000)*(p_tba[i,:]<-0.75))*(n_s1[i]==1)*(n_s2[i]==0)
            if add_wfm and save_avg_wfm:
                plotyn = add_wfm # in avg wfm mode, plot the events which will go into the average
                avg_wfm += ch_data_phdPerSample[-1,i-j*block_size,:]
                n_wfms_summed += 1

            # Both S1 and S2 condition
            s1s2 = (n_s1[i] == 1)*(n_s2[i] == 1)

            if inn == 's': sys.exit()
            if not inn == 'q' and plotyn and plot_event_ind == i:

                fig = pl.figure()
                ax = pl.gca()
                pl.plot(x*tscale, ch_data_phdPerSample[-1,i-j*block_size,:],color='black',lw=0.7, label = "Summed All" )
                pl.plot(x[:-1]*tscale, np.diff(ch_data_phdPerSample[-1,i-j*block_size,:]),"blue", label="Derivative")
                pl.xlabel(r"Time [$\mu$s]")
                pl.ylabel("phd/sample")
                pl.title("Event {}".format(i))
                for ps in range(n_pulses[i]):
                    pl.axvspan(tscale*p_start[i,ps],tscale*p_end[i,ps],alpha=0.3,color=pulse_class_colors[p_class[i,ps]],zorder=0)
                    pl.axvline(tscale*p_afs_1[i,ps],color=pulse_class_colors[p_class[i,ps]],zorder=0,linestyle='--')
                    ax.text((p_end[i,ps]) * tscale, (0.94-ps*0.2) * ax.get_ylim()[1], '{:.2f} phd'.format(p_area[i, ps]),
                            fontsize=9, color=pulse_class_colors[p_class[i, ps]])
                    ax.text((p_end[i,ps]) * tscale, (0.9-ps*0.2) * ax.get_ylim()[1], 'TBA={:.2f}'.format(p_tba[i, ps]),
                        fontsize=9, color=pulse_class_colors[p_class[i, ps]])
                    ax.text((p_end[i,ps]) * tscale, (0.86-ps*0.2) * ax.get_ylim()[1], 'Rise={:.2f} us'.format(afs50_2[ps]),
                        fontsize=9, color=pulse_class_colors[p_class[i, ps]])
                    #ax.text((end_times[ps]) * tscale, (0.82-ps*0.2) * ax.get_ylim()[1], 'Check={}'.format(temp_condition[ps]),
                    #    fontsize=9, color=pulse_class_colors[p_class[i, ps]])    
                pl.legend()
                pl.grid(which="both",axis="both",linestyle="--")
                pl.xlim(0,event_window)
                pl.show()
                inn = input("Press enter to continue, q to stop plotting, evt # to skip to # (forward only)")
                #fig.clf()

            counter += 1
                
        # end of loop over events

        if save_avg_wfm:
            avg_wfm /= n_wfms_summed
            np.savetxt(data_dir+'average_waveform.txt',avg_wfm)
            print("Average waveform saved")

        n_events = i
        t_end = time.time()
        print("total number of events processed:", n_events)
        #print("Time used: {}".format(t_end-t_start))
        print("empty events: {0}".format(np.sum(empty_evt_ind>0)))

        j += 1

    # end of loop over compressed files


    # ==========================================================================================================================
    # Final step: save rq's

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
    list_rq['drift_Time_max'] = drift_Time_max
    list_rq['index_max_s1'] = index_max_s1
    list_rq['index_max_s2'] = index_max_s2
    list_rq['p_max_height'] = p_max_height
    list_rq['p_min_height'] = p_min_height
    list_rq['p_width'] = p_width
    list_rq['p_afs_1'] = p_afs_1
    list_rq['p_afs_2l'] = p_afs_2l
    list_rq['p_afs_10'] = p_afs_10
    list_rq['p_afs_50'] = p_afs_50
    list_rq['p_afs_90'] = p_afs_90
    list_rq['p_afs_99'] = p_afs_99   
    list_rq['p_hfs_10l'] = p_hfs_10l
    list_rq['p_hfs_50l'] = p_hfs_50l
    list_rq['p_hfs_90l'] = p_hfs_90l
    list_rq['p_hfs_10r'] = p_hfs_10r
    list_rq['p_hfs_50r'] = p_hfs_50r
    list_rq['p_hfs_90r'] = p_hfs_90r   
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
    list_rq['ev_time_s'] = rq_ev_time_s
    list_rq['p_max_height_ch'] = p_max_height_ch
    list_rq['right_area'] = right_area
    #list_rq[''] =    #add more rq

    #remove zeros in the end of each RQ array. 
    for rq in list_rq.keys():
        if rq != 'n_events':
            list_rq[rq] = list_rq[rq][:n_events]

    if filtered:
        save_name = "/rq_filtered.npy"
    else:
        save_name = "/rq.npy"
    rq = open(data_dir + save_name,'wb')
    np.savez(rq, **list_rq)
    rq.close()

  

def main():
    with open(sys.path[0]+"/path.txt", 'r') as path:
        data_dir = path.read()
    make_rq(data_dir)

if __name__ == "__main__":
    main()