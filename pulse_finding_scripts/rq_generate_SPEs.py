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


def find_SPEs(data_dir, handscan=False, max_pulses=10, filtered=True, correct_swap=False):
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


    #spe_area = []
    #spe_height = []
    #spe_rms = []
    #spe_width = []


    spe_area = np.zeros((n_sipms,n_events,max_pulses))
    spe_height = np.zeros((n_sipms,n_events,max_pulses))
    spe_rms = np.zeros((n_sipms,n_events,max_pulses))
    spe_width = np.zeros((n_sipms,n_events,max_pulses))
    spe_start = np.zeros((n_sipms,n_events,max_pulses))
    spe_end = np.zeros((n_sipms,n_events,max_pulses))

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
                
                keep_looking = True
                while keep_looking:

                    peak = max(ch_data_mV[k,i-j*block_size,:])
                    peak_i = np.argmax(ch_data_mV[k,i-j*block_size,:])
                    if peak < 0 or peak > 3 or peak_i < 205 or peak_i > wsize-805: break

                
                    beg = peak_i - 120
                    end = peak_i + 280
                    spe_start[k,i,0] = beg
                    spe_end[k,i,0] = end
                    spe_area[k,i,0] = np.sum(tscale*1000*ch_data_mV[k,i-j*block_size,beg:end])
                    spe_height[k,i,0] = peak
                    spe_rms[k,i,0] = get_rms(ch_data_mV[k,i-j*block_size,beg:end])

                  

                    if avg_wf and spe_height[k,i,0] > 0.16 and spe_height[k,i,0] < spe_upper_height:
                        #print(ch_data_mV[k,i-j*block_size,peak_i-200:peak_i+500].shape, spe_avg[k,:].shape)
                        spe_sum[k,:] += ch_data_mV[k,i-j*block_size,peak_i-200:peak_i+800]
                        spe_nSaved[k] += 1


                    lh_cut = beg - 500
                    rh_cut = end + 500

                    if lh_cut > 205:

                        peak = max(ch_data_mV[k,i-j*block_size,:lh_cut])
                        peak_i = np.argmax(ch_data_mV[k,i-j*block_size,:lh_cut])
                        if peak < 0 or peak > 3: break

                        beg = peak_i - 60
                        end = peak_i + 110
                        spe_start[k,i,1] = beg
                        spe_end[k,i,1] = end
                        spe_area[k,i,1] = np.sum(tscale*1000*ch_data_mV[k,i-j*block_size,beg:end])
                        spe_height[k,i,1] = peak
                        spe_rms[k,i,1] = get_rms(ch_data_mV[k,i-j*block_size,beg:end])

                        """
                        if avg_wf and spe_height[k,i,0] > 0.2 and spe_height[k,i,0] < 0.45:
                            #print(ch_data_mV[k,i-j*block_size,peak_i-200:peak_i+500].shape, spe_avg[k,:].shape)
                            spe_sum[k,:] += ch_data_mV[k,i-j*block_size,peak_i-200:peak_i+500]
                            spe_nSaved[k] += 1
                        """

                    if rh_cut < wsize - 505:

                        peak = max(ch_data_mV[k,i-j*block_size,rh_cut:])
                        peak_i = np.argmax(ch_data_mV[k,i-j*block_size,rh_cut:])
                        if peak < 0 or peak > 3: break

                        beg = peak_i - 60 + rh_cut
                        end = peak_i + 110 + rh_cut
                        spe_start[k,i,2] = beg
                        spe_end[k,i,2] = end
                        spe_area[k,i,2] = np.sum(tscale*1000*ch_data_mV[k,i-j*block_size,beg:end])
                        spe_height[k,i,2] = peak
                        spe_rms[k,i,2] = get_rms(ch_data_mV[k,i-j*block_size,beg:end])

                        """
                        if avg_wf and spe_height[k,i,0] > 0.2 and spe_height[k,i,0] < 0.45:
                            #print(ch_data_mV[k,i-j*block_size,peak_i-200:peak_i+500].shape, spe_avg[k,:].shape)
                            spe_sum[k,:] += ch_data_mV[k,i-j*block_size,peak_i-200:peak_i+500]
                            spe_nSaved[k] += 1
                        """




                    keep_looking = False




            # Find max thing in pulse
            # If greater than some value, find pulse bounds
            # Cut left and right and repeat

            
 






          

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
            #plotyn = np.any(spe_area[]) #True #np.any( a2 > 7  )
            #plotyn = np.any( (se_area[i,:] > 0)&(se_area[i,:] < 2)&(se_width[i,:] > 1.3/tscale)&(se_coincidence[i,:] == 3)  )
            plotyn=False #np.any(spe_height[0,i,:] > 0.19)
            if not inn == 'q' and plotyn and plot_event_ind == i:
                #print((se_area[i,:] > 0)&(se_area[i,:] < 2)&(se_width[i,:] > 1.3/tscale)&(se_coincidence[i,:] == 3))

 
                fig = pl.figure()
                ax = pl.gca()
                #pl.pcolormesh(tspect, fspect, np.log10(Sxx), shading="gouraud",)
                #pl.plot(x*tscale, ch_data_mV[-1,i-j*block_size,:],color='black',lw=1.2, label = "Summed All" )
                for ch in range(0,1):
                    pl.plot(x*tscale, ch_data_mV[ch,i-j*block_size,:] )
                    for ps in range(3):
                        pl.axvspan(tscale*spe_start[ch,i,ps], tscale*spe_end[ch,i,ps],alpha=0.2,color="blue",zorder=0  )

  
  
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
        print(np.count_nonzero(se_coincidence))
        j += 1

    # end of loop over compressed files


    # ==========================================================================================================================
    # Final step: save rq's

    #create a dictionary with all RQs
    list_rq = {}
    list_rq['spe_area'] = spe_area
    list_rq['spe_height'] = spe_height
    list_rq['spe_rms'] = spe_rms
    list_rq['spe_width'] = spe_width
    list_rq['spe_sum'] = spe_sum 
    list_rq['spe_nSaved'] = spe_nSaved
   


    if filtered:
        save_name = "/rq_SPE_filtered_new.npy"
    else:
        save_name = "/rq_SPE_test.npy"
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



def get_rms(wf):

    return np.sqrt(np.sum(np.power(wf,2) )/wf.size )




def main():
    #with open(sys.path[0]+"/path.txt", 'r') as path:
    #    data_dir = path.read()
   
    #data_dir = "/media/xaber/G-Drive2/crystalize_data/data-202305/20230524/20230524-1122_0.5DR_10mVtrig_50us_8203.0C_8003.0G_1000A_54SiPM_1.68bar_-151.12ICVbot_2fold_plainMesh_liquid_10msDelay_1min/"


    #data_dir = "/media/xaber/G-Drive2/crystalize_data/data-202305/20230524/20230524-1639_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_54SiPM_1.68bar_-151.12ICVbot_2fold_plainMesh_liquid_BaOCVTop_delay500us_1min/"

    #data_dir = "/media/xaber/G-Drive2/crystalize_data/data-202305/20230524/20230524-1713_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_53SiPM_1.68bar_-151.12ICVbot_2fold_plainMesh_liquid_BaOCVTop_delay500us_1min/"

    #data_dir = "/media/xaber/G-Drive2/crystalize_data/data-202305/20230524/20230524-1717_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_52SiPM_1.68bar_-151.12ICVbot_2fold_plainMesh_liquid_BaOCVTop_delay500us_1min/"

    #data_dir = "/media/xaber/G-Drive2/crystalize_data/data-202305/20230524/20230524-1721_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_51SiPM_1.68bar_-151.12ICVbot_2fold_plainMesh_liquid_BaOCVTop_delay500us_1min/"

    #data_dir = "/media/xaber/G-Drive2/crystalize_data/data-202305/20230524/20230524-1726_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_50SiPM_1.68bar_-151.12ICVbot_2fold_plainMesh_liquid_BaOCVTop_delay500us_1min/"

    #data_dir = "/media/xaber/G-Drive2/crystalize_data/data-202305/20230524/20230524-1730_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_49SiPM_1.68bar_-151.12ICVbot_2fold_plainMesh_liquid_BaOCVTop_delay500us_1min/"

    #data_dir_list = glob.glob("/media/xaber/G-Drive2/crystalize_data/data-202305/20230524/*54SiPM*5min/")

    #data_dir_list = glob.glob("/media/xaber/G-Drive2/crystalize_data/data-202305/20230531/*_0.5DR*54SiPM*SPE*/")

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202305/20230531/20230531-1147_0.5DR_10mVtrig_20us_5202.0C_5002.0G_500A_54SiPM_1.7bar_-151.12ICVbot_2fold_SPE_500usDelay_plainMesh_liquid_BaTop_20min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202305/20230525/20230525-1317_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_51SiPM_1.68bar_-151.12ICVbot_2fold_plainMesh_liquid_BaOCVTop_delay500us_30min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202305/20230525/20230525-1455_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_52SiPM_1.68bar_-151.12ICVbot_2fold_plainMesh_liquid_BaOCVTop_delay500us_30min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202305/20230525/20230525-1616_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_40SiPM_1.68bar_-151.12ICVbot_2fold_plainMesh_liquid_BaOCVTop_delay500us_10min/"]

    #data_dir_list = glob.glob("/media/xaber/G-Drive2/crystalize_data/data-202305/20230525/*BaOCVTop_delay500us_30min/")

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230605/20230605-0907_0.5DR_10mVtrig_20us_5202.0C_5002.0G_500A_54SiPM_1.68bar_-151.12ICVbot_2fold_SPEtest_noAmpCh0_plainMesh_liquid_5min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230605/20230605-1024_0.5DR_10mVtrig_20us_5202.0C_5002.0G_500A_50SiPM_1.7bar_-151.12ICVbot_2fold_SPEtest_noAmpCh0_plainMesh_liquid_5min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230605/20230605-1030_0.5DR_10mVtrig_20us_5202.0C_5002.0G_500A_51SiPM_1.72bar_-150.82ICVbot_2fold_SPEtest_noAmpCh0_plainMesh_liquid_5min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230606/20230606-1204_0.5DR_10mVtrig_20us_5202.0C_5002.0G_500A_50SiPM_1.65bar_-151.12ICVbot_2fold_BaTop_noAmp_SPEdelay500us_plainMesh_liquid_5min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230605/20230605-1037_0.5DR_10mVtrig_20us_5202.0C_5002.0G_500A_52SiPM_1.7bar_-151.12ICVbot_2fold_SPEtest_noAmpCh0_plainMesh_liquid_5min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202305/20230531/20230531-1147_0.5DR_10mVtrig_20us_5202.0C_5002.0G_500A_54SiPM_1.7bar_-151.12ICVbot_2fold_SPE_500usDelay_plainMesh_liquid_BaTop_20min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202305/20230524/20230524-1124_0.5DR_10mVtrig_50us_8203.0C_8003.0G_1000A_54SiPM_1.68bar_-151.12ICVbot_2fold_plainMesh_liquid_10msDelay_20min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230622/20230622-1023_2DR_10mVtrig_6us_3202.0C_3001.0G_0A_50SiPM_1.89bar_-149.61ICVbot_2fold_SPE_contTrig_noAmp_plainMesh_liquid_1min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230622/20230622-1034_2DR_10mVtrig_30us_3202.0C_3001.0G_0A_52SiPM_1.91bar_-149.91ICVbot_2fold_SPE_delay500us_noAmp_plainMesh_liquid_1min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230622/20230622-1040_2DR_10mVtrig_30us_3202.0C_3001.0G_0A_52SiPM_1.89bar_-149.91ICVbot_2fold_SPE_delay500us_noAmp_plainMesh_liquid_2min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230622/20230622-1104_2DR_10mVtrig_30us_3702.0C_3502.0G_0A_54SiPM_1.82bar_-149.91ICVbot_2fold_SPE_delay500us_noAmp_plainMesh_liquid_2min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230622/20230622-1123_2DR_10mVtrig_30us_3702.0C_3501.0G_0A_54SiPM_1.58bar_-149.91ICVbot_2fold_SPE_delay500us_noAmp_plainMesh_liquid_2min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230622/20230622-1132_2DR_10mVtrig_30us_5202.0C_5002.0G_500A_54SiPM_1.58bar_-149.61ICVbot_2fold_SPE_delay500us_noAmp_plainMesh_liquid_2min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230622/20230622-1140_2DR_10mVtrig_30us_5202.0C_5002.0G_500A_54SiPM_1.55bar_-149.91ICVbot_2fold_SPE_delay500us_noAmp_plainMesh_liquid_2min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230622/20230622-1159_0.5DR_10mVtrig_20us_5202.0C_5002.0G_500A_54SiPM_1.58bar_-149.91ICVbot_2fold_SPE_delay500us_noAmp_plainMesh_liquid_2min/"]



    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230622/20230622-1555_0.5DR_10mVtrig_20us_5203.0C_5003.0G_500A_54SiPM_1.5bar_-149.61ICVbot_2fold_SPE_delay500us_noAmp_plainMesh_liquid_1min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230622/20230622-1628_0.5DR_10mVtrig_20us_5203.0C_5003.0G_500A_54SiPM_1.51bar_-149.61ICVbot_2fold_SPE_delay500us_noAmp_plainMesh_liquid_2min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230622/20230622-1646_0.5DR_10mVtrig_20us_5203.0C_5003.0G_500A_49SiPM_1.48bar_-149.91ICVbot_2fold_SPE_delay500us_noAmp_plainMesh_liquid_2min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230623/20230623-1457_0.5DR_10mVtrig_20us_5203.0C_5002.0G_500A_50SiPM_1.48bar_-149.91ICVbot_2fold_SPE_delay500us_BaTop_noAmp_plainMesh_liquid_30min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230627/20230627-1149_0.5DR_10mVtrig_20us_5203.0C_5003.0G_500A_51SiPM_1.43bar_-149.91ICVbot_2fold_BaTop_SPE_delay500us_noAmp_plainMesh_liquid_30min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230628/20230628-0735_0.5DR_10mVtrig_20us_5203.0C_5002.0G_500A_52SiPM_1.53bar_-149.91ICVbot_2fold_BaTop_SPE_delay500us_noAmp_plainMesh_liquid_30min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230628/20230628-1731_0.5DR_10mVtrig_20us_5203.0C_5003.0G_500A_53SiPM_1.5bar_-149.61ICVbot_2fold_BaTop_SPE_delay500us_noAmp_plainMesh_liquid_30min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202306/20230629/20230629-1241_0.5DR_10mVtrig_20us_5203.0C_5003.0G_500A_54SiPM_1.51bar_-149.91ICVbot_2fold_BaTop_SPE_delay500us_noAmp_plainMesh_liquid_30min/"]

    print("Sleeping before rq-ing")
    time.sleep(4*60*60)
    data_dir_list = glob.glob("/media/xaber/G-Drive2/crystalize_data/data-202307/20230711/*52SiPM*SPE*60min/")
    print("\n")
    for i in data_dir_list: print(i)
    print("\n")
    

    for data_dir in data_dir_list:
        print(data_dir)
        #compression(data_dir, save_everything = True, debug=False)
        find_SPEs(data_dir)

if __name__ == "__main__":
    main()