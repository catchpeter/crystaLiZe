import numpy as np
import matplotlib.pyplot as pl
import os
import sys

from read_settings import get_event_window, get_vscale



def quickDraw(data_dir, howMany = 10, showFig = False):

    # Some important quantites
    vscale = get_vscale(data_dir)
    event_window = get_event_window(data_dir)
    wsize = int(500 * event_window)  # samples per waveform # 12500 for 25 us
    tscale = (8.0/4096.0)     # = 0.002 µs/sample, time scale
    n_sipms = 32    
    n_channels = n_sipms + 1 # include sum

    pulse_class_colors = np.array(['blue', 'green', 'red', 'magenta', 'darkorange'])

    # These values are close enough for plotting...
    # SPE sizes in SOLID, Oct 2022
    spe_sizes_0 = np.array([90.836,93.329,90.721,90.831,93.071,91.682,93.485,95.265,88.747,91.275,88.771,89.520,93.875,94.136,94.966,94.632])
    spe_sizes_1 = np.array([94.666,95.533,93.915,99.042,97.783,94.895,97.134,97.501])
    spe_sizes_2 = np.array([94.553,95.514,96.554,96.465,96.711,96.920,95.460,95.705])
    spe_sizes = np.concatenate((spe_sizes_0,spe_sizes_1,spe_sizes_2))

    # Load txt file of events to plot
    #try: 
    eventsToDraw = np.loadtxt(data_dir+"handscan.txt", dtype=int)
    #except:
    #    print("Error in loading handscan.txt")
    #    return
 

    # Load rq file
    try:
        listrq = np.load(data_dir+"rq_filtered.npy")
    except:
        print("Error in loading rq file")
        return
    p_area = listrq['p_area']
    n_pulses = listrq['n_pulses']
    p_start = listrq['p_start']
    p_end = listrq['p_end']
    p_class = listrq['p_class']
    p_afs_1 = listrq['p_afs_1']
    p_tba = listrq['p_tba']
    p_afs_50 = listrq['p_afs_50']

    
    # Create handscan directory
    command = "mkdir "+data_dir+"/handscan"
    try:
        os.system(command)
    except:
        print("Handscan directory already exists")

    # Determine which compressed files these events are in
    block_size = int(1500*15/event_window) # Number of events loaded, then saved per compressed file
    compressedFilesToLoad = np.floor(eventsToDraw/block_size).astype("int")
    whereInCompressedFile = np.mod(eventsToDraw,block_size).astype("int")

    # Loop over events
    for i in range(eventsToDraw.size):

        if i >= howMany: continue
        print("Event "+str(eventsToDraw[i]))

        # Load compressed data
        #try:
        with np.load(data_dir+"compressed_filtered_data/compressed_filtered_"+str(compressedFilesToLoad[i])+".npz") as data:
            ch_data = data["arr_0"]
        #except:
        #    print("Error in loading compressed_"+str(compressedFilesToLoad[i])+".npy")
        #    continue

        # Create array of events
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

        # Get specific waveform make a plot
        waveform = v_bls_matrix_all_ch[-1,whereInCompressedFile[i],:]
        t = np.arange(0,wsize)*tscale
        evN = eventsToDraw[i]

        pl.rcParams.update({'font.size': 28})
        fig = pl.figure(figsize=(30,10))
        ax = pl.gca()

        pl.plot(t,waveform,color='black',lw=0.7, label = "Summed All")
        
        pl.xlabel(r"Time [$\mu$s]")
        pl.ylabel("phd/sample")
        pl.title("Event {}".format(evN))

        for ps in range(n_pulses[evN]):
                
            pl.axvspan(tscale*p_start[evN,ps],tscale*p_end[evN,ps],alpha=0.3,color=pulse_class_colors[p_class[evN,ps]],zorder=0)
            pl.axvline(tscale*p_afs_1[evN,ps],color=pulse_class_colors[p_class[evN,ps]],zorder=0,linestyle='--')
                    
            ax.text((p_end[evN,ps]) * tscale, (0.94-ps*0.2) * ax.get_ylim()[1], '{:.1f} phd'.format(p_area[evN, ps]),
                fontsize=22, color=pulse_class_colors[p_class[evN, ps]])
            ax.text((p_end[evN,ps]) * tscale, (0.9-ps*0.2) * ax.get_ylim()[1], 'TBA={:.1f}'.format(p_tba[evN, ps]),
                fontsize=22, color=pulse_class_colors[p_class[evN, ps]])
            ax.text((p_end[evN,ps]) * tscale, (0.86-ps*0.2) * ax.get_ylim()[1], 'Rise={:.1f} us'.format(p_afs_50[evN,ps]),
                fontsize=22, color=pulse_class_colors[p_class[evN, ps]])
        

        pl.grid(which="both",axis="both",linestyle="--")
        pl.xlim(0,event_window)
        pl.legend()

        pl.savefig(data_dir+"handscan/Event_"+str(eventsToDraw[i])+".png")
        if showFig: pl.show() 
        pl.close()
    


    return



def main():

    #data_dir = sys.argv[1]
    #"/media/xaber/extradrive4/crystalize_data/data-202211/20221108/20221108-1327_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_54SiPM_1.41bar_-101.59ICVbot_2fold_degraded_CsSide_120min/"
    #data_dir = "/media/xaber/extradrive2/crystalize_data/data-202303/20230318/20230318-0642_2DR_10mVtrig_20us_5202.0C_5002.0G_500A_54SiPM_1.51bar_78.42ICVbot_2fold_degradedNew_60min/"
    data_dir = "/media/xaber/G-Drive2/crystalize_data/data-202304/20230417/20230417-0933_2DR_10mVtrig_15us_2701.0C_2501.0G_500A_54SiPM_0.67bar_-119.71ICVbot_2fold_plainMesh_solid_CoOCVTop_10min/"
    data_dir = "/media/xaber/G-Drive2/crystalize_data/data-202304/20230417/20230417-1015_2DR_10mVtrig_15us_3202.0C_3001.0G_500A_54SiPM_0.69bar_-119.71ICVbot_2fold_plainMesh_solid_CoOCVTop_20min/"
    quickDraw(data_dir, howMany=100, showFig=True)

    return



if __name__ == "__main__":
    main()