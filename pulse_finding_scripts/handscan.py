import numpy as np 
import matplotlib.pyplot as pl
import sys
import glob

from read_settings import get_event_window, get_vscale, get_sipm_bias


def handscan(data_dir, correct_swap=False):

    print(f"\nStarting handscan in: {data_dir}\n")

    # Some important quantites
    vscale = get_vscale(data_dir)
    event_window = get_event_window(data_dir)
    wsize = int(500 * event_window)  # samples per waveform # 12500 for 25 us
    tscale = (8.0/4096.0)     # = 0.002 µs/sample, time scale
    t = np.arange(0,wsize)*tscale # time array for plotting
    n_sipms = 32 # number of SiPMs
    n_channels = n_sipms + 1 # include sum
    block_size = int(1500*15/event_window) # number of events per compressed file

    # Load in SPE values
    spe_dir = "/home/xaber/crystalize/Analysis/spe_calibration/202306/"
    sipm_bias = get_sipm_bias(data_dir)
    spe_sizes = np.loadtxt(spe_dir+f"spes_{sipm_bias}.txt", dtype='float')

    # Loop over compressed files
    n_compressed_files = len(glob.glob(data_dir+"/compressed_filtered_data/compressed_filtered_*.npz"))
    for i in range(n_compressed_files):
        
        # Load data
        print(f"  Loading compressed file {i}")
        with np.load(data_dir+f"/compressed_filtered_data/compressed_filtered_{i}.npz") as data:
            ch_data_adcc = data["arr_0"]

        # Convert from ADCC to phd/sample and get summed waveform
        n_tot_samp_per_ch = int( (ch_data_adcc.size)/n_sipms )
        n_events_b = int((ch_data_adcc.size)/(n_sipms*wsize))
        ch_data_adcc = np.concatenate((ch_data_adcc, np.zeros(n_tot_samp_per_ch) ))
        ch_data_mV = vscale*np.reshape(ch_data_adcc, (n_channels,n_events_b,wsize))
        ch_data_phdPerSample = np.zeros_like(ch_data_mV) 
        for ch in range(n_sipms):
            if ch == 27: continue # dead SiPM
            if correct_swap:
                if ch == 0: print("    Correcting SiPM swap")
                if ch==7:
                    ch_data_phdPerSample[ch,:,:] = ch_data_mV[31,:,:]*tscale*(1000)/spe_sizes[ch]
                elif ch==31:
                    ch_data_phdPerSample[ch,:,:] = ch_data_mV[7,:,:]*tscale*(1000)/spe_sizes[ch]
                else:
                    ch_data_phdPerSample[ch,:,:] = ch_data_mV[ch,:,:]*tscale*(1000)/spe_sizes[ch]  
            else:
                ch_data_phdPerSample[ch,:,:] = ch_data_mV[ch,:,:]*tscale*(1000)/spe_sizes[ch]
        ch_data_phdPerSample[-1,:,:] = np.sum(ch_data_phdPerSample, axis=0)   


        # Loop over events
        for j in range(n_events_b):
            ev_num = j+i*block_size
            print(f"    Plotting event {ev_num}")

            # Make a plot
            fig = pl.figure()
            ax = pl.gca()
            for k in range(32):
                pl.plot(t, ch_data_phdPerSample[k,j,:],lw=0.7)
            pl.xlabel(r"Time [$\mu$s]")
            pl.ylabel("phd/sample")
            pl.title(f"Event {ev_num}")
            #pl.legend()
            pl.grid(which="both",axis="both",linestyle="--")
            pl.xlim(0,event_window)
            pl.show()


        # end of loop over events
    

    # end of loop over compressed files

    return



def main():

    data_dir = sys.argv[1]
    handscan(data_dir)

    return



if __name__ == "__main__":
    main()