import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mpl
import time
import sys
from scipy import stats
from scipy.optimize import curve_fit

import PulseFinderScipy as pf
import PulseQuantities as pq
import PulseClassification as pc

import PulseFinderSPE_crystalize as pfSPE

from read_settings import get_event_window, get_vscale


def gaussian(x,a,mu,sigma):
    return a*np.exp(-(x-mu)**2/(2*sigma**2))

"""

#data_dir = "/home/xaber/Data/20220304/202203041344SPE/"

data_dir = "/home/xaber/Data/data-202205/20220524/20220524-1624_0.5DR_NAmVtrig_2us_5202.0C_5002.0G_1000A_54SiPM_1.2bar_-102.19ICVbot_SPE/"

plotyn = False # Waveform plotting
saveplot = True # Save RQ plots

# set plotting style
#mpl.rcParams['font.size']=12
#mpl.rcParams['axes.labelsize']=12
#mpl.rcParams['legend.fontsize']='small'
#mpl.rcParams['figure.autolayout']=True
#mpl.rcParams['figure.figsize']=[8.0,6.0]

# ==================================================================




load_dtype = "int16"
array_dtype = "float32"




list_rq = {}



# DAQ parameters
wsize = 1000+8  # samples per waveform # 12500 for 25 us
vscale = (500.0/16384.0) # vertical scale, = 0.031 mV/ADCC for 0.5 Vpp dynamic range
tscale = (8.0/4096.0)     # = 0.002 µs/sample, time scale

# SiPM ordering
n_boards = 3
n_sipms = [16,8,8]
n_channels = np.sum(n_sipms)
n_order = 0

max_pulses = 4

inn=""

# Loop over boards
for bd in range(2,3): #range(n_boards):
    # Loop over channels
    for ch in range(6,7): #range(n_sipms[bd]):

       



        print("Channel "+str(n_order+1)+"/32")

        # Open data
        filename = data_dir+"b"+str(bd)+"/"+"waveforms_0_"+str(ch)+".dat"
        try:
            ch_data = np.fromfile(filename, dtype=load_dtype)
        except:
            print(filename)
            print("no data from bd "+str(bd)+" ch "+str(ch))
            continue 

        # Make nice np array and scale
        n_events = int(ch_data.size/wsize)
        ch_data = np.reshape(ch_data, (n_events, wsize))
        ch_data = ch_data[:,8:wsize]
        ch_data = np.multiply(ch_data, vscale)

        # Time array, units of us
        t = tscale*np.arange(0,wsize-8)

        # Initialize rq's
        p_start = np.zeros((n_events,max_pulses), dtype=int)
        p_end = np.zeros((n_events,max_pulses), dtype=int)
        p_width = np.zeros((n_events,max_pulses), dtype=int)
        p_area = np.zeros((n_events,max_pulses))
        p_sarea = np.zeros(n_events)
        p_max_height = np.zeros((n_events,max_pulses))
        rms = np.zeros(n_events)

        # Loop over events
        for ev in range(n_events):

            # Preliminary baseline subtraction
            ch_data[ev,:] -= np.mean(ch_data[ev,0:100])

            # Pulse finding
            #start_times, end_times, peaks, baselines, data_conv, baselines_start, baselines_end = pfSPE.findDarkSPEs(ch_data[ev,:])

            # Loop over found pulses, calculate rq's
            #n_pulses = len(start_times)
            #for ps in range(min(n_pulses,max_pulses)):
                #if start_times[ps] >= end_times[ps]: continue
                #p_start[ev,ps] = start_times[ps]
                #p_end[ev,ps] = end_times[ps]
                #p_width[ev,ps] =  end_times[ps] - start_times[ps]
                #p_area[ev,ps] = pq.GetPulseArea(start_times[ps], end_times[ps], ch_data[ev,:]-baselines[ps])
                
                #p_max_height[ev,ps] = pq.GetPulseMaxHeight(start_times[ps], end_times[ps], ch_data[ev,:]-baselines[ps])
            
            l_bound = 380
            r_bound = 550
            

            rms[ev] = np.sqrt(np.sum(np.power(ch_data[ev,l_bound:r_bound],2) )/(r_bound-l_bound) )

            if rms[ev] < 0.52: continue #0.425: continue


            p_sarea[ev] = np.sum(ch_data[ev,l_bound:r_bound])
            
            # Pulse = 0.0011


            # Event plotter
            if not inn == "q" and p_sarea[ev] > 9999999.00: #plotyn and ch>3 and p_sarea[ev] < 0.1/tscale and p_sarea[ev] > 0.01/tscale:
                print(p_sarea[ev])
                print(rms)
                print()
                fig = pl.figure(1,figsize=(10, 7))
                pl.grid(b=True,which='major',color='lightgray',linestyle='--')
                pl.plot(t, ch_data[ev,:], color="blue")
                #pl.plot(t, data_conv, color="red")
                #for ps in range(min(n_pulses,max_pulses)):
                #    pl.axvspan(start_times[ps]*tscale, end_times[ps]*tscale, alpha=0.25, color="green", label="Pulse area = {:0.3f} mV*ns".format(p_area[ev,ps]*tscale*1000))
                #    pl.axvspan(baselines_start[ps]*tscale, baselines_end[ps]*tscale, alpha=0.25, color='purple', label="baseline")
                pl.hlines(0,0,2,color="black")
                pl.axvline(l_bound*tscale,0,1,color="black")
                pl.axvline(r_bound*tscale,0,1,color="black")
                pl.xlim([0,2])
                pl.title("Board "+str(bd)+", Channel "+str(ch)+", Event "+str(ev))  
                pl.xlabel('Time (us)')
                pl.ylabel('mV')
                pl.legend(loc="upper right")
                pl.draw()
                #pl.ion()
                pl.show()
                pl.show(block=0)
                inn = input("Press enter to continue, q to stop plotting, evt # to skip to # (forward only)")
                pl.close()
                #fig.clf()
            


        # Cleaning up, plotting, and saving
        #clean_cut = (p_area > 0)*(p_width > 0.05/tscale)


        #clean_area = p_area[clean_cut]*tscale 
        #clean_width = p_width[clean_cut]*tscale 

        clean_sarea = p_sarea*tscale*1000

        image_basename = "_"+str(bd)+"_"+str(ch)+".png"

        full_range = (-100,300)
        full_bins = 200

        if bd == 0:
            fit_range = (60,105)
        else:
            fit_range = (65,100) #(55,100)
        fit_bins = int( full_bins*(fit_range[1] - fit_range[0])/(full_range[1] - full_range[0]) )

    
       
        vals, bins = np.histogram(clean_sarea, bins=fit_bins, range=fit_range)
        binsC = np.array([0.5*(bins[i] + bins[i+1]) for i in range(len(bins)-1)])
        popt, pcov = curve_fit(gaussian, binsC, vals, p0=(100,85,20))

        x = np.linspace(fit_range[0],fit_range[1],10000)
        y = gaussian(x, *popt)

        spe_h = popt[0]

        pl.figure()
        pl.hist(clean_sarea, bins=full_bins, range=full_range, histtype="step")
        pl.plot(x,y,"r")
        pl.annotate("$\mu$ = ("+str( round(popt[1],3) )+"$\pm$"+str( round(np.sqrt(pcov[1,1]),3)) +") mV*ns", xy=(100,spe_h*1.5))
        pl.annotate("$\sigma$ = ("+str( round(np.abs(popt[2]),3) )+"$\pm$"+str( round(np.sqrt(pcov[2,2]),3)) +") mV*ns", xy=(100,spe_h*1.4))
        #pl.xlim([0,10])
        pl.ylim([0,spe_h*2])
        pl.title("Board "+str(bd)+", Channel "+str(ch))
        #pl.yscale("log")
        pl.grid("both","both")
        pl.xlabel("Pulse area [mV*ns]")
        pl.savefig(data_dir+"area"+image_basename)
        pl.close()


        pl.figure()
        pl.hist2d(clean_sarea, rms, bins=full_bins, range=(full_range,[0.42,1]) ,  cmin=1, norm=mpl.colors.LogNorm() )
        pl.grid("both","both")
        pl.savefig(data_dir+"area_rms"+image_basename)
        pl.close()

        #pl.figure()
        #pl.hist(clean_width.flatten(), bins=100)
        #pl.savefig(data_dir+"width"+image_basename)
        #pl.close()

        #pl.figure()
        #pl.scatter(clean_area.flatten(), clean_width.flatten())
        #pl.xlim([0,6])
        #pl.savefig(data_dir+"area_vs_width"+image_basename)
        #pl.close()


        n_order += 1

"""

def spe_calibrate(data_dir, bd, ch, rms_cut=0.42):

    saveplot = True # Save RQ plots

    load_dtype = "int16"
    array_dtype = "float32"

    # Open data
    filename = data_dir+"b"+str(bd)+"/"+"waveforms_0_"+str(ch)+".dat"
    try:
        ch_data = np.fromfile(filename, dtype=load_dtype)
    except: 
        #print(filename)
        print("no data from bd "+str(bd)+" ch "+str(ch))
        return

    # DAQ parameters
    vscale = get_vscale(data_dir)
    event_window = get_event_window(data_dir)
    if event_window < 0: 
        print("Invalid event window")
        return
    wsize = int(500 * event_window) + 8  # samples per waveform # 12500 for 25 us
    tscale = (8.0/4096.0)     # = 0.002 µs/sample, time scale
    t = tscale*np.arange(0,wsize-8)

    # Make nice np array and scale
    n_events = int(ch_data.size/wsize)
    ch_data = np.reshape(ch_data, (n_events, wsize))
    ch_data = ch_data[:,8:wsize]
    ch_data = np.multiply(ch_data, vscale)

    # Initialize rq's
    rms = np.zeros(n_events)
    p_sarea = np.zeros(n_events)


    # Loop over events
    for ev in range(n_events):

        # Baseline subtraction
        ch_data[ev,:] = np.subtract(ch_data[ev,:], np.mean(ch_data[ev,:100]) )

        # Look at window by trigger
        l_bound = 380
        r_bound = 550
        rms[ev] = np.sqrt(np.sum(np.power(ch_data[ev,l_bound:r_bound],2) )/(r_bound-l_bound) )
        p_sarea[ev] = np.sum(ch_data[ev,l_bound:r_bound])


    # Scale data, [mV*ns]
    p_sarea *= 1000*tscale

    # Make a cut
    cut = rms > rms_cut
    clean_rms = rms[cut]
    clean_p_sarea = p_sarea[cut]

    # Find bounds for fit to single
    hist_l, bin_edges_l = np.histogram(clean_p_sarea, bins=20, range=(40,90))
    #hist_r, bin_edges_r = np.histogram(clean_p_sarea, bins=20, range=(90,120))
    ind_l = np.argmin(hist_l)
    #ind_r = np.argmin(hist_r)
    l_bound = bin_edges_l[ind_l]
    r_bound = 120 # bin_edges_r[ind_r[0]]

    # Bounds and ranges for fit 
    full_range = (-100,300)
    full_bins = 100 #200
    fit_range = (l_bound, r_bound)
    fit_bins = int( full_bins*(fit_range[1] - fit_range[0])/(full_range[1] - full_range[0]) )

    # Fit the single
    vals, bins = np.histogram(clean_p_sarea, bins=fit_bins, range=fit_range)
    binsC = np.array([0.5*(bins[i] + bins[i+1]) for i in range(len(bins)-1)])
    popt, pcov = curve_fit(gaussian, binsC, vals, p0=(100,85,20), maxfev=10000)
    x = np.linspace(fit_range[0],fit_range[1],10000)
    y = gaussian(x, *popt)
    spe_h = popt[0]

    # Make some plots
    image_basename = "_"+str(bd)+"_"+str(ch)+".png"
    pl.figure()
    pl.hist(clean_p_sarea, bins=full_bins, range=full_range, histtype="step")
    pl.plot(x,y,"r")
    pl.annotate("$\mu$ = ("+str( round(popt[1],3) )+"$\pm$"+str( round(np.sqrt(pcov[1,1]),3)) +") mV*ns", xy=(100,spe_h*1.5))
    pl.annotate("$\sigma$ = ("+str( round(np.abs(popt[2]),3) )+"$\pm$"+str( round(np.sqrt(pcov[2,2]),3)) +") mV*ns", xy=(100,spe_h*1.4))
    #pl.xlim([0,10])
    pl.ylim([0,spe_h*2])
    pl.title("Board "+str(bd)+", Channel "+str(ch))
    #pl.yscale("log")
    pl.grid("both","both")
    pl.xlabel("Integrated area [mV*ns]")
    pl.savefig(data_dir+"area"+image_basename)
    pl.close()

    pl.figure()
    pl.hist2d(clean_p_sarea, clean_rms, bins=full_bins, range=(full_range,[0.42,1]) ,  cmin=1, norm=mpl.colors.LogNorm() )
    pl.grid("both","both")
    pl.title("Board "+str(bd)+", Channel "+str(ch))
    pl.xlabel("Integrated area [mV*ns]")
    pl.ylabel("RMS [mV]")
    pl.savefig(data_dir+"area_rms"+image_basename)
    pl.close()

    
    return


             
def main():

    

    #data_dir = "/media/xaber/extradrive1/crystalize_data/data-202209/20220926/20220926-1629_0.5DR_NAmVtrig_2us_3201.0C_3001.0G_500A_54SiPM_1.32bar_-94.94ICVbot_SPE/"

    #data_dir = "/media/xaber/extradrive1/crystalize_data/data-202209/20220926/20220926-1810_0.5DR_NAmVtrig_2us_3202.0C_3001.0G_500A_54SiPM_1.36bar_-94.94ICVbot_SPE/"

    #data_dir = "/media/xaber/extradrive1/crystalize_data/data-202209/20220928/20220928-1151_0.5DR_NAmVtrig_2us_3202.0C_3001.0G_500A_54SiPM_1.2bar_-97.36ICVbot_SPE/"

    #data_dir = "/media/xaber/extradrive1/crystalize_data/data-202209/20220928/20220928-1437_0.5DR_NAmVtrig_2us_3202.0C_3001.0G_500A_54SiPM_1.22bar_-97.06ICVbot_SPE/"

    data_dir = "/media/xaber/extradrive1/crystalize_data/data-202209/20220928/20220928-1726_0.5DR_NAmVtrig_2us_3201.0C_3001.0G_500A_54SiPM_1.2bar_-97.06ICVbot_SPE/"

    b0_ch = np.arange(0,16)
    b0_rms_cut = [0.5,0.4,0.4,0.45,0.45,0.45,0.45,0.44,0.4,0.4,0.4,0.48,0.47,0.45,0.46,0.44]
    #0.42*np.ones(16)
    # 8,9 need re-doing
    for ch in [8,9]: #range(16):
        print("Board 0, channel "+str(ch) )
        #spe_calibrate(data_dir, bd=0, ch=ch, rms_cut=b0_rms_cut[ch])

    b1_ch = np.arange(0,8)
    #b1_rms_cut = 0.42*np.ones(8)
    b1_rms_cut = [0.47,0.5,0.5,0.4,0.4,0.5,0.47,0.4] 
    for ch in range(8):
        print("Board 1, channel "+str(ch) )
        spe_calibrate(data_dir, bd=1, ch=ch, rms_cut=b1_rms_cut[ch])

    b2_ch = np.arange(0,8)
    b2_rms_cut = 0.42*np.ones(8)
    #b2_rms_cut = [0.4,0.46,0.4]
    for ch in range(8):
        print("Board 2, channel "+str(ch) )
        spe_calibrate(data_dir, bd=2, ch=ch, rms_cut=b2_rms_cut[ch])


    return 



if __name__ == "__main__":
    main()