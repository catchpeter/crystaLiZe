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


def gaussian(x,a,mu,sigma):
    return a*np.exp(-(x-mu)**2/(2*sigma**2))



#data_dir = "/home/xaber/Data/20220304/202203041344SPE/"

data_dir = "/home/xaber/Data/20220321/202203211637SPE/"

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
wsize = 3000+8  # samples per waveform # 12500 for 25 us
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
for bd in range(n_boards):
    # Loop over channels
    for ch in range(n_sipms[bd]):

        print("Channel "+str(n_order+1)+"/32")

        # Open data
        filename = data_dir+"waveforms_"+str(bd)+"_"+str(ch)+".dat"
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
            
            p_sarea[ev] = np.sum(ch_data[ev,1350:1500])

            
            # Event plotter
            if not inn == "q" and plotyn and ch>3 and p_sarea[ev] < 0.1/tscale and p_sarea[ev] > 0.01/tscale:
                fig = pl.figure(1,figsize=(10, 7))
                pl.grid(b=True,which='major',color='lightgray',linestyle='--')
                pl.plot(t, ch_data[ev,:], color="blue")
                #pl.plot(t, data_conv, color="red")
                #for ps in range(min(n_pulses,max_pulses)):
                #    pl.axvspan(start_times[ps]*tscale, end_times[ps]*tscale, alpha=0.25, color="green", label="Pulse area = {:0.3f} mV*ns".format(p_area[ev,ps]*tscale*1000))
                #    pl.axvspan(baselines_start[ps]*tscale, baselines_end[ps]*tscale, alpha=0.25, color='purple', label="baseline")
                pl.hlines(0,0,6,color="black")
                pl.xlim([0,6])
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

        clean_sarea = p_sarea*tscale

        image_basename = "_"+str(bd)+"_"+str(ch)+".png"


        nBins = int(100*(0.060-0.025)/(0.2+0.02))

        vals, bins = np.histogram(clean_sarea, bins=nBins, range=(0.025,0.060) )
        binsC = np.array([0.5*(bins[i] + bins[i+1]) for i in range(len(bins)-1)])
        popt, pcov = curve_fit(gaussian, binsC, vals, p0=(100,0.04,0.02))

        x = np.linspace(0.025,0.060,10000)
        y = gaussian(x, *popt)


        pl.figure()
        pl.hist(clean_sarea, bins=100, range=(-0.02,0.20), histtype="step")
        pl.plot(x,y,"r")
        pl.annotate("$\mu$ = ("+str( round(popt[1],3) )+"$\pm$"+str( round(np.sqrt(pcov[1,1]),3)) +") mV*$\mu$s", xy=(0.1,1000))
        pl.annotate("$\sigma$ = ("+str( round(np.abs(popt[2]),3) )+"$\pm$"+str( round(np.sqrt(pcov[2,2]),3)) +") mV*$\mu$s", xy=(0.1,600))
        #pl.hist(clean_sarea, bins=200, range=(-0.015,0.040), histtype="step")
        #pl.xlim([0,10])
        #pl.ylim([0,500])
        pl.title("Board "+str(bd)+", Channel "+str(ch))
        pl.yscale("log")
        pl.grid("both","both")
        pl.xlabel("Pulse area [mV*us]")
        pl.savefig(data_dir+"area"+image_basename)
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

             



                
  