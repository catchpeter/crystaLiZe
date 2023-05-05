import numpy as np
import matplotlib.pyplot as pl
import matplotlib as mpl
from scipy.optimize import curve_fit

from read_settings import get_event_window, get_vscale



# set plotting style
mpl.rcParams['font.size']=12
mpl.rcParams['axes.labelsize']=12
mpl.rcParams['legend.fontsize']='small'
mpl.rcParams['figure.autolayout']=True
mpl.rcParams['figure.figsize']=[8.0,6.0]

rms_cut_default = 0.4


def gaussian(x,a,mu,sigma):
    return a*np.exp(-(x-mu)**2/(2*sigma**2))


def spe_calibrate(data_dir, bd, ch, rms_cut=rms_cut_default, plot_wf=False, save_plots=True):


    load_dtype = "int16"

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


    inn = ""
    # Loop over events
    for ev in range(n_events):

        # Baseline subtraction
        ch_data[ev,:] = np.subtract(ch_data[ev,:], np.mean(ch_data[ev,:100]) )

        # Look at window by trigger
        l_bound = 380
        r_bound = 550
        rms[ev] = np.sqrt(np.sum(np.power(ch_data[ev,l_bound:r_bound],2) )/(r_bound-l_bound) )
        p_sarea[ev] = np.sum(ch_data[ev,l_bound:r_bound])

        # Option to plot
        if not inn == "q" and plot_wf:
            fig = pl.figure(1,figsize=(10, 7))
            pl.grid(b=True,which='major',color='lightgray',linestyle='--')
            pl.plot(t,ch_data[ev,:],"b")
            pl.axvline(l_bound*tscale,0,1,color="black")
            pl.axvline(r_bound*tscale,0,1,color="black")
            pl.hlines(0,0,2,color="black")
            pl.title("Event "+str(ev))
            pl.xlabel("Time [us]")
            pl.ylabel("mV")
            pl.show() 
            inn = input("Press enter to continue, q to stop plotting")


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
    if save_plots: pl.savefig(data_dir+"area"+image_basename)
    pl.close()

    pl.figure()
    pl.hist2d(clean_p_sarea, clean_rms, bins=full_bins, range=(full_range,[0.4,1]) ,  cmin=1, norm=mpl.colors.LogNorm() )
    pl.grid("both","both")
    pl.title("Board "+str(bd)+", Channel "+str(ch))
    pl.xlabel("Integrated area [mV*ns]")
    pl.ylabel("RMS [mV]")
    if save_plots: pl.savefig(data_dir+"area_rms"+image_basename)
    pl.close()

    
    return


             
def main():

    #data_dir = "/media/xaber/extradrive1/crystalize_data/data-202209/20220928/20220928-1726_0.5DR_NAmVtrig_2us_3201.0C_3001.0G_500A_54SiPM_1.2bar_-97.06ICVbot_SPE/"
    data_dir = "/media/xaber/G-Drive2/crystalize_data/data-202304/20230417/20230417-0838_0.5DR_NAmVtrig_2us_0.0C_0.0G_500A_54SiPM_0.69bar_-119.71ICVbot_SPE/"
    
    b0_ch = np.arange(0,16)
    b0_rms_cut = rms_cut_default*np.ones(16)
    #b0_rms_cut = [0.5,0.4,0.4,0.45,0.45,0.45,0.45,0.44,0.4,0.4,0.4,0.48,0.47,0.45,0.46,0.44]
    for ch in b0_ch:
        print("Board 0, channel "+str(ch) )
        spe_calibrate(data_dir, bd=0, ch=ch, rms_cut=b0_rms_cut[ch], plot_wf=False)

    b1_ch = np.arange(0,8)
    b1_rms_cut = rms_cut_default*np.ones(8)
    #b1_rms_cut = [0.47,0.5,0.5,0.4,0.4,0.5,0.47,0.4] 
    for ch in b1_ch:
        print("Board 1, channel "+str(ch) )
        spe_calibrate(data_dir, bd=1, ch=ch, rms_cut=b1_rms_cut[ch], plot_wf=False)

    b2_ch = np.arange(0,8)
    b2_rms_cut = rms_cut_default*np.ones(8)
    #b2_rms_cut = [0.4,0.46,0.4]
    for ch in b2_ch:
        print("Board 2, channel "+str(ch) )
        spe_calibrate(data_dir, bd=2, ch=ch, rms_cut=b2_rms_cut[ch], plot_wf=False)


    return 



if __name__ == "__main__":
    main()