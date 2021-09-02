#This will produce S1 peak 

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
from scipy.optimize import curve_fit

def peak_dt(data_dir, start_dt, end_dt, no_data, fname):
	sum_index = [0, 1, 3, 4, 5, 6]



	#load the RQs
	listrq = np.load(data_dir+'rq.npz')

	n_events = listrq['n_events'][()]
	n_pulses = listrq['n_pulses']
	n_s1 = listrq['n_s1']
	n_s2 = listrq['n_s2']
	p_area = listrq['p_area']
	p_class = listrq['p_class']
	drift_Time = listrq['drift_Time']
	p_tba = listrq['p_tba']
	p_start = listrq['p_start']

	listrq.close()

	cut_dict = {}
	cut_dict['S1'] = (p_class == 1) + (p_class == 2)
	
	# fitting algorithm
	def gaussian(x, mean, std,a):
		return a*np.exp(-((x-mean)/std)**2)

	# iterative Gaussian fit; first fit full range, then just mean +/- range_sig*sigma
	def gauss_fit(data, bins, p0=[100,30,300], range_sig=1.5):
		# data: bin heights from histogram
		# bins: bin edges (as output by plt.hist())
		bin_centers = bins[:-1] + np.diff(bins)/2
		start_bin = 0
		end_bin = -1

		popt, pcov = curve_fit(gaussian,bin_centers[start_bin:end_bin],data[start_bin:end_bin], p0=p0)

		start = popt[0]-range_sig*popt[1]
		end = popt[0]+range_sig*popt[1]
		start_bin = np.digitize(start, bins)
		end_bin = np.digitize(end, bins)
		popt, pcov = curve_fit(gaussian,bin_centers[start_bin:end_bin],data[start_bin:end_bin], p0=[popt[0],popt[1],popt[2]])
		return popt, pcov, start_bin, end_bin

	# plot S1 area at different drift times, set the range of drift time to plot
	drift_time_range = np.linspace(start_dt, end_dt, no_data)

	# save drift time and peak positon
	dt = np.zeros(no_data)
	peak = np.zeros(no_data)

	plot_col = 3
	plot_row = int(np.ceil(no_data/plot_col))


	fig = plt.figure(figsize = (30,35))
	gs = fig.add_gridspec(plot_row, plot_col, hspace=0.02, wspace=0.02)
	ax = gs.subplots(sharex='col', sharey=False)

	for i in range(plot_row):
		for j in range(plot_col):
			try: 
				min_dt = drift_time_range[j+i*plot_col]
				max_dt = drift_time_range[j+i*plot_col+1]
			except IndexError:
				break

			bin_n = 50
			range_set = (0, 300)
			cut_dict['s1_dt'] = np.transpose(np.tile((drift_Time>min_dt)*(drift_Time<max_dt), (4, 1)))*cut_dict['S1']
			cleanS1Area_SS = p_area[cut_dict['s1_dt']].flatten()
			h, b, _= ax[i, j].hist(cleanS1Area_SS, bins=bin_n, histtype = 'step', range = range_set)

			try:
				# iterative Gaussian fit to subtracted data
				popt, pcov, start_bin, end_bin = gauss_fit(h, b, p0=[170,30,100], range_sig=2.0) # p0=[mean, sigma, height]
				perr = np.sqrt(np.diag(pcov))

				bin_centers = b[:-1] + np.diff(b)/2
				ax[i, j].plot(bin_centers[start_bin:end_bin], gaussian(bin_centers[start_bin:end_bin], *popt), label ="{0:.1f}<DT<{1:.1f} Co".format(min_dt, max_dt)+'\nfit, after bkg sub\nmean: {0:.1f} +/- {1:.1f}'.format(popt[0],perr[0]))
				#save the drift time and peak position
				dt[j+i*plot_col] = (min_dt + max_dt)/2.
				peak[j+i*plot_col] = popt[0]
			except ValueError:
				pass
			except TypeError:
				pass

			ax[i, j].legend(fontsize = "x-small")
			#ax[i,j].xlabel("Pulse area (phd)")
			
			
	fig.suptitle("Co S1", fontsize=24)
	plt.savefig(fname)

	return peak, dt