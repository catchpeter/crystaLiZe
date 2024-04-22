import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os
import glob


#data_dir = '/Users/peter/Public/data/20240314-095146/'
#data_dir = './20240314-185842/' # garbage
data_dir = '/Users/peter/Public/data/20240314-190415/'
	
if 1: # 
	# Load headers and calculate event time
	try:
		h_file = np.load(data_dir+"/compressed_filtered_data/headers.npz")
		h_array = h_file["arr_0"]
		h_n_events = int(np.floor(h_array.size/8))
		h_array = np.reshape(h_array,(h_n_events,8))
	except:
		print("Error in loading header file!")
		h_n_events = n_events
		h_array = np.zeros((h_n_events,8))

	# Calculate event time
	# Precision up to 0.5 ms. To-do: get precision to 16 ns
	second_16 = h_array[:,5]
	second_16[second_16 < 0] = second_16[second_16 < 0] + 2**15
	second_16_next = np.zeros(h_n_events,dtype=int)
	for i in range(1,h_n_events):
		if second_16[i] - second_16[i-1] < 0:
			second_16_next[i:] += 1
	ev_time_s = 16*(second_16 + second_16_next*2**15)*(10**-9 * 2**15)
	
if 1: # specific to 3-fold cascade trigger
	shape1 = 4 # number of cascades in trigger
	shape0 = int(np.floor(ev_time_s.shape[0]/shape1))
	cascade_times_s = np.reshape(ev_time_s[0:shape0*shape1],( shape0 , shape1 ) )
	pl.figure(9);pl.clf()
	pl.plot(ev_time_s,'k.')
# 	ind=np.array([32,33]);pl.plot(ind,ev_time_s[ind],'r.')
# 	ind=np.array([46,47,48]);pl.plot(ind,ev_time_s[ind],'r.')
# 	ind=np.array([85,86]);pl.plot(ind,ev_time_s[ind],'r.')

if 1:
	diffs = np.diff(cascade_times_s)
	pl.figure(1);pl.clf()
	binn = np.arange(0,0.1,0.0001)
	pl.hist(diffs[:,0],binn,color='r')
	pl.hist(diffs[:,1],binn,color='g')
	pl.hist(diffs[:,2],binn,color='b')
	pl.yscale('log')
	
	
	
if 1:
	compressed_file_list = glob.glob(data_dir+"./compressed_filtered_data/c*.npz")
	for compressed_file in compressed_file_list:
		#print(f"Loading compressed file {j}/{len(compressed_file_list)-1}")
		# load data
		try:
			with np.load(compressed_file) as data:
				ch_data_adcc = data["arr_0"]
		except:
			print("Error in loading "+compressed_file)
			continue
		
		# Some important quantities
		vscale = (500.0/16384.0) #get_vscale(data_dir)
		event_window = 100 #get_event_window(data_dir)
		cascade = np.array([0.5,1.0,5])
		
		wsize = int(500 * event_window)  # samples per waveform # 12500 for 25 us
		tscale = (8.0/4096.0)     # = 0.002 µs/sample, time scale
		t = np.arange(0,wsize)*tscale # time array for plotting
		n_sipms = 32 # number of SiPMs
		n_channels = n_sipms + 1 # include sum
		block_size = int(1500*15/event_window) # number of events per compressed file
		
		n_tot_samp_per_ch = int( (ch_data_adcc.size)/n_sipms )
		n_events_b = int((ch_data_adcc.size)/(n_sipms*wsize)) # n events per compressed file (same as block_size)

		# Convert from ADCC to phd/sample and get summed waveform
		ch_data_adcc = np.concatenate((ch_data_adcc, np.zeros(n_tot_samp_per_ch) ))
		ch_data_mV = vscale*np.reshape(ch_data_adcc, (n_channels,n_events_b,wsize))
		#ch_data_mV = vscale*np.reshape(ch_data_adcc[0:n_channels*n_events_b*wsize], (n_channels,n_events_b,wsize))
		#ch_data_mV.shape = (33, 166, 10000) = (ch,evt,samples)

		#spe_dir = "/home/xaber/crystalize/Analysis/spe_calibration/202403/50V_3-6-2024.txt"
		# mV ns
		spe_dir = "/Users/peter/Desktop/crops/50V_3-6-2024.txt"
		spe_sizes = np.loadtxt(spe_dir, dtype='float')
		spe_sizes[spe_sizes==0] = 1e9
		
		np.mean(spe_sizes[spe_sizes>0])
		np.std(spe_sizes[spe_sizes>0])



for n in np.arange(0,int(ch_data_mV.shape[1]),4):

	print(n)
	pl.figure(8);pl.clf()
	
	pl.subplot(1,4,1)
	for i in range(0,32):
		pl.plot(t,ch_data_mV[i,n+0,:],'-',linewidth=0.5)
	pl.xlabel(r'$\mu$s');pl.ylabel('mV');pl.title('event %d'%n);
	ax=pl.gca();ax.set_ylim([-10,100])#;ax.set_xlim([0,20])
	
	pl.subplot(1,4,2)
	for i in range(0,32):
		pl.plot(t,ch_data_mV[i,n+1,:],'-',linewidth=0.5);
	pl.xlabel(r'$\mu$s');pl.ylabel('mV');pl.title('%1.1f ms delayed'%cascade[0]);
	ax=pl.gca();ax.set_ylim([-0.2,1])
	
	pl.subplot(1,4,3)
	for i in range(0,32):
		pl.plot(t,ch_data_mV[i,n+2,:],'-',linewidth=0.5);
	pl.xlabel(r'$\mu$s');pl.ylabel('mV');pl.title('%1.1f ms delayed'%cascade[1]);
	ax=pl.gca();ax.set_ylim([-0.2,1])
	
	pl.subplot(1,4,4)
	for i in range(0,32):
		pl.plot(t,ch_data_mV[i,n+3,:],'-',linewidth=0.5);
	pl.xlabel(r'$\mu$s');pl.ylabel('mV');pl.title('%1.1f ms delayed'%cascade[2]);
	ax=pl.gca();ax.set_ylim([-0.2,1])
	
	
	pl.show();pl.pause(0.1)	
	input('press any key to continue...')