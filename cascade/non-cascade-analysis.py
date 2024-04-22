import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob
from scipy import stats as st

### sandboxing
if (socket.gethostname()[0]=='b') :
	data_dir = '/Users/peter/Public/data/20240315-150207/' # 133Ba
	spe_dir = "/Users/peter/Desktop/crops/50V_3-6-2024.txt"
else:
	data_dir = '/media/xaber/extradrive1/crystalize_data/data-202403/20240315/20240315-150207/'
	spe_dir = "/home/xaber/crystalize/Analysis/spe_calibration/202403/50V_3-6-2024.txt"

aa_dir = data_dir + "aa/" # output rqs "aa" here
try:
	os.mkdir(aa_dir)
except:
	print('')
	
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
	
	if 0: # specific to 3-fold cascade trigger
		shape1 = 4 # number of cascades in trigger
		shape0 = int(np.floor(ev_time_s.shape[0]/shape1))
		cascade_times_s = np.reshape(ev_time_s[0:shape0*shape1],( shape0 , shape1 ) )
	
	
	spe_sizes = np.loadtxt(spe_dir, dtype='float') # mV ns
	spe_sizes[spe_sizes==0] = 1e9

ff = 0	
compressed_file_list = glob.glob(data_dir+"./compressed_filtered_data/c*.npz")
compressed_file_list = sorted( compressed_file_list )
for compressed_file in compressed_file_list:
	if 1:
		print("loading: ... %s"%compressed_file_list[ff][42:])
		try:
			with np.load(compressed_file) as data:
				ch_data_adcc = data["arr_0"]
		except:
			print("Error in loading "+compressed_file)
			continue
	
		# Some important quantities
		vscale = (500.0/16384.0) #get_vscale(data_dir)
		event_window = 20 #get_event_window(data_dir)
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


	# buffer some memory for the RQs
	aa = np.zeros([32, ch_data_mV.shape[1]]) # areas of cathode_echo
	s1 = np.zeros([32, ch_data_mV.shape[1]]) # areas of s1
	s2 = np.zeros([32, ch_data_mV.shape[1]]) # areas of s1
	ee = np.zeros([32, ch_data_mV.shape[1]]) # areas of SE
	ge = np.zeros([32, ch_data_mV.shape[1]]) # areas of gate echo
	is1 = np.zeros([ch_data_mV.shape[1]]) # i_s1
	is2 = np.zeros([ch_data_mV.shape[1]]) # i_s2
	s2f = np.zeros([ch_data_mV.shape[1]]) # s2_frac

	print("analyzing compressed data...")

	sipms = np.concatenate((np.arange(0,27),np.arange(28,32))) # skip 27

	n_start = 0
	n_stop = ch_data_mV.shape[1]
	for n in np.arange(n_start,n_stop,1):
		i_max = int(np.mean(np.argmax(ch_data_mV[sipms,n+0,:],axis=1)))
		s=np.zeros(32,dtype=np.int16)
		for i in range(0,32): # refined edge-finding
			while (ch_data_mV[i,n+0,i_max-s[i]]>0.05): 
				s[i]=s[i]+1
				if s[i]>i_max:
					break
		i_s2 = int(i_max-np.mean(s))			
		i_s1 = np.argmax(np.sum(ch_data_mV[0:32,n+0,0:i_s2-150],axis=0),axis=0)-40 # 40 is the rise time
		i_ge = i_s2 + 1150 # by eye, index of max of gate echo
		i_ce = i_s2 + 3250 # by eye, index of max of cathode echo
		i_ee = i_s1 + 1125 # by eye
		i_w = 325

		is1[n+0] = i_s1
		is2[n+0] = i_s2
		s2f[n+0] = np.sum(np.sum(ch_data_mV[:,n+0,i_s2:i_s2+1000])) / np.sum(np.sum(ch_data_mV[:,n+0,:])) # approx S2 fraction

		for i in range(0,32): # now go an integrate individual channel areas
			y0 = ch_data_mV[i,n+0,:]

			beans = np.arange(-0.2,0.2,0.001)
			binc = (beans[0:-1] + beans[1:])/2
			blr = np.arange(0,2000) # only for initiating event in cascade
			bstd0 = np.std(y0[blr]) # get the std of the triggering event baseline
	
			# need weighted mean for baseline
			[counts0,bine] = np.histogram(y0[blr],beans)
			wm0 = np.average(binc,axis=0,weights=counts0)
			y0z = y0-wm0 # apply additional baseline-correction to data	

			
			aa[i,n+0] = np.sum(y0z[i_ce-i_w:i_ce+i_w]) * 2 / spe_sizes[i] # rough, but should be close
			s1[i,n+0] = np.sum(y0z[i_s1:i_s1+i_w]) * 2 / spe_sizes[i] # rough, but should be close
# 			s2[i,n+0] = np.sum(y0z[i_s2:i_s2+i_w]) * 2 / spe_sizes[i] # rough, but should be close
			s2[i,n+0] = np.sum(y0z[i_s2:i_ge-i_w]) * 2 / spe_sizes[i] # rough, but should be close
			ge[i,n+0] = np.sum(y0z[i_ge-i_w:i_ge+i_w]) * 2 / spe_sizes[i] # rough, but should be close
			if i_s2>(i_ee+i_w) and i_s1<(i_ee-i_w):
				ee[i,n+0] = np.sum(y0z[i_ee-i_w:i_ee+i_w]) * 2 / spe_sizes[i] # maybe SE
		print("S2 frac = %1.2f"%s2f[n+0])	
		print("S2 = %1.0f, S2ce = %1.0f, ratio = %1.3f"%(np.sum(s2[:,n]), np.sum(aa[:,n+0]), np.sum(aa[:,n+0])/np.sum(s2[:,n])) )
		print("S1 = %1.1f, S1 echo = %1.1f"%(np.sum(s1[:,n+0]),np.sum(ee[:,n+0])))
# 		if (np.sum(aa[:,n+0],axis=0)<1e4):#((n<10) or (n>210)):
		if 0:
			pl.figure(9);pl.clf()
			for i in range(0,32):
				pl.plot(ch_data_mV[i,n+0,:],'-',linewidth=0.5)			
			pl.plot(i_s2,0,'ko',markerfacecolor='None')
			pl.plot(i_s2+i_w,0,'ks',markerfacecolor='None')
			pl.plot(i_s1,0,'ko',markerfacecolor='None')
			pl.plot(i_s1+i_w,0,'ks',markerfacecolor='None')
			pl.plot(np.array([1,1])*(i_ce-i_w),np.array([0,5]),'k-')
			pl.plot(np.array([1,1])*(i_ce+i_w),np.array([0,5]),'k-')
			pl.plot(np.array([1,1])*(i_ge-i_w),np.array([0,10]),'k-')
			pl.plot(np.array([1,1])*(i_ge+i_w),np.array([0,10]),'k-')
			pl.plot(np.array([1,1])*(i_ee-i_w),np.array([0,1]),'k-')
			pl.plot(np.array([1,1])*(i_ee+i_w),np.array([0,1]),'k-')
			pl.xlabel(r'$\mu$s');pl.ylabel('mV');pl.title('event %d'%n);
			#pl.ylim([-1,5]);pl.xlim([1.5e3,9e3])

			pl.show();pl.pause(0.1)	
			input('press any key to continue...')

		# check for pileup / problems
# 		a = np.array([np.sum(aa[:,n+0]),np.sum(aa[:,n+1]),np.sum(aa[:,n+2]),np.sum(aa[:,n+3])])
# 		#print(a)
# 		if (a[0]>a[1]):
# 			print("processed: n==%d"%n)
# 		if (a[1]>a[0]):
# 			print('detected: pileup (or pileup trigger?) on n==%d'%(n+1))
# 			aa[:,n+0]=0;aa[:,n+1]=0;aa[:,n+2]=0;aa[:,n+3]=0;
# 			#break
# 		if (a[2]>a[1]):
# 			print('detected: pileup (or pileup trigger?) on n==%d'%(n+2))
# 			aa[:,n+0]=0;aa[:,n+1]=0;aa[:,n+2]=0;aa[:,n+3]=0;
# 			#break
# 		if (a[3]>a[2]):
# 			print('detected: pileup (or pileup trigger?) on n==%d'%(n+3))
# 			aa[:,n+0]=0;aa[:,n+1]=0;aa[:,n+2]=0;aa[:,n+3]=0;
			#break
# 		input('press any key to continue...')
	print('saving events %d:%d'%(n_start,n_stop))
# 	aa_file = compressed_file[0:-29]+"aa_"+("%02d"%ff)+".npz"		
	aa_file = aa_dir+compressed_file[-9:] # write aa to npz
	np.savez(aa_file,aa,ee,s1,is1,is2,s2f,s2,ge)
	ff = ff + 1
	
	
	if 0:
		pl.figure(2);pl.clf();
		pl.subplot(1,2,1)
		pl.hist(y0[blr],beans,histtype='step')
		pl.hist(y1[:],beans,histtype='step',label='y1')
		pl.plot(np.ones(2)*0,np.array([1,8e3]),'k--')
		pl.plot(np.ones(2)*wm1,np.array([1,8e3]),'k-',label='pbz')
		pl.legend()
		pl.yscale('log')
		pl.title('evt %d ch %d'%(n,i))
	
	
	
		pl.subplot(1,2,2)
		pl.plot(y1,'.',markersize=3,label='as-is');
		pl.plot(y1z,'k.',markersize=3,label='pbz');
		pl.plot(np.array([0,5e4]),np.ones(2)*bstd*4,'k--',label='4$\sigma$')
		pl.plot(np.array([0,5e4]),-np.ones(2)*bstd*4,'k--')
		y1[y1>bstd*3]
		pl.ylim([-0.3,0.5])
		pl.legend()
		pl.show();pl.pause(0.1)	
		input('press any key to continue...')


