import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob
from scipy import stats as st

### sandboxing
if (socket.gethostname()[0]=='b') :
	data_dir = '/Users/peter/Public/data/20240314-190415/' # first good cascade data
#	data_dir = '/Users/peter/Public/data/20240319-150343/' # after open to air
	spe_dir = "/Users/peter/Desktop/crops/50V_3-6-2024.txt"
else:
	data_dir = '/media/xaber/extradrive1/crystalize_data/data-202403/20240314/20240314-190415/'
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
	
	if 1: # specific to 3-fold cascade trigger
		shape1 = 4 # number of cascades in trigger
		shape0 = int(np.floor(ev_time_s.shape[0]/shape1))
		cascade_times_s = np.reshape(ev_time_s[0:shape0*shape1],( shape0 , shape1 ) )
	
	
	spe_sizes = np.loadtxt(spe_dir, dtype='float') # mV ns
	spe_sizes[spe_sizes==0] = 1e9

ff = 0	
compressed_file_list = glob.glob(data_dir+"./compressed_filtered_data/c*.npz")
compressed_file_list = sorted( compressed_file_list )
for compressed_file in compressed_file_list:
	if 0:
		print("loading: ... %s"%compressed_file_list[ff][42:])
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


	# buffer some memory for the RQs
	aa = np.zeros([32, ch_data_mV.shape[1]]) # areas of cathode_echo, delay_window_1, delay_window_2, delay_window_3
	s1 = np.zeros([32, ch_data_mV.shape[1]]) # areas of s1, 0, 0, 0
	ee = np.zeros([32, ch_data_mV.shape[1]]) # areas of SE, 0, 0, 0
	is1 = np.zeros([ch_data_mV.shape[1]]) # i_s1, 0,0,0
	is2 = np.zeros([ch_data_mV.shape[1]]) # i_s2, 0,0,0
	s2f = np.zeros([ch_data_mV.shape[1]]) # s2_frac, 0,0,0
	eec = np.zeros([ch_data_mV.shape[1]]) # s2_frac, 0,0,0

	print("analyzing compressed data...")
# 	if ff==0:
# 		n_start = 0 # start from 0 in ideal world, but offset if a cascade was split
# 	else:
# 		n_start = 4-n_oops # this is n_oops from the previous file
# 	# then we define n_oops for the current file...
# 	n_oops = int(4*(np.abs(np.floor((ch_data_mV.shape[1]-n_start)/4)-(ch_data_mV.shape[1]-n_start)/4)))
# 	print('file has %d events'%ch_data_mV.shape[1])
# 	print('will ignore first %d triggers'%n_start)
# 	print('will ignore last %d triggers'%(n_oops))

	sipms = np.concatenate((np.arange(0,27),np.arange(28,32))) # skip 27

	# rough check on cascade ordering/status
	a = np.array([ np.sum(np.sum(ch_data_mV[:,0,:],axis=0)) , np.sum(np.sum(ch_data_mV[:,1,:],axis=0)) , np.sum(np.sum(ch_data_mV[:,2,:],axis=0)) , np.sum(np.sum(ch_data_mV[:,3,:],axis=0)) ])
	n_start = int(np.argmax(a))
	n_stop = int(np.floor((ch_data_mV.shape[1]-n_start)/4)*4)
	print('rough check indicates start with event %d (and so skip last %d events)'%(n_start,(ch_data_mV.shape[1]-n_stop)))
# 	for n in np.arange(n_start,ch_data_mV.shape[1]-n_oops,4):
	for n in np.arange(n_start,n_stop,4):
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

		y0s = np.zeros(ch_data_mV.shape[2])
		y1s = np.zeros(ch_data_mV.shape[2])
		y2s = np.zeros(ch_data_mV.shape[2])
		y3s = np.zeros(ch_data_mV.shape[2])
		
		s8 = 0.1 # parameter for baseline zero-suppression. was previously using 4*bstd0
		y0z = np.zeros((ch_data_mV.shape[0],ch_data_mV.shape[2]))
		y1z = np.zeros((ch_data_mV.shape[0],ch_data_mV.shape[2]))
		y2z = np.zeros((ch_data_mV.shape[0],ch_data_mV.shape[2]))
		y3z = np.zeros((ch_data_mV.shape[0],ch_data_mV.shape[2]))

		print('\nEvent %d'%n)		
		for i in range(0,32): # now go an integrate individual channel areas
			y0 = ch_data_mV[i,n+0,:]
			y1 = ch_data_mV[i,n+1,:]
			y2 = ch_data_mV[i,n+2,:]
			y3 = ch_data_mV[i,n+3,:]

			beans = np.arange(-0.2,0.2,0.001)
			binc = (beans[0:-1] + beans[1:])/2
			blr = np.arange(0,2000) # only for initiating event in cascade
			bstd0 = np.std(y0[blr]) # get the std of the triggering event baseline
	
			# need weighted mean for baseline
			[counts0,bine] = np.histogram(y0[blr],beans)
			wm0 = np.average(binc,axis=0,weights=counts0)
			y0z[i,:] = y0-wm0 # apply additional baseline-correction to data	
			y0s[y0z[i,:]>s8] = y0s[y0z[i,:]>s8] + y0z[i,(y0z[i,:]>s8)]

			[counts1,bine] = np.histogram(y1,beans)
			wm1 = np.average(binc,axis=0,weights=counts1)
			y1z[i,:] = y1-wm1 # apply additional baseline-correction to data	
			y1s[y1z[i,:]>s8] = y1s[y1z[i,:]>s8] + y1z[i,(y1z[i,:]>s8)]

			[counts2,bine] = np.histogram(y2,beans)
			wm2 = np.average(binc,axis=0,weights=counts2)
			y2z[i,:] = y2-wm2 # apply additional baseline-correction to data	
# 			y2s[y2z>(4*bstd0)] = y2s[y2z>(4*bstd0)] + y2z[y2z>(4*bstd0)]
			y2s[y2z[i,:]>s8] = y2s[y2z[i,:]>s8] + y2z[i,(y2z[i,:]>s8)]
			
			[counts3,bine] = np.histogram(y3,beans)
			wm3 = np.average(binc,axis=0,weights=counts3)
			y3z[i,:] = y3-wm3 # apply additional baseline-correction to data	
			y3s[y3z[i,:]>s8] = y3s[y3z[i,:]>s8] + y3z[i,(y3z[i,:]>s8)]
			
			aa[i,n+0] = np.sum(y0z[i,i_ce-i_w:i_ce+i_w]) * 2 / spe_sizes[i] # rough, but should be close
# 			aa[i,n+2] = np.sum(y2z[y2z>s8]) * 2 / spe_sizes[i] # factor x2 for 2 ns samples
			aa[i,n+1] = np.sum(y1z[i,y1z[i,:]>s8]) * 2 / spe_sizes[i] # factor x2 for 2 ns samples
			aa[i,n+2] = np.sum(y2z[i,y2z[i,:]>s8]) * 2 / spe_sizes[i] # factor x2 for 2 ns samples
			aa[i,n+3] = np.sum(y3z[i,y3z[i,:]>s8]) * 2 / spe_sizes[i] # factor x2 for 2 ns samples
			
			s1[i,n+0] = np.sum(y0z[:,i_s1:i_s1+i_w]) * 2 / spe_sizes[i] # rough, but should be close
			if i_s2>(i_ee+i_w) and i_s1<(i_ee-i_w):
				ee[i,n+0] = np.sum(y0z[i,i_ee-i_w:i_ee+i_w]) * 2 / spe_sizes[i] # maybe SE in S1 echo

		# find SE in summed cascade waveforms
		it=0
		if (np.sum(aa[:,n+1],axis=0)<1000): # ad-hoc threshold, above which the waveform is too messy to try to find SE
			while (np.max(y1s)>2): # then likely SE (based on handscanning)
				for i in range(0,32):
					# this next line is not ideal as it co-adds all the SE areas (need to create a new variable if we want to collect info on multiple SE areas per trigger)
					ee[i,n+1] = ee[i,n+1] + np.sum(y1z[i,np.argmax(y1s)-i_w:np.argmax(y1s)+i_w]) * 2 / spe_sizes[i]
				y1s[np.argmax(y1s)-i_w:np.argmax(y1s)+i_w] = -y1s[np.argmax(y1s)-i_w:np.argmax(y1s)+i_w] # invert the region so we can look for more SE (a lazy trick)
				eec[n+1] = eec[n+1] + 1				
				print("..................................found SE (cascade 1) %2.1f"%np.sum(ee[:,n+1]))
				it=it+1
				if (it>10):
					break
		it=0
		if (np.sum(aa[:,n+2],axis=0)<1000): # ad-hoc threshold, above which the waveform is too messy to try to find SE
			while (np.max(y2s)>2): # then likely SE (based on handscanning)
				for i in range(0,32):
					# this next line is not ideal as it co-adds all the SE areas (need to create a new variable if we want to collect info on multiple SE areas per trigger)
					ee[i,n+2] = ee[i,n+2] + np.sum(y2z[i,np.argmax(y2s)-i_w:np.argmax(y2s)+i_w]) * 2 / spe_sizes[i]
				y2s[np.argmax(y2s)-i_w:np.argmax(y2s)+i_w] = -y2s[np.argmax(y2s)-i_w:np.argmax(y2s)+i_w] # invert the region so we can look for more SE (a lazy trick)
				eec[n+2] = eec[n+2] + 1
				print("..................................found SE (cascade 2) %2.1f"%np.sum(ee[:,n+2]))
				it=it+1
				if (it>10):
					break
		it=0
		if (np.sum(aa[:,n+3],axis=0)<1000): # ad-hoc threshold, above which the waveform is too messy to try to find SE
			while (np.max(y3s)>2): # then likely SE (based on handscanning)
				for i in range(0,32):
					# this next line is not ideal as it co-adds all the SE areas (need to create a new variable if we want to collect info on multiple SE areas per trigger)
					ee[i,n+3] = ee[i,n+3] + np.sum(y3z[i,np.argmax(y3s)-i_w:np.argmax(y3s)+i_w]) * 2 / spe_sizes[i]
				y3s[np.argmax(y3s)-i_w:np.argmax(y3s)+i_w] = -y3s[np.argmax(y3s)-i_w:np.argmax(y3s)+i_w] # invert the region so we can look for more SE (a lazy trick)
				eec[n+3] = eec[n+3] + 1
				print("..................................found SE (cascade 3) %2.1f"%np.sum(ee[:,n+3]))
				it=it+1
				if (it>10):
					break
			
		print("s2 frac = %1.2f"%s2f[n+0])	
		print("cascade sums:%1.0f,%1.0f,%1.0f,%1.0f"%(np.sum(aa[:,n+0],axis=0),np.sum(aa[:,n+1],axis=0),np.sum(aa[:,n+2],axis=0),np.sum(aa[:,n+3],axis=0) ))
		print("S1 = %1.1f, S1 echo = %1.1f"%(np.sum(s1[:,n+0]),np.sum(ee[:,n+0])))
#		if (np.sum(aa[:,n+0],axis=0)<2e3):#((n<10) or (n>210)):
		if 1:
			pl.figure(8);pl.clf()
	
			pl.subplot(1,4,1)
			pl.plot(t,y0s+1,'-',linewidth=0.5)
			for i in range(0,32):
				pl.plot(t,y0z[i,:],'-',linewidth=0.5)
			
			pl.plot(t[i_s1],1,'ko',markerfacecolor='None')
			pl.plot(t[i_s2],1,'ko',markerfacecolor='None')
			pl.plot(t[np.array([1,1])*(i_ce-i_w)],np.array([0,10]),'k-')
			pl.plot(t[np.array([1,1])*(i_ce+i_w)],np.array([0,10]),'k-')
			pl.plot(t[np.array([1,1])*(i_ge-i_w)],np.array([0,10]),'k-')
			pl.plot(t[np.array([1,1])*(i_ge+i_w)],np.array([0,10]),'k-')
			pl.plot(t[np.array([1,1])*(i_ee-i_w)],np.array([0,2]),'k-')
			pl.plot(t[np.array([1,1])*(i_ee+i_w)],np.array([0,2]),'k-')
			pl.xlabel(r'$\mu$s');pl.ylabel('mV');pl.title('event %d'%n);
# 			pl.xlim([4,12])
# 			pl.ylim([-1,5])
	
			pl.subplot(1,4,2)
# 			pl.plot(t,np.sum(ch_data_mV[:,n+1,:],axis=0),'-',linewidth=0.5);
			pl.plot(t,y1s+1,'-',linewidth=0.5)
			for i in range(0,32):
				pl.plot(t,y1z[i,:],'-',linewidth=0.5);
			pl.xlabel(r'$\mu$s');pl.ylabel('mV');pl.title('%1.1f ms delayed'%cascade[0]);
			pl.ylim([-1,5])
	
			pl.subplot(1,4,3)
# 			pl.plot(t,np.sum(ch_data_mV[:,n+2,:],axis=0),'-',linewidth=0.5);
			pl.plot(t,y2s+1,'-',linewidth=0.5)
			for i in range(0,32):
				pl.plot(t,y2z[i,:],'-',linewidth=0.5);
			pl.xlabel(r'$\mu$s');pl.ylabel('mV');pl.title('%1.1f ms delayed'%cascade[1]);
			pl.ylim([-1,5])#pl.ylim([-0.2,3])
	
			pl.subplot(1,4,4)
			pl.plot(t,y3s+1,'-',linewidth=0.5)
			for i in range(0,32):
				pl.plot(t,y3z[i,:],'-',linewidth=0.5);
# 				pl.plot(t,ch_data_mV[i,n+3,:],'-',linewidth=0.5);
			pl.xlabel(r'$\mu$s');pl.ylabel('mV');pl.title('%1.1f ms delayed'%cascade[2]);
			pl.ylim([-1,5])
	
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
				pl.ylim([-1,5]);pl.xlim([1.5e3,9e3])

			pl.show();pl.pause(0.1)	
			input('press any key to continue...')

		# check for pileup / problems
		a = np.array([np.sum(aa[:,n+0]),np.sum(aa[:,n+1]),np.sum(aa[:,n+2]),np.sum(aa[:,n+3])])
		#print(a)
		if (a[0]>a[1]):
			print("processed: n==%d"%n)
		if (a[1]>a[0]):
			print('detected: pileup (or pileup trigger?) on n==%d'%(n+1))
			aa[:,n+0]=0;aa[:,n+1]=0;aa[:,n+2]=0;aa[:,n+3]=0;
			#break
		if (a[2]>a[1]):
			print('detected: pileup (or pileup trigger?) on n==%d'%(n+2))
			aa[:,n+0]=0;aa[:,n+1]=0;aa[:,n+2]=0;aa[:,n+3]=0;
			#break
		if (a[3]>a[2]):
			print('detected: pileup (or pileup trigger?) on n==%d'%(n+3))
			aa[:,n+0]=0;aa[:,n+1]=0;aa[:,n+2]=0;aa[:,n+3]=0;
			#break
# 		input('press any key to continue...')
	print('saving events %d:%d'%(n_start,n+4))
	aa = aa[:,n_start:n+4] # do I have the indexing right??
	s1 = s1[:,n_start:n+4] 
	ee = ee[:,n_start:n+4] 
	is1 = is1[n_start:n+4] 
	is2 = is2[n_start:n+4] 
	s2f = s2f[n_start:n+4] 
# 	aa_file = compressed_file[0:-29]+"aa_"+("%02d"%ff)+".npz"		
	aa_file = aa_dir+compressed_file[-9:] # write aa to npz
	np.savez(aa_file,aa,ee,s1,is1,is2,s2f,eec)
	ff = ff + 1
	
	
	if 0:
		pl.figure(2);pl.clf();
		for i in range(0,32): # now go an integrate individual channel areas
			y2 = ch_data_mV[i,n+2,:]
			[counts2,bine] = np.histogram(y2,beans)
			wm2 = np.average(binc,axis=0,weights=counts2)
			pl.hist(y2-wm2,beans,histtype='step')
# 		pl.plot(np.ones(2)*0,np.array([1,8e3]),'k--')
# 		pl.plot(np.ones(2)*wm1,np.array([1,8e3]),'k-',label='pbz')
# 		pl.legend()
# 		pl.yscale('log')
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


