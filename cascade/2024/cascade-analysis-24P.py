import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob
from scipy import stats as st
#from scipy import signal

# 2024-06-13 - pfs cleaned up and improved algorithm
# 2024-06-14 - pfs decreased "ad-hoc threshold" from 1000 to 800 based on 
#					0.5us LED data, should re-run all analysis even though 
#					I expect no change
# 2024-09-24 - pfs forked off from 24H
#
# 2024-11-06 - pfs forked off from 24L, now have only a single PMT lower instead of 16 sipms
# 2024-12-17 - pfs forked off from 24O

plotz=0


def linfit(X,Y):
	# following "Error Analysis" by Taylor 2nd Ed. p181 -pfs
	# returns coefficients a and b for the equation Y = a + bX
	N = len(X)
	delta = N*np.sum(X**2) - (np.sum(X))**2
	a = ( np.sum(X**2) * np.sum(Y) - np.sum(X) * sum(X*Y) ) / delta
	b = ( N*sum(X*Y) - sum(X) * sum(Y) ) / delta
	sigma_y = np.sqrt( (1/(N-2)) * sum((Y-a-b*X)**2) )
	sigma_a = sigma_y * np.sqrt(np.sum(X**2) / delta)
	sigma_b = sigma_y * np.sqrt(N/delta)
	
	return (a, b, sigma_a, sigma_b)

#https://drive.google.com/drive/folders/1DjYev4OSSBny3FIhqXFQoQJvNfHxt7ow
### sandboxing #(socket.gethostname()[0]=='b') :



if 1: # 24P single PMT bottom, cathode alphas (?)
	data_dir = '/Users/peter/Public/data/20241219-121719/' # -85C xenon gas

spe_dir = "/Users/peter/Dropbox/GitHub/crystaLiZe/cascade/50V_3-6-2024-plus-pmt.txt"

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
		print("loaded header file")
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
	spe_sizes = spe_sizes[0:17]
	
	
ff = 0	
compressed_file_list = glob.glob(data_dir+"./compressed_filtered_data/c*.npz")
compressed_file_list = sorted( compressed_file_list )

sipms = np.arange(0,17) #
nch = sipms.shape[0]


print("looking in: %s"%compressed_file_list[0][0:41])
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
		event_window = 100 #get_event_window(data_dir)
		cascade = np.array([0.5,1.0,5])
	
		wsize = int(500 * event_window)  # samples per waveform # 12500 for 25 us
		tscale = (8.0/4096.0)     # = 0.002 µs/sample, time scale
		t = np.arange(0,wsize)*tscale # time array for plotting
		n_sipms = 17 # number of SiPMs
		n_channels = n_sipms # no sum
		block_size = int(1500*15/event_window) # number of events per compressed file
	
		n_tot_samp_per_ch = int( (ch_data_adcc.size)/n_sipms )
		n_events_b = int((ch_data_adcc.size)/(n_sipms*wsize)) # n events per compressed file (same as block_size)

		# Convert from ADCC to mV
		ch_data_mV = vscale*np.reshape(ch_data_adcc, (n_channels,n_events_b,wsize))
		ch_data_mV[16,:,:] = -ch_data_mV[16,:,:] # invert PMT

		ds_pmt = 5 #spe_sizes[16] / np.mean(spe_sizes[0:16])
		ch_data_mV[16,:,:] = ch_data_mV[16,:,:]/ds_pmt # downsize the PMT to make it SiPM size


	# buffer some memory for the RQs
	nSE = 20 # max to find
	nspe = 10 # max to find
	aa = np.zeros([nch, ch_data_mV.shape[1]]) # areas of delay_window_0[0:20us], delay_window_1, delay_window_2, delay_window_3
	s1 = np.zeros([nch, ch_data_mV.shape[1]]) # areas of s1, 0, 0, 0
	ee = np.zeros([nch, ch_data_mV.shape[1],nSE]) # areas of SE, 0, 0, 0
	is1 = np.zeros([ch_data_mV.shape[1]]) # i_s1, 0,0,0
	is2 = np.zeros([ch_data_mV.shape[1]]) # i_s2, 0,0,0
	s2cf = np.zeros([ch_data_mV.shape[1]]) # s2_frac, 0,0,0
	ss = np.zeros([nch, ch_data_mV.shape[1],nspe])  # areas of spe

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

	# rough check on cascade ordering/status
	a = np.array([ np.sum(np.sum(ch_data_mV[:,0,:],axis=0)) , np.sum(np.sum(ch_data_mV[:,1,:],axis=0)) , np.sum(np.sum(ch_data_mV[:,2,:],axis=0)) , np.sum(np.sum(ch_data_mV[:,3,:],axis=0)) ])
	n_start = int(np.argmax(a))
	n_stop = int(np.floor((ch_data_mV.shape[1]-n_start)/4)*4)
	print('rough check indicates start with event %d (and so skip last %d events)'%(n_start,(ch_data_mV.shape[1]-n_stop)))
# 	for n in np.arange(n_start,ch_data_mV.shape[1]-n_oops,4):
	for n in np.arange(n_start,n_stop,4):
		i_max = int(np.mean(np.argmax(ch_data_mV[sipms,n+0,:],axis=1)))
		s=np.zeros(nch,dtype=np.int16)
		for i in range(0,nch): # refined edge-finding
			while (ch_data_mV[i,n+0,i_max-s[i]]>0.05): 
				s[i]=s[i]+1
				if s[i]>i_max:
					break
		i_s2 = int(i_max-np.mean(s))			
		i_s1 = np.argmax(np.sum(ch_data_mV[0:nch,n+0,0:i_s2-150],axis=0),axis=0)-40 # 40 is the rise time
		i_ge = i_s2 + 1150 # by eye, index of max of gate echo
		i_ce = i_s2 + 3250 # by eye, index of max of cathode echo
		i_ee = i_s1 + 1125 # by eye
		i_w = 375 # was 325 prior to June 13

		is1[n+0] = i_s1
		is2[n+0] = i_s2
# 		s2cf[n+0] = np.sum(np.sum(ch_data_mV[:,n+0,i_s2:i_s2+1000])) / np.sum(np.sum(ch_data_mV[:,n+0,:])) # approx S2 center fraction
		
		# these are all in mV
		y0s = np.zeros(ch_data_mV.shape[2])
		y1s = np.zeros(ch_data_mV.shape[2])
		y2s = np.zeros(ch_data_mV.shape[2])
		y3s = np.zeros(ch_data_mV.shape[2])

		# this is in spe
		y0spe = np.zeros(ch_data_mV.shape[2])
		y1spe = np.zeros(ch_data_mV.shape[2])
		y2spe = np.zeros(ch_data_mV.shape[2])
		y3spe = np.zeros(ch_data_mV.shape[2])
		
		s8 = 0.1 # parameter for baseline zero-suppression. was previously using 4*bstd0
		y0z = np.zeros((ch_data_mV.shape[0],ch_data_mV.shape[2]))
		y1z = np.zeros((ch_data_mV.shape[0],ch_data_mV.shape[2]))
		y2z = np.zeros((ch_data_mV.shape[0],ch_data_mV.shape[2]))
		y3z = np.zeros((ch_data_mV.shape[0],ch_data_mV.shape[2]))
		# last half of 1st trigger
		y0p5z = np.zeros((ch_data_mV.shape[0],int(ch_data_mV.shape[2]/2)))

		print('\nEvent %d'%n)		
		for i in range(0,nch): # now go an integrate individual channel areas
			y0 = ch_data_mV[i,n+0,:]
			y1 = ch_data_mV[i,n+1,:]
			y2 = ch_data_mV[i,n+2,:]
			y3 = ch_data_mV[i,n+3,:]
			xx = np.arange(0,y0.shape[0],1)
			
			beans = np.arange(-0.2,0.2,0.001)
			binc = (beans[0:-1] + beans[1:])/2
			blr = np.arange(0,2000) # only for initiating event in cascade
			bstd0 = np.std(y0[blr]) # get the std of the triggering event baseline
	
			# need weighted mean for baseline -- or do I? see below, try fitting baseline
			[counts0,bine] = np.histogram(y0[blr],beans)
			try:
				wm0 = np.average(binc,axis=0,weights=counts0)
			except:
				wm0 = 0
			y0z[i,:] = y0-wm0 # apply additional baseline-correction to data	
			y0s[y0z[i,:]>s8] = y0s[y0z[i,:]>s8] + y0z[i,(y0z[i,:]>s8)]
			y0spe[y0z[i,:]>s8] = y0spe[y0z[i,:]>s8] + y0z[i,(y0z[i,:]>s8)]*2/spe_sizes[i]

			[counts1,bine] = np.histogram(y1,beans)
			(a, b, sigma_a, sigma_b) = linfit(xx,y1) # baseline often sloping, so fit it
			y1z[i,:] = y1-(a+b*xx)
			y1s[y1z[i,:]>s8] = y1s[y1z[i,:]>s8] + y1z[i,(y1z[i,:]>s8)]
			y1spe[y1z[i,:]>s8] = y1spe[y1z[i,:]>s8] + y1z[i,(y1z[i,:]>s8)]*2/spe_sizes[i] # also want the 1st cascade summed waveform

			[counts2,bine] = np.histogram(y2,beans)
			(a, b, sigma_a, sigma_b) = linfit(xx,y2) # baseline often sloping, so fit it
			y2z[i,:] = y2-(a+b*xx)
			y2s[y2z[i,:]>s8] = y2s[y2z[i,:]>s8] + y2z[i,(y2z[i,:]>s8)]
			y2spe[y2z[i,:]>s8] = y2spe[y2z[i,:]>s8] + y2z[i,(y2z[i,:]>s8)]*2/spe_sizes[i] # also want the 1st cascade summed waveform
			
			[counts3,bine] = np.histogram(y3,beans)
			(a, b, sigma_a, sigma_b) = linfit(xx,y3) # baseline often sloping, so fit it
			y3z[i,:] = y3-(a+b*xx)
			y3s[y3z[i,:]>s8] = y3s[y3z[i,:]>s8] + y3z[i,(y3z[i,:]>s8)]
			y3spe[y3z[i,:]>s8] = y3spe[y3z[i,:]>s8] + y3z[i,(y3z[i,:]>s8)]*2/spe_sizes[i] # also want the 1st cascade summed waveform

			# add special processing for last 50% of first trigger
			(a, b, sigma_a, sigma_b) = linfit(xx[25000:50000],y0[25000:50000]) # baseline often sloping, so fit it
			y0p5z[i,:] = y0[25000:50000]-(a+b*xx[25000:50000])
		
			# counting photoelectrons in aggregate
			aa[i,n+0] = np.sum(y0p5z[i,y0p5z[i,:]>s8]) * 2 / spe_sizes[i] 
			aa[i,n+1] = np.sum(y1z[i,y1z[i,:]>s8]) * 2 / spe_sizes[i] # factor x2 for 2 ns samples
			aa[i,n+2] = np.sum(y2z[i,y2z[i,:]>s8]) * 2 / spe_sizes[i] # factor x2 for 2 ns samples
			aa[i,n+3] = np.sum(y3z[i,y3z[i,:]>s8]) * 2 / spe_sizes[i] # factor x2 for 2 ns samples
			
			# no threshold here since all channels are enormous
			s1[i,n+0] = np.sum(y0z[i,0:10000]) * 2 / spe_sizes[i] 
# 			s1[i,n+0] = np.sum(y0z[:,i_s1:i_s1+i_w]) * 2 / spe_sizes[i] # rough, but should be close

			if i_s2>(i_ee+i_w) and i_s1<(i_ee-i_w):
				ee[i,n+0] = np.sum(y0z[i,i_ee-i_w:i_ee+i_w]) * 2 / spe_sizes[i] # maybe SE in S1 echo
							
		p_area_ch = np.sum(y0z[:,i_s2:i_s2+1000],axis=1)
# 		(bot_x, bot_y, top_x, top_y) = GetCentroids(p_area_ch)
# 		s2cf[n+0] = np.sqrt(bot_x**2 + bot_y**2)

# top sipm map
#
# 00 01 04 05
# 03 02 07 06
# 12 13 08 09
# 15 14 11 10
#

# 		print("s2 radial = %1.2f"%s2cf[n+0])	
		print("s1 = %1.2f"% np.sum(s1[:,n+0],axis=0))	
		print("cascade sums:%1.0f,%1.0f,%1.0f,%1.0f"%(np.sum(aa[:,n+0],axis=0),np.sum(aa[:,n+1],axis=0),np.sum(aa[:,n+2],axis=0),np.sum(aa[:,n+3],axis=0) ))
		if plotz: # 
			pl.figure(8);pl.clf()
	
			pl.subplot(1,4,1)
# 			pl.plot(t,y0s+1,'-',linewidth=0.5,color='grey')
			for i in range(0,nch-1):
				pl.plot(t,y0z[i,:],'-',linewidth=0.5)
# 				pl.plot(t[25000:],2+y0p5z[i,:],'-',linewidth=0.5)
			pl.plot(t,1+y0z[16,:],'k-',linewidth=0.5)
			
			pl.plot(t[i_s1],1,'ko',markerfacecolor='None')
			pl.plot(t[i_s2],1,'ko',markerfacecolor='None')
			pl.xlabel(r'$\mu$s');pl.ylabel('mV');pl.title('event %d'%n);
# 			pl.xlim([4,20])
			pl.ylim([-10,500])
			pl.ylim([-1,5])
	
			pl.subplot(1,4,2)
			for i in range(0,nch-1):
				pl.plot(t,y1z[i,:],'-',linewidth=0.5);
			try:
				pl.plot(t[argmax1[1:]],-0.2*np.ones(t[argmax1[1:]].shape),'k^')
			except:
				print('')
			pl.plot(t,1+y1z[16,:],'k-',linewidth=0.5);
			pl.xlabel(r'$\mu$s');pl.ylabel('mV');#pl.title('%1.1f ms delayed'%cascade[0]);
			pl.ylim([-1,5])
	
			pl.subplot(1,4,3)
			for i in range(0,nch-1):
				pl.plot(t,y2z[i,:],'-',linewidth=0.5);
			try:
				pl.plot(t[argmax2[1:]],-0.2*np.ones(t[argmax2[1:]].shape),'k^')
			except:
				print('')
			pl.plot(t,1+y2z[16,:],'k-',linewidth=0.5);
			pl.xlabel(r'$\mu$s');pl.ylabel('mV');#pl.title('%1.1f ms delayed'%cascade[1]);
			pl.ylim([-1,5])#pl.ylim([-0.2,3])
	
			pl.subplot(1,4,4)
			for i in range(0,nch-1):
				pl.plot(t,y3z[i,:],'-',linewidth=0.5);
			try:
				pl.plot(t[argmax3[1:]],-0.2*np.ones(t[argmax3[1:]].shape),'k^')
			except:
				print('')
			pl.plot(t,1+y3z[16,:],'k-',linewidth=0.5);
			pl.xlabel(r'$\mu$s');pl.ylabel('mV');#pl.title('%1.1f ms delayed'%cascade[2]);
			pl.ylim([-1,5])
	
			pl.show();pl.pause(0.1)	


		# get single phd areas for PMT
		def getss(wf,nspe):
			spes = np.zeros(nspe) # same size as ss defined above
			for it in range(0,spes.shape[0]): # find several pulses
				mxv = np.max(wf)
				amx = np.argmax(wf)
				kr = amx; kl = amx
				thr = 0.1
				while (wf[kr]>thr):
					kr+=1
					if kr>=wf.shape[0]:
						break
				while (wf[kl]>thr):
					kl-=1
					if kl<=0:
						break
				spes[it] = np.sum(wf[kl:kr])
				wf[kl:kr] = 0
			return spes

		for ch in range(0,nch-1): #ch = 16 # to save spe data (just one channel? jeez...)
			if plotz:
				pl.figure(9);pl.clf()	
				pl.plot(t[25000:],y0p5z[ch,:],'-',linewidth=0.5)
				pl.ylim([-0.1,0.5])	
				pl.ylim([-0.1,0.5])	
				pl.xlabel(r'$\mu$s');pl.ylabel('mV');pl.title('event %d'%n);
				pl.show();pl.pause(0.1)	
			ss[ch,n+0,:] = getss(y0p5z[ch,:],nspe)*2 # mV ns
			ss[ch,n+1,:] = getss(y1z[ch,:],nspe)*2
			ss[ch,n+2,:] = getss(y2z[ch,:],nspe)*2
			ss[ch,n+3,:] = getss(y3z[ch,:],nspe)*2
# 			print("spes in 50-100 us of initiating trigger:")
# 			print(ss[:,n+0])

		if plotz:
			input('press any key to continue...')


	print('saving events %d:%d'%(n_start,n+4))
	aa = aa[:,n_start:n+4] # do I have the indexing right??
	ee = ee[:,n_start:n+4] 
	s1 = s1[:,n_start:n+4] 
	is1 = is1[n_start:n+4] 
	is2 = is2[n_start:n+4] 
	s2cf = s2cf[n_start:n+4]
	ss = ss[:,n_start:n+4]


# 	aa_file = compressed_file[0:-29]+"aa_"+("%02d"%ff)+".npz"		
	aa_file = aa_dir+compressed_file[-9:] # write aa to npz
	np.savez(aa_file,aa,ee,s1,is1,is2,s2cf,ss)
	ff = ff + 1
	






