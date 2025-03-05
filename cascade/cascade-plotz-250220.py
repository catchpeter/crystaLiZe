import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob



pmt_ch = 16
ds_pmt = 5 # pmt signal was scaled down for the analysis, put it back
if 1: # PTFE, Xe liquid, alphas from cathode pointing down I think, PMT bottom, compact cylinder, not really a TPC
	data_folders = np.array(['20250220-094512']) 
	data_folders = np.array(['20250220-131431']) 

	data_folders = np.array(['20250226-202029']) # liquid Xe, new PMT

	data_folders = np.array(['20250227-132142']) # 300K no Xe

	colorz = np.array(['gray','steelblue','olivedrab','goldenrod','firebrick','sienna'])
	labl = np.array(['','','',''])

dtt = 0.01
tt = np.arange(dtt,1.5,dtt)
fitdp = np.array([1.0e-4*tt**-1,1.2e-4*tt**-1])



af = np.zeros((5,data_folders.shape[0]))

for ii in range(0,1):#data_folders.shape[0]):
	data_dir = '/Users/peter/Public/data/'+data_folders[ii]+'/'
	aa_file_list = glob.glob(data_dir+"./aa/*v2.npz")
	print('looking in: %s'%data_folders[ii])
	print('found %d files'%len(aa_file_list))
	h_file = np.load(data_dir+"/compressed_filtered_data/headers.npz")
	h_array = h_file["arr_0"]
	h_n_events = int(np.floor(h_array.size/8))

	# check data_size
	check_data = np.load(aa_file_list[0])
	nch = check_data["arr_2"].shape[0]
	nss = check_data["arr_3"].shape[2]
	aa = np.zeros([nch,h_n_events])
	ss = np.zeros([nch,h_n_events,nss])
	s1 = np.zeros([nch,h_n_events])
	s1ap = np.zeros([nch,h_n_events])
	aa_last = 0
	for aa_file in aa_file_list:
		#print('load file %s'%aa_file)
	#	np.savez(aa_file,aa,ee,s1,is1,is2,s2f)
		with np.load(aa_file) as data:
			a = data["arr_0"].shape[1]
			aa[:,aa_last:(aa_last+a)] = data["arr_0"]
			s1[:,aa_last:(aa_last+a)] = data["arr_1"]
			s1ap[:,aa_last:(aa_last+a)] = data["arr_2"]
			ss[:,aa_last:(aa_last+a)] = data["arr_3"]
		aa_last = aa_last + a

	ei = (h_n_events-1) # end index
	ei = int(np.floor(ei/4)*4)

	print('multiplying back the PMT down-size factor x%d'%ds_pmt)
	ss[pmt_ch,:,:] = ss[pmt_ch,:,:]*ds_pmt
	aa[pmt_ch,:] = aa[pmt_ch,:]*ds_pmt
	s1[pmt_ch,:] = s1[pmt_ch,:]*ds_pmt
	s1ap[pmt_ch,:] = s1ap[pmt_ch,:]*ds_pmt

	# get gains
	gains = np.zeros(nch)
	for ch in range (0,nch):
		db=1
		beans = np.arange(0,150,db)
		beanc = (beans[1:]+beans[0:-1])/2
		[cts0,beans] = np.histogram(ss[ch,np.arange(0,ei,4),:],beans); cts0[0]=0 
		[cts1,beans] = np.histogram(ss[ch,np.arange(1,ei,4),:],beans); cts1[0]=0 
		if (ch<pmt_ch):
			gains[ch] = np.average(beanc[7:33],axis=0,weights=cts0[7:33]) # 7 picked to stay above noise bkg
		else:
			gains[ch] = np.average(beanc[7:],axis=0,weights=cts1[7:]) # 7 picked to stay above noise bkg
		
		if 1:
			pl.figure(6);pl.clf();
			pl.errorbar(beanc,cts0,yerr=np.sqrt(cts0),xerr=db/2,fmt='k.')

			pl.errorbar(beanc,cts1,yerr=np.sqrt(cts1),xerr=db/2,fmt='r.')
			
			[cts,beans] = np.histogram(ss[ch,np.arange(2,ei,4),:],beans); cts[0]=0 
			pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='b.')
			[cts,beans] = np.histogram(ss[ch,np.arange(3,ei,4),:],beans); cts[0]=0 
			pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='g.')
	
			pl.plot(np.ones(2)*gains[ch],np.array([1,max(cts0)]),'m:')

			pl.yscale('log')
			pl.title('ch %d, gain %2.1f'%(ch,gains[ch]))
			pl.show();pl.pause(0.1)	
# 			input('press any key')
			
		if (ch==pmt_ch) & (data_folders[0][-6:]=='103023'):
			print('assume PMT gain measurement failed, and use previous value ==32')
			gains[pmt_ch] = 32
		s1[ch,:] = s1[ch,:] / gains[ch]
		s1ap[ch,:] = s1ap[ch,:] / gains[ch]
		aa[ch,:] = aa[ch,:] / gains[ch]

#		input('paused...')


	asum = np.sum(aa,axis=0)
	a0t = np.sum( aa[0:nch-1,np.arange(0,ei,4)] ,axis=0 )
	a1t = np.sum( aa[0:nch-1,np.arange(1,ei,4)] ,axis=0 )
	a2t = np.sum( aa[0:nch-1,np.arange(2,ei,4)] ,axis=0 )
	a3t = np.sum( aa[0:nch-1,np.arange(3,ei,4)] ,axis=0 )

	a0b = ( aa[nch-1,np.arange(0,ei,4)])
	a1b = ( aa[nch-1,np.arange(1,ei,4)])
	a2b = ( aa[nch-1,np.arange(2,ei,4)])
	a3b = ( aa[nch-1,np.arange(3,ei,4)])

	
	print('NOTE: code assumes last channel in the array is the PMT')
	s1bot = (s1[pmt_ch,np.arange(0,ei,4)])
	s1bot_ap = (s1ap[pmt_ch,np.arange(0,ei,4)])
	s1top = np.sum(s1[0:pmt_ch,np.arange(0,ei,4)],axis=0)
	s1top_ap = np.sum(s1ap[0:pmt_ch,np.arange(0,ei,4)],axis=0)
	s1top = s1top + s1top_ap
	s10 = s1top + s1bot
	s1tba = (s1top-s1bot)/(s1top+s1bot)
	
	if 1: ## quality cuts
		pl.figure(10);pl.clf();
		db=50
		beans = np.arange(0,15e3,db)
		beanc = (beans[1:]+beans[0:-1])/2
		[cts,beans] = np.histogram(s1top,beans); cts[0]=0
		pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='c+',label='sipm top')
		[cts,beans] = np.histogram(s1bot,beans); cts[0]=0
		pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='b.',label='pmt bot')
		pl.title(data_folders[ii])
		pl.xlabel('phd')
		pl.legend()
		pl.minorticks_on()
		
		if 0:
			pl.figure(11);pl.clf();
			pl.plot(a0t,'o')
			pl.plot(a1t,'*')
			pl.plot(a2t,'^')
			pl.plot(a3t,'+')
	
	
		pl.figure(12);pl.clf();
		pl.plot((s1top+s1bot),s1tba,'k.')
# 		pl.plot((s1top),s1tba,'c+')
		pl.xlim([-100,12e3])
		pl.ylim([-1,1])

		# not a cut, rather a calibration of afterpulse as proxy for (saturated) S1
		pl.figure(13);pl.clf();
		pl.plot(s1bot,s1bot_ap,'o',color='powderblue',markersize=2)
		db=2
		beans = np.arange(0,3500,db)
		beanc = (beans[1:]+beans[0:-1])/2
		[cts,beans] = np.histogram(s1bot,beans); cts[0]=0
		pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='k.')
		
		[cts,beans] = np.histogram(s1bot_ap,beans); cts[0]=0
		pl.errorbar(cts,beanc,xerr=np.sqrt(cts),yerr=db/2,fmt='.',color='grey')
		
		if 0:
			pl.figure(14);pl.clf();
			db=0.1
			beans = np.arange(0,10,db)
			beanc = (beans[1:]+beans[0:-1])/2
			[cts,beans] = np.histogram(s1bot/s1bot_ap,beans); cts[0]=0
			pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='k.')
	
		if 0:
			pl.figure(15);pl.clf();
			pl.plot(s1bot,s1top,'o',color='gray',markersize=1)


	###
	cut = (s10>1000) & (s10<10000) \
		& (a1t<a0t) & (a2t<a0t) & (a3t<a0t) & (a1b<a0b) & (a2b<a0b) & (a3b<a0b)

	if (data_folders[0][-15:-7]=='20250227'): # then looking at no-xenon data
		cut = (s10>300) & (s10<1000)
	
# 	cut = (s10>1000) & (s10<5000) \
# 		& (a0t<50) & (a1t<20) & (a2t<20) & (a3t<20) \
# 		& (s1tba>-0.25) & (s1tba<0.25) \
		

	nz=np.nonzero(cut); nz = nz[0]
	N = np.sum(cut)
	print('cut keeps %d events'%N)
	###

			



	### aggregate the data
	

# 	ct = np.array([0.01, 0.075, 0.5+0.05, 1.0+0.05, 5.0+0.05]) # trigger cascade (set by hardware) -- used for ever all earlier data
	ct = np.array([0.01, 0.075, 0.2+0.05, 0.5+0.05, 1.0+0.05]) # trigger cascade (set by hardware) -- afternoon Feb 20 +
	if (data_folders[0][-6:]=='094512'):
		ct = np.array([0.01, 0.075, 0.3+0.05, 0.5+0.05, 1.0+0.05]) # trigger cascade (set by hardware) -- briefly used before noon on Feb 20

	# Rce = 0.0063 # small-s2 limit of ratio S2ce/S2 -- should triple check
	### define the number of detected photons. subtract the number of phd identified as single e-
	dpt = np.array([ np.sum(s1top[cut])/N , np.sum(a0t[cut])*2/N , np.sum(a1t[cut])/N , np.sum(a2t[cut])/N , np.sum(a3t[cut])/N ])
	print('\n*** NOTE: accounting for PMT ADC saturation factor 5/3, obtained from afterpulse size')
	dpb = np.array([ np.sum(s1bot[cut])/N*(5/3) , np.sum(a0b[cut])*2/N , np.sum(a1b[cut])/N , np.sum(a2b[cut])/N , np.sum(a3b[cut])/N ])
# 	af[:,ii] = detp#/triggerS1
	
	if 1:
		pbw = 0.1 # plot bin width, fixed to cascade event window!
		pl.figure(8);#pl.clf()
		pl.errorbar(ct,dpt,yerr=np.sqrt(dpt*N)/N,xerr=np.array([0.0025,0.025,0.05,0.05,0.05]),fmt='o',color='grey',markersize=5,markerfacecolor='white',label='sipm')
		pl.errorbar(ct,dpb,yerr=np.sqrt(dpb*N)/N,xerr=np.array([0.0025,0.025,0.05,0.05,0.05]),fmt='o',color='powderblue',markersize=5,markerfacecolor='white',label='pmt')

		# photons	
		#pl.plot(tt,1.2e-4*tt**-1.3,'k-',linewidth=0.5) # PTFE
		#pl.plot(tt,0.8e-4*tt**-1.3,'k:',linewidth=0.5) # Aluminum
		pl.plot(tt,np.sum(s1top[cut])/N*1.2e-4*tt**-1.3,'k-',linewidth=0.5,label='2024 PTFE') # PTFE
		
		if (int(data_folders[0][-15:-7]) == 20250220): # original PMT		
			pl.plot(tt,1.7*tt**-1.0,'k--',linewidth=0.5) 
			pl.plot(tt,0.8*tt**-1.5,'--',color='powderblue',linewidth=0.5) 
		if (int(data_folders[0][-15:-7]) == 20250226): # new PMT		
			pl.plot(tt,1.8*tt**-1.0,'k--',linewidth=0.5) 
			pl.plot(tt,0.6*tt**-1.5,'--',color='powderblue',linewidth=1) 

		if (int(data_folders[0][-15:-7]) == 20250227): # new PMT, no Xe		
			pl.plot(tt,1.8*tt**-1.0,'k--',linewidth=0.5) 
			pl.plot(tt,0.15*0.6*tt**-1.5,'--',color='powderblue',linewidth=1) 
		
# 		pl.plot(tt,fitdp[ii],'-',linewidth=0.5,color=colorz[ii],label=pwrlabl[ii])
				
# 		pl.plot(np.array([1e-3,1e-1]),np.array([1,1]),'k:',label='progenitor window')	
		pl.xlabel('time (ms)')
		pl.ylabel('photons / 0.1 ms')
# 		pl.ylim([1e-5,2])
		pl.xscale('log')
		pl.yscale('log')
		pl.legend()
		pl.title(data_dir[-16:-1])
		pl.ylim([0.1,1e4])

	if 0:
		for ch in range (0,nch):
			pl.figure(7);pl.clf();
			pl.loglog(s1[ch,np.arange(0,ei,4)],aa[ch,np.arange(1,ei,4)],'o',color='gray',markersize=1)
			pl.loglog(s1[ch,np.arange(0,ei,4)],aa[nch-ch-1,np.arange(1,ei,4)],'o',color='blue',markersize=1)
			pl.xlim([0.3,3e3])
			pl.ylim([0.3,1e2])
			pl.title('ch %1.0f' % ch)
			pl.show();pl.pause(0.1)	
# 			input('press any key')


if 0:
	print('cascade times:')
	print(dpt)
	print('detected photons:')
	print(detp)
# 	print('yerr:')
# 	print(np.sqrt(detp*N)/N/triggerS1)
	
	if 0:
		Ef = np.array([0,1,2,3,4,5]) # electric field
		pl.figure(7);pl.clf()
		pl.plot(Ef,af[4,:],'ro')


