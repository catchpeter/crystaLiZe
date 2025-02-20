import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob



pmt_ch = 16
ds_pmt = 5 # pmt signal was scaled down for the analysis, put it back
if 1: # PTFE, Xe liquid, alphas from cathode pointing down I think, PMT bottom
#	data_folders = np.array(['20250114-211818'])
#	data_folders = np.array(['20250115-151222']); ds_pmt = 10 # pmt signal was scaled down for the analysis, put it back
	data_folders = np.array(['20250115-211732'])
	data_folders = np.array(['20250115-213907'])
	data_folders = np.array(['20250129-182353'])

	data_folders = np.array(['20250205-155211'])

	data_folders = np.array(['20250205-174037'])


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
		[cts,beans] = np.histogram(ss[ch,np.arange(0,ei,4),:],beans); cts[0]=0 
		if 1:
			pl.figure(6);pl.clf();
			pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='k.')

			[cts,beans] = np.histogram(ss[ch,np.arange(1,ei,4),:],beans); cts[0]=0 
			pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='r.')
			[cts,beans] = np.histogram(ss[ch,np.arange(2,ei,4),:],beans); cts[0]=0 
			pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='b.')
			[cts,beans] = np.histogram(ss[ch,np.arange(3,ei,4),:],beans); cts[0]=0 
			pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='g.')

			pl.yscale('log')
			pl.title(ch)
			pl.show();pl.pause(0.1)	

		gains[ch] = np.average(beanc[7:],axis=0,weights=cts[7:]) # 7 picked to stay above noise bkg
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
	s10 = s1top + s1bot
	s1tba = (s1top-s1bot)/(s1top+s1bot)
	
	if 1: ## quality cuts
		pl.figure(10);pl.clf();
		db=25
		beans = np.arange(0,3e3,db)
		beanc = (beans[1:]+beans[0:-1])/2
		[cts,beans] = np.histogram(s1top,beans); cts[0]=0
		pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='c+',label='sipm top')
		[cts,beans] = np.histogram(s1bot,beans); cts[0]=0
		pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='b.',label='pmt bot')
		pl.title(data_folders[ii])
		pl.xlabel('phd')
		pl.legend()
		pl.minorticks_on()
		
		pl.figure(11);pl.clf();
		pl.plot(a0t,'o')
		pl.plot(a1t,'o')
		pl.plot(a2t,'o')
		pl.plot(a3t,'o')
	
	
		pl.figure(12);pl.clf();
		pl.plot((s1top+s1bot),s1tba,'k.')
		pl.plot((s1top),s1tba,'c+')
		pl.xlim([-100,5000])
		pl.ylim([-1,1])

		# not a cut, rather a calibration of afterpulse as proxy for (saturated) S1
		pl.figure(13);pl.clf();
		pl.plot(s1bot,s1bot_ap,'o',color='powderblue')
		
		if 0:
			pl.figure(14);pl.clf();
			db=0.1
			beans = np.arange(0,10,db)
			beanc = (beans[1:]+beans[0:-1])/2
			[cts,beans] = np.histogram(s1bot/s1bot_ap,beans); cts[0]=0
			pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='k.')
	
	###
	cut = (s10>1000) & (s10<5000)

# 	cut = (s10>1000) & (s10<5000) \
# 		& (a0t<50) & (a1t<20) & (a2t<20) & (a3t<20) \
# 		& (s1tba>-0.25) & (s1tba<0.25) \
		

	nz=np.nonzero(cut); nz = nz[0]
	N = np.sum(cut)
	print('cut keeps %d events'%N)
	###

			



	### aggregate the data
	

	ct = np.array([0.01, 0.075, 0.5+0.05, 1.0+0.05, 5.0+0.05]) # trigger cascade (set by hardware)
	# Rce = 0.0063 # small-s2 limit of ratio S2ce/S2 -- should triple check
	### define the number of detected photons. subtract the number of phd identified as single e-
	dpt = np.array([ np.sum(s1top[cut])/N , np.sum(a0t[cut])*2/N , np.sum(a1t[cut])/N , np.sum(a2t[cut])/N , np.sum(a3t[cut])/N ])
	dpb = np.array([ np.sum(s1bot[cut])/N , np.sum(a0b[cut])*2/N , np.sum(a1b[cut])/N , np.sum(a2b[cut])/N , np.sum(a3b[cut])/N ])
# 	af[:,ii] = detp#/triggerS1
	
	if 1:
		pbw = 0.1 # plot bin width, fixed to cascade event window!
		pl.figure(8);#pl.clf()
		pl.errorbar(ct,dpt,yerr=np.sqrt(dpt*N)/N,xerr=np.array([0.0025,0.025,0.05,0.05,0.05]),fmt='o',color='grey',markersize=5,markerfacecolor='white',label=(labl[ii]))
		pl.errorbar(ct,dpb,yerr=np.sqrt(dpb*N)/N,xerr=np.array([0.0025,0.025,0.05,0.05,0.05]),fmt='o',color='powderblue',markersize=5,markerfacecolor='white',label=(labl[ii]))

		# photons	
		#pl.plot(tt,1.2e-4*tt**-1.3,'k-',linewidth=0.5) # PTFE
		#pl.plot(tt,0.8e-4*tt**-1.3,'k:',linewidth=0.5) # Aluminum
		pl.plot(tt,1.0*tt**-1.0,'k--',linewidth=0.5) 
		
# 		pl.plot(tt,fitdp[ii],'-',linewidth=0.5,color=colorz[ii],label=pwrlabl[ii])
				
# 		pl.plot(np.array([1e-3,1e-1]),np.array([1,1]),'k:',label='progenitor window')	
		pl.xlabel('time (ms)')
		pl.ylabel('photons / 0.1 ms')
# 		pl.ylim([1e-5,2])
		pl.xscale('log')
		pl.yscale('log')
		pl.legend()
		pl.title(data_dir[-16:-1])
		pl.ylim([0.1,5e3])


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


