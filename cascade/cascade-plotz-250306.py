import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob
from matplotlib.patches import Rectangle
from matplotlib import rc
rc('font', **{'family':'sans-serif','serif':['Palatino']})
rc('text', usetex=True)
rc('xtick', labelsize=14)
rc('ytick', labelsize=14)



pmt_ch = 16
ds_pmt = 5 # pmt signal was scaled down for the analysis, put it back
if 1: # PTFE, Xe liquid, alphas from cathode pointing down I think, PMT bottom, compact cylinder, not really a TPC

	data_folders = np.array(['20250226-202029']) # liquid Xe, new PMT
# 	data_folders = np.array(['20250227-070842']) # liquid Xe, new PMT

#	data_folders = np.array(['20250311-161601']) # 165 K no Xe, sipms at 50 V, line trigger

	colorz = np.array(['blue','dodgerblue','olivedrab','goldenrod','firebrick','sienna'])
	labl = np.array(['','','',''])


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

dtt = 0.01
tt = np.arange(0.012,1.5,dtt)
fitdp = np.array([1.0e-4*tt**-1,1.2e-4*tt**-1])



af = np.zeros((5,data_folders.shape[0]))
cutvals = np.arange(1000,5000,500) # for PMT
#cutvals = np.arange(1500,11000,1000) # for sipm
cutvals = np.array([1500,2500,3500,4500,5500,6500])
mvals = np.zeros(len(cutvals))
avals = np.zeros(len(cutvals))
bvals = np.zeros(len(cutvals))
sumdp = np.zeros(len(cutvals))
mvals1 = np.zeros(len(cutvals))
avals1 = np.zeros(len(cutvals))
bvals1 = np.zeros(len(cutvals))
sumdp1 = np.zeros(len(cutvals))

# for cc in range(0,len(cutvals)):
for cc in range(0,1):
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
# 			print(aa_file)
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
			[cts2,beans] = np.histogram(ss[ch,np.arange(2,ei,4),:],beans); cts2[0]=0 
			[cts3,beans] = np.histogram(ss[ch,np.arange(3,ei,4),:],beans); cts3[0]=0 
			cts0123 = cts0+cts1+cts2+cts3
			cts123 = cts1+cts2+cts3
		
			lr = 7  # lower fit range for finding gain
			ur = 55 # upper fit range for finding gain
			if (ch<pmt_ch):
				gains[ch] = np.average(beanc[lr:ur],axis=0,weights=cts0123[lr:ur]) # 7 picked to stay above noise bkg
			else:
				gains[ch] = np.average(beanc[lr:],axis=0,weights=cts1[lr:]) # 7 picked to stay above noise bkg
		
			if 0:
				pl.figure(6);pl.clf();
				pl.errorbar(beanc,cts0,yerr=np.sqrt(cts0),xerr=db/2,fmt='k.')
				pl.errorbar(beanc,cts1,yerr=np.sqrt(cts1),xerr=db/2,fmt='r.')
				pl.errorbar(beanc,cts2,yerr=np.sqrt(cts2),xerr=db/2,fmt='b.')
				pl.errorbar(beanc,cts3,yerr=np.sqrt(cts3),xerr=db/2,fmt='g.')

				pl.errorbar(beanc,cts0123,yerr=np.sqrt(cts0123),xerr=db/2,fmt='.',color='gray')
	
				pl.plot(np.ones(2)*gains[ch],np.array([1,max(cts0)]),':',color='gray')
				pl.plot(np.ones(2)*beanc[lr],np.array([1,max(cts0)]),'-',color='gray')
				pl.plot(np.ones(2)*beanc[ur],np.array([1,max(cts0)]),'-',color='gray')
				

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
		onech=0
		if onech:
			a0t = ( aa[8,np.arange(0,ei,4)] )
			a1t = ( aa[8,np.arange(1,ei,4)] )
			a2t = ( aa[8,np.arange(2,ei,4)] )
			a3t = ( aa[8,np.arange(3,ei,4)] )
		else:
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
		if (int(data_folders[0][-15:-7])>20250227) :
			print('identified as no Xe data')
			xe=0
		else:
			s1bot = s1bot*(5/3)
			print('\n*** Assuming alphas in Xe (ADC saturated): accounting for PMT ADC saturation factor 5/3, obtained from afterpulse size')
			xe=1

		s1bot_ap = (s1ap[pmt_ch,np.arange(0,ei,4)])
		
		if onech:
			s1top = (s1[8,np.arange(0,ei,4)])
			s1top_ap = (s1ap[8,np.arange(0,ei,4)])		
		else:
			s1top = np.sum(s1[0:pmt_ch,np.arange(0,ei,4)],axis=0)
			s1top_ap = np.sum(s1ap[0:pmt_ch,np.arange(0,ei,4)],axis=0)

		s1top = s1top + s1top_ap
		s10 = s1top + s1bot
		s1tba = (s1top-s1bot)/(s1top+s1bot)
	

		###
		if xe: # for UCLA cut on 4500 top and bot
			cut = (s1bot>0) & (s1bot<4500) \
				& (a0b<100) & (a1b<100) & (a2b<100) & (a3b<100) \
				& (s1top>0) & (s1top<5500)
# 				& (s1top>cutvals[cc-1]) & (s1top<cutvals[cc]) \
# 				& (s1top>000) & (s1top<4500) \
	# 	 		& (s1bot_ap<300) # used 300 for UCLA, but don't need; previous line takes care of that data issue
	# 	 		& (a1t<(a0t*2)) # used for UCLA but maybe don't need

	# 			& (a1t<a0t) & (a2t<a0t) & (a3t<a0t) & (a1b<a0b) & (a2b<a0b) & (a3b<a0b)
	
		if (int(data_folders[0][-6:])==161601): # the line trigger dataset
			cut = (s10>0) & (s10<100)

		nz=np.nonzero(cut); nz = nz[0]
		N = np.sum(cut)
		print('\n*** cut keeps %d events ***\n'%N)

		### aggregate the data
	

	# 	ct = np.array([0.01, 0.075, 0.5+0.05, 1.0+0.05, 5.0+0.05]) # trigger cascade (set by hardware) -- used for ever all earlier data
		ct = np.array([0.01, 0.075, 0.2+0.05, 0.5+0.05, 1.0+0.05]) # trigger cascade (set by hardware) -- afternoon Feb 20 +
		if (data_folders[0][-6:]=='094512'):
			ct = np.array([0.01, 0.075, 0.3+0.05, 0.5+0.05, 1.0+0.05]) # trigger cascade (set by hardware) -- briefly used before noon on Feb 20

		# Rce = 0.0063 # small-s2 limit of ratio S2ce/S2 -- should triple check
		### define the number of detected photons. subtract the number of phd identified as single e-
		dpt = np.array([ np.sum(s1top[cut])/N , np.sum(a0t[cut])*2/N , np.sum(a1t[cut])/N , np.sum(a2t[cut])/N , np.sum(a3t[cut])/N ])
		dpb = np.array([ np.sum(s1bot[cut])/N , np.sum(a0b[cut])*2/N , np.sum(a1b[cut])/N , np.sum(a2b[cut])/N , np.sum(a3b[cut])/N ])
	# 	af[:,ii] = detp#/triggerS1
	
 
		if 1:
			pbw = 0.1 # plot bin width, fixed to cascade event window!
			pl.figure(8);pl.clf();ax=pl.gca()
			pl.errorbar(ct[1:],dpt[1:],yerr=np.sqrt(dpt[1:]*N)/N,xerr=np.array([0.025,0.05,0.05,0.05]),fmt='s',color='grey',markersize=5,markerfacecolor='white',lw=0.5,label='S13371 SiPM')
			pl.errorbar(ct[1:],dpb[1:],yerr=np.sqrt(dpb[1:]*N)/N,xerr=np.array([0.025,0.05,0.05,0.05]),fmt='o',color=colorz[ii],markersize=5,markerfacecolor='white',lw=0.5,label='R8778 PMT')
		
			if (int(data_folders[0][-6:])!=161601):
				pl.plot(ct[0],dpt[0],'s',color='grey',markersize=7,markerfacecolor='white')
				pl.plot(ct[0],dpb[0],'o',color=colorz[ii],markersize=7,markerfacecolor='white')

	# 		pl.text(0.015,dpb[0],('Xe scintillation trigger $\mu=%1.0f$'%dpb[0]),color='blue',verticalalignment='center')
			pl.text(0.012,dpb[0],(r'Xe scintillation pulse $\bar{a}$'),color='k',verticalalignment='center',size=14)

			if (int(data_folders[0][-15:-7]) == 20250226) | (int(data_folders[0][-15:-7]) == 20250227): # new PMT		
	# 			pl.plot(tt,1.6*tt**-1.0,'-',color='grey',linewidth=1) 
	# 			pl.plot(tt,0.7*tt**-1.3,'-',color='blue',linewidth=1) 
			
				(a, b, sigma_a, sigma_b) = linfit(np.log10(ct[2:5]),np.log10(dpb[2:5]))
	# 			(a, b, sigma_a, sigma_b) = linfit(np.log10(ct[1:5]),np.log10(dpb[1:5]))
				print(b)
				print(1/(10**a/dpb[0]))
				pl.plot(tt[0:],(10**a)*tt[0:]**b,'--',color=colorz[ii],linewidth=0.5)#,label='fit to PMT') 
# 				pl.plot(tt[0:18],(10**a)*tt[0:18]**b,'--',color=colorz[ii],linewidth=0.5,label=None) 
				pl.text(0.012,dpb[0]/150,('delayed photons $d_p(t)$'),rotation=-24,color='k',verticalalignment='center',size=14)

			pl.text(0.012,0.3*1.5,('random photon background'),color='k',verticalalignment='center',size=14)
		
			sig = 0.04; mu_t = 0.45*1.5; mu_b = 0.19*1.5 # factor x1.5 is data-driven
			pl.plot(np.array([1e-2,2]),np.ones(2)*mu_b,':',color='grey',linewidth=1)
			ax.add_patch(Rectangle((1e-2, mu_b-sig), 2, sig*2,facecolor=colorz[ii],alpha=0.25,edgecolor='None'))
			pl.plot(np.array([1e-2,2]),np.ones(2)*mu_t,':',color='grey',linewidth=1) 
			ax.add_patch(Rectangle((1e-2, mu_t-sig), 2, sig*2,facecolor='grey',alpha=0.25,edgecolor='None'))

			if 1: # top sipm array
				(a1, b1, sigma_a, sigma_b) = linfit(np.log10(ct[1:5]),np.log10(dpt[1:5]))
				print(b1)
				print(1/(10**a1/dpb[0]))
				pl.plot(tt[0:],(10**a1)*tt[0:]**b1,'--',color='grey',linewidth=0.5)#,label='fit to sipms') 
		
				
	# 		pl.plot(np.array([1e-3,1e-1]),np.array([1,1]),'k:',label='progenitor window')	
			pl.xlabel('time (ms)',fontsize=14)
			pl.ylabel('photon counts',fontsize=14)
	# 		pl.ylim([1e-5,2])
			pl.xscale('log')
			pl.yscale('log')
			pl.legend(fontsize=14)
	# 		pl.title(data_dir[-16:-1])
			pl.ylim([0.1,1e4])

	if cutvals[cc]==5500:
		print('saving fig1')
		pl.savefig('fig1.png',dpi=300)

	if 1:
		mvals[cc] = dpb[0] # mean value of bottom pmt hit as function of cut value on top sipm array
		avals[cc] = 1/(10**a/dpb[0])
		bvals[cc] = b
		sumdp[cc] = np.sum(mvals[cc]/avals[cc]*tt[0:99]**bvals[cc])

		mvals1[cc] = dpt[0] # mean value of sipm hit as function of cut value on top sipm array
		avals1[cc] = 1/(10**a1/dpt[0])
		bvals1[cc] = b1
		sumdp1[cc] = np.sum(mvals1[cc]/avals1[cc]*tt[0:99]**bvals1[cc])

# 		input('press a key')

# with bot pmt seeing<4500, cut on top sipm array seeing < 1500,2000,2500,etc
# remove data point from <8500
#tslt = np.array([1500,2000,2500,3000,3500,4000,4500,5000,5500,6000,6500,7000,7500,8000,9000,9500,10000,10500]) # top sipm less than tslt
#mm =  np.array([878,1030,1243,1542,1693,1920,2091,2430,2704,2735,3068,3393,3755,4194,4460,4895,5300,5888])
#aa = np.array([3153,3731,4006,4656,4830,5099,5186,5387,5518,5515,5580,5581,5422,5160,4843,4592,4302,3718])
#bb = -np.array([0.776,0.838,0.862,0.951,0.990,1.081,1.095,1.180,1.260,1.260,1.323,1.362,1.402,1.421,1.456,1.457,1.481,1.472])



pl.figure(88);pl.clf();ax=pl.gca()
pl.plot(mvals1[1:],sumdp1[1:]/np.max(sumdp1),'o',label=r'$\Sigma(at^{-b}$) in sipm')
#pl.plot(mvals[1:],avals[1:],'o',label='1/a in PMT')
pl.plot(mvals1[1:],-bvals1[1:],'o',label='-b in sipm')

pl.plot(mvals1[1:],sumdp[1:]/np.max(sumdp),'o',label=r'$\Sigma(at^{-b}$) in pmt')
#pl.plot(mvals[1:],avals[1:],'o',label='1/a in PMT')
pl.plot(mvals1[1:],-bvals[1:],'o',label='-b in pmt')

pl.xlabel('S13371 average photons per pulse ')
pl.ylabel('see legend')
pl.legend()
#(a, b, sigma_a, sigma_b) = linfit(mm[6:],sumdp[6:])
#pl.plot(np.arange(0,6000,10),a+b*np.arange(0,6000,10),'k--')
#pl.xlim([0,6600])
#pl.ylim([0,1100])

#pl.figure(89);pl.clf();ax=pl.gca()
#pl.plot(mvals,bvals,'o')



# quality cuts below
if 0: ## quality cuts
	pl.figure(10);pl.clf();
	db=50
	if xe==1:
		db=50
		beans = np.arange(0,17e3,db)
	else:
		db=5
		beans = np.arange(0,500,db)		
	beanc = (beans[1:]+beans[0:-1])/2
	[cts,beans] = np.histogram(s1top,beans); cts[0]=0
	pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='c+',label='sipm top')
	[cts,beans] = np.histogram(s1bot,beans); cts[0]=0
	pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='b.',label='pmt bot')
	[cts,beans] = np.histogram(s10,beans); cts[0]=0
	pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='k.',label='S10')
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


	if 0:
		pl.figure(12);pl.clf();
		pl.plot((s1top+s1bot),s1tba,'k.')
# 		pl.plot((s1top),s1tba,'c+')
		pl.xlim([-100,12e3])
		pl.ylim([-1,1])

	# not a cut, rather a calibration of afterpulse as proxy for (saturated) S1
	if 0:
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

	if 1:
		pl.figure(15);pl.clf();
		pl.plot(s1bot,s1top,'o',color='gray',markersize=1)

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






