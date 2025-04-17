import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob
from matplotlib.patches import Rectangle



if 1: # PTFE, Xe liquid

	data_folders = np.array(['20250416-175203']) # liquid Xe, windowless sipms S13370
	data_folders = np.array(['20250416-183422']) # liquid Xe, windowless sipms S13370

	colorz = np.array(['gray','steelblue','olivedrab','goldenrod','firebrick','sienna'])
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
print('NOTE: code assumes windowless sipms are channels 0,5,10,15 (physical channels +1)')
ww = 10 # which channel to look at

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
		print(aa_file)
		with np.load(aa_file) as data:
			a = data["arr_0"].shape[1]
			aa[:,aa_last:(aa_last+a)] = data["arr_0"]
			s1[:,aa_last:(aa_last+a)] = data["arr_1"]
			s1ap[:,aa_last:(aa_last+a)] = data["arr_2"]
			ss[:,aa_last:(aa_last+a)] = data["arr_3"]
		aa_last = aa_last + a

	ei = (h_n_events-1) # end index
	ei = int(np.floor(ei/4)*4)

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
		
		lr = 30  # lower fit range for finding gain
		ur = 100 # upper fit range for finding gain
		try:
			gains[ch] = np.average(beanc[lr:ur],axis=0,weights=cts0123[lr:ur]) # 7 picked to stay above noise bkg
		except:
			print('not enough counts for average')
					
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
			
		s1[ch,:] = s1[ch,:] / gains[ch]
# 		s1ap[ch,:] = s1ap[ch,:] / gains[ch]
# 		aa[ch,:] = aa[ch,:] / gains[ch]

#		input('paused...')

	s10 = s1[:,np.arange(0,ei,4)]
	
	ss0 = np.sum( ss[:,np.arange(0,ei,4),:] ,axis=2 )
	ss1 = np.sum( ss[:,np.arange(1,ei,4),:] ,axis=2 )
	ss2 = np.sum( ss[:,np.arange(2,ei,4),:] ,axis=2 )
	ss3 = np.sum( ss[:,np.arange(3,ei,4),:] ,axis=2 )

# 	asum = np.sum(aa,axis=0)
# 	a0t = np.sum( aa[0:nch-1,np.arange(0,ei,4)] ,axis=0 )
# 	a1t = np.sum( aa[0:nch-1,np.arange(1,ei,4)] ,axis=0 )
# 	a2t = np.sum( aa[0:nch-1,np.arange(2,ei,4)] ,axis=0 )
# 	a3t = np.sum( aa[0:nch-1,np.arange(3,ei,4)] ,axis=0 )
# 
# 	a0b = ( aa[nch-1,np.arange(0,ei,4)])
# 	a1b = ( aa[nch-1,np.arange(1,ei,4)])
# 	a2b = ( aa[nch-1,np.arange(2,ei,4)])
# 	a3b = ( aa[nch-1,np.arange(3,ei,4)])

	
# 	s1bot = (s1[pmt_ch,np.arange(0,ei,4)])
# 	if    (data_folders[0][-15:-7]=='20250227') \
# 		| (data_folders[0][-15:-7]=='20250305') \
# 		| (data_folders[0][-15:-7]=='20250306') \
# 		| (data_folders[0][-15:-7]=='20250307'): # then looking at no-xenon data
# 	if (int(data_folders[0][-15:-7])>20250227) :
# 		print('identified as no Xe data')
# 		xe=0
# 	else:
# 		s1bot = s1bot*(5/3)
# 		print('\n*** Assuming alphas in Xe (ADC saturated): accounting for PMT ADC saturation factor 5/3, obtained from afterpulse size')
# 		xe=1

# 	s1bot_ap = (s1ap[pmt_ch,np.arange(0,ei,4)])
# 	s1top = np.sum(s1[0:pmt_ch,np.arange(0,ei,4)],axis=0)
# 	s1top_ap = np.sum(s1ap[0:pmt_ch,np.arange(0,ei,4)],axis=0)
# 	s1top = s1top + s1top_ap
# 	s10 = s1top + s1bot
# 	s1tba = (s1top-s1bot)/(s1top+s1bot)
	
	if 0: ## quality cuts
		pl.figure(10);pl.clf();
		db=50
		db=50
		beans = np.arange(0,17e3,db)
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


	###
	cut = (s10[ww,:]>1500) & (s10[ww,:]<3000)
	
	if (int(data_folders[0][-6:])==161601): # the line trigger dataset
		cut = (s10>0) & (s10<100)

	nz=np.nonzero(cut); nz = nz[0]
	N = np.sum(cut)
	print('\n*** cut keeps %d events ***\n'%N)

	### aggregate the data
	ct = np.array([0.01, 0.075, 0.2+0.05, 0.5+0.05, 1.0+0.05]) # trigger cascade (set by hardware) -- afternoon Feb 20 +

	### define the number of detected photons. subtract the number of phd identified as single e-
	dpb = np.array([ np.sum(s10[ww,cut])/N , np.sum(ss0[ww,cut]/gains[ww])*2/N , np.sum(ss1[ww,cut]/gains[ww])/N , np.sum(ss2[ww,cut]/gains[ww])/N , np.sum(ss3[ww,cut]/gains[ww])/N ])
	
	if 1:
		pbw = 0.1 # plot bin width, fixed to cascade event window!
		pl.figure(8);pl.clf();ax=pl.gca()
		pl.plot(ct[0],dpb[0],'o',color='blue',markersize=7,markerfacecolor='white')
		pl.errorbar(ct[1:],dpb[1:],yerr=np.sqrt(dpb[1:]*N)/N,xerr=np.array([0.025,0.05,0.05,0.05]),fmt='o',color='blue',markersize=5,markerfacecolor='white',label='S13370')
		

		pl.text(0.012,dpb[0],('Xe scintillation trigger'),color='blue',verticalalignment='center')

		if (int(data_folders[0][-15:-7]) == 20250226): # new PMT		
# 			pl.plot(tt,1.6*tt**-1.0,'-',color='grey',linewidth=1) 
# 			pl.plot(tt,0.7*tt**-1.3,'-',color='blue',linewidth=1) 
			
			(a, b, sigma_a, sigma_b) = linfit(np.log10(ct[2:5]),np.log10(dpb[2:5]))
#			(a, b, sigma_a, sigma_b) = linfit(np.log10(ct[1:5]),np.log10(dpb[1:5]))
			print(b)
			print(1/(10**a/dpb[0]))
			pl.plot(tt[18:],(10**a)*tt[18:]**b,'-',color='blue',linewidth=1)#,label='fit to PMT') 
			pl.plot(tt[0:18],(10**a)*tt[0:18]**b,'--',color='blue',linewidth=1,label=None) 
# 			pl.plot(tt,(dpb[0]/8000)*tt**-1.5,'-',color='blue',linewidth=1) 
			pl.text(0.012,dpb[0]/200,('delayed photons $d_p(t)$'),rotation=-22,color='blue',verticalalignment='center')
# 			pl.plot(tt,(dpb[0]/4000)*tt**-1.5,'-',color='blue',linewidth=1) 

		pl.plot(tt,(dpb[0]/5313)*tt**-1.097,'--',color='blue',linewidth=1,label='R8778') # fitting last 3 points, AP<300

		pl.text(0.012,0.3*1.5,('random single photon background'),color='blue',verticalalignment='center')
		
		sig = 0.04; mu_t = 0.45*1.5; mu_b = 0.19*1.5 # factor x1.5 is data-driven
		pl.plot(np.array([1e-2,2]),np.ones(2)*mu_b,':',color='blue',linewidth=1)
		ax.add_patch(Rectangle((1e-2, mu_b-sig), 2, sig*2,color='powderblue',alpha=0.5))
		pl.plot(np.array([1e-2,2]),np.ones(2)*mu_t,':',color='grey',linewidth=1) 
		ax.add_patch(Rectangle((1e-2, mu_t-sig), 2, sig*2,color='grey',alpha=0.5))
		
				
# 		pl.plot(np.array([1e-3,1e-1]),np.array([1,1]),'k:',label='progenitor window')	
		pl.xlabel('time (ms)')
		pl.ylabel('photon counts')
# 		pl.ylim([1e-5,2])
		pl.xscale('log')
		pl.yscale('log')
		pl.legend()
# 		pl.title(data_dir[-16:-1])
		pl.ylim([0.1,1e4])

#pl.savefig('fig1.png',dpi=300)

