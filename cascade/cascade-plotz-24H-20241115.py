import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob

# 20241115 - pfs - rewritten to tell the LED story
#
#
#

pl.figure(9, figsize=(10,4));pl.clf()
for ds in range(0,3):
	# 24B
	#data_dir = '/Users/peter/Public/data/20240314-190415/' # first good cascade data, 0.5,1.0,5.0 ms

	# 24C
	#data_dir = '/Users/peter/Public/data/20240319-150343/' # after open to air, 0.5,1.0,5.0 ms

	# 24D
	#data_dir = '/Users/peter/Public/data/20240411-103253/' # after 100C bake overnight 0.5,1.0,5.0 ms

	# 24G
	#data_dir = '/Users/peter/Public/data/20240604-163805/' # regular cascade

	# 24H
	if ds==0:
		data_dir = '/Users/peter/Public/data/20240611-173543/' # 133Ba cascade
	if ds==1:
		data_dir = '/Users/peter/Public/data/20240612-173658/' # 133Ba+LED0.5us cascade	
	if ds==2:
		data_dir = '/Users/peter/Public/data/20240613-155550/' # LED0.5us cascade	

	# 24L			S1 from BG but similar trigger to 24G
	#data_dir = '/Users/peter/Public/data/20240917-154415/' # Vc = 0
	#data_dir = '/Users/peter/Public/data/20240917-193113/' # Vc = 1 kV
	#data_dir = '/Users/peter/Public/data/20240917-203328/' # Vc = 2 kV

	aa_file_list = glob.glob(data_dir+"./aa/*v1.npz")
	print('found %d files'%len(aa_file_list))
	h_file = np.load(data_dir+"/compressed_filtered_data/headers.npz")
	h_array = h_file["arr_0"]
	h_n_events = int(np.floor(h_array.size/8))

	aa = np.zeros([32,h_n_events])
	ee = np.zeros([32,h_n_events,20])
	s1 = np.zeros([32,h_n_events])
	is1 = np.zeros([h_n_events])
	is2 = np.zeros([h_n_events])
	s2cf = np.zeros([h_n_events])
	# ss = np.zeros([25000,h_n_events])
	aa_last = 0
	for aa_file in aa_file_list:
		#print('load file %s'%aa_file)
	#	np.savez(aa_file,aa,ee,s1,is1,is2,s2f)
		with np.load(aa_file) as data:
			a = data["arr_0"].shape[1]
			aa[:,aa_last:(aa_last+a)] = data["arr_0"]
			ee[:,aa_last:(aa_last+a),:] = data["arr_1"]
			s1[:,aa_last:(aa_last+a)] = data["arr_2"]
			is1[aa_last:(aa_last+a)] = data["arr_3"]
			is2[aa_last:(aa_last+a)] = data["arr_4"]
			s2cf[aa_last:(aa_last+a)] = data["arr_5"]
		aa_last = aa_last + a

	asum = np.sum(aa,axis=0)
	ei = (asum.shape[0]-1) # end index
	ei = int(np.floor(ei/4)*4)

	# parse the photon and electron counts from the 4-step cascade
	a0 = ( asum[np.arange(0,ei,4)] )
	a1 = ( asum[np.arange(1,ei,4)] )
	a2 = ( asum[np.arange(2,ei,4)] )
	a3 = ( asum[np.arange(3,ei,4)] )
	s10 = np.sum( s1[0:32,np.arange(0,ei,4)] ,axis=0 )

	# parse the electron counts
	coin = np.sum(ee[:,:,:]>(1/3),axis=0)
	ees = np.sum(ee[:,:,:],axis=0)
	ees = ees * (coin>13)
	ee1 = ees[np.arange(1,ei,4),:]
	ee2 = ees[np.arange(2,ei,4),:]
	ee3 = ees[np.arange(3,ei,4),:]

	# a few more variables that aren't being used
	# s2cf = s2cf[np.arange(0,ei,4)]
	# s1bot = np.sum(s1[16:32,np.arange(0,ei,4)],axis=0)
	# s1top = np.sum(s1[0:16,np.arange(0,ei,4)],axis=0)
	

	dtt = 0.1
	tt = np.arange(0.001,1000,dtt)

	###
# 	cut = (s10>1e5) & (s10<2.5e5) & (a1>0)# seem to need this for 20240611-173543
	cut = (s10>10) & (s10<2.5e5) & (a1>0)# seem to need this for 20240611-173543

	nz=np.nonzero(cut); nz = nz[0]
	N = np.sum(cut)
	print('cut keeps %d events'%N)


	###

	pl.figure(10);pl.clf();
	db=3
	beans = np.arange(0,100,db)
	beanc = (beans[1:]+beans[0:-1])/2
	[cts,beans] = np.histogram(np.concatenate([ee1[cut,:],ee2[cut,:],ee3[cut,:]]),beans); cts[0]=0
	pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='k.')
	se = 25
	#se = np.average(beanc,axis=0,weights=cts) # over-estimate

	pl.figure(11);pl.clf();
	db=5e3
	beans = np.arange(0,5e5,db)
	beanc = (beans[1:]+beans[0:-1])/2
	[cts,beans] = np.histogram(s10,beans); cts[0]=0
	pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='k.')

	# get the top-bot ratio of the S2 from single e-
	eetop = np.sum(ee[0:16,:,:],axis=0) * (coin>13)
	eebot = np.sum(ee[16:32,:,:],axis=0) * (coin>13)
	eeta = eetop[np.arange(1,ei,4),:] + eetop[np.arange(2,ei,4),:] + eetop[np.arange(3,ei,4),:]
	eeba = eebot[np.arange(1,ei,4),:] + eebot[np.arange(2,ei,4),:] + eebot[np.arange(3,ei,4),:]
	eebfrac = eeba/(eeta+eeba)
	pl.figure(12);pl.clf();
	#pl.plot((eeta[cut,:]+eeba[cut,:]),eebfrac[cut,:],'k.')
	db = 0.03
	beans = np.arange(0,1,db)
	beanc = (beans[1:]+beans[0:-1])/2
	[cts,beans] = np.histogram(eebfrac[cut,:],beans); cts[0]=0
	pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='k.')
	bfrac = 0.44 # previous estimate was 0.33 from higher energy data

	### aggregate the data
	dpt = np.array([0.01, 0.55, 1.05, 5.05]) # trigger cascade (set by hardware)

	### define the number of detected photons. subtract the number of phd identified as single e-
	detp = np.array([ np.sum(s10[cut])/N/bfrac , np.sum(a1[cut]-np.sum(ee1[cut,:],axis=1))/N , np.sum(a2[cut]-np.sum(ee2[cut,:],axis=1))/N , np.sum(a3[cut]-np.sum(ee3[cut,:],axis=1))/N ])

	### number of detected electrons in the cascades
	dete = np.array([ detp[0] , np.sum(ee1[cut,:])/N , np.sum(ee2[cut,:])/N , np.sum(ee3[cut,:])/N ]) / se 

	ptrainf = detp[1]/detp[0]*100
	print("photon train fraction in 100 us window at 500 us delayt is %1.4f"%ptrainf)

	# pl.figure(100);pl.clf();pl.plot(s10[cut],a1[cut]/s10[cut],'ko')


	if 1:
		pbw = 0.1 # plot bin width, fixed to cascade event window!
		pl.figure(9);
		pl.subplot(1,2,1) # photons
		pl.errorbar(dpt,detp/detp[0],xerr=np.array([0.05,0.05,0.05,0.05]),yerr=np.sqrt(detp*N)/N/detp[0],fmt='o',color='steelblue',markersize=5,markerfacecolor='white') #,label=(r'photon signal')
		if ds==0:
			pl.plot(tt,1e-4/tt,'-',linewidth=0.5,color='steelblue',label='S2')#,label=(r'$t^{-1}$ ($\gamma$)')
		if ds==1:
			pl.plot(tt,10e-4*tt**-1.3,'--',linewidth=1,color='steelblue',label='S2+LED')#,label=(r'$t^{-1.3}$ ($\gamma$)')
		if ds==1:
			pl.plot(tt,100e-4*tt**-1.3,':',linewidth=0.5,color='steelblue',label='LED only')#,label=(r'$t^{-1.3}$ ($\gamma$)')

		pl.xscale('log')
		pl.yscale('log')
		pl.axis([0.5e-3,10,0.1e-5,15e-1])
		pl.legend()
		pl.xlabel('time (ms)')
		pl.ylabel('detected photons / 0.1 ms')

		pl.subplot(1,2,2) # electrons
		pl.errorbar(dpt,dete/detp[0],xerr=np.array([0.05,0.05,0.05,0.05]),yerr=np.sqrt(dete*N)/N/detp[0],fmt='o',color='darkorange',markersize=5,markerfacecolor='white')#,label=(r'electron signal')
		if ds==0:
			pl.plot(tt,0.6e-5*tt**-0.5,'-',linewidth=0.5,color='darkorange',label='S2')#,label=(r'$t^{-%1.1f}$ (e-)'%0.5)
		if ds==1:
			pl.plot(tt,2e-5*tt**-0.5,'--',linewidth=1,color='darkorange',label='S2+LED')#,label=(r'$t^{-%1.1f}$ (e-)'%0.5)
		if ds==2:
			pl.plot(tt,0.22e-5*tt**-0.2,':',linewidth=0.5,color='darkorange',label='LED only')#,label=(r'$t^{-%1.1f}$ (e-)'%0.5)
# 			pl.plot(tt,0.2e-5*np.ones(len(tt)),':',linewidth=0.5,color='darkorange',label='LED only')#,label=(r'$t^{-%1.1f}$ (e-)'%0.5)
			
		pl.xscale('log')
		pl.yscale('log')
		pl.axis([0.5e-3,10,0.1e-5,15e-1])
		pl.legend()
		pl.xlabel('time (ms)')
		pl.ylabel('detected electrons / 0.1 ms')

#		pl.title(data_dir[-16:-1])






