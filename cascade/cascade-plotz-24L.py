import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob

# 24L			S1 from BG but similar trigger to 24G

#data_folders = np.array(['20240917-154415','20240917-193113','20240917-203328','20240918-065259/','20240918-090427','20240918-100716'])
# data_dir = '/Users/peter/Public/data/20240917-154415/' # Vc = 0
# data_dir = '/Users/peter/Public/data/20240917-193113/' # Vc = 1 kV
# data_dir = '/Users/peter/Public/data/20240917-203328/' # Vc = 2 kV
# data_dir = '/Users/peter/Public/data/20240918-065259/' # Vc = 3 kV
# data_dir = '/Users/peter/Public/data/20240918-090427/' # Vc = 4 kV
# data_dir = '/Users/peter/Public/data/20240918-100716/' # Vc = 5 kV

if 0: # PTFE TPC
	data_folders = np.array(['20240919-072902','20240919-082943','20240919-093025','20240919-103105','20240919-113144','20240919-123224'])
	colorz = np.array(['gray','steelblue','olivedrab','goldenrod','firebrick','sienna'])
	labl = np.array(['0 kV','1 kV','2 kV','3 kV','4 kV','5 kV'])
	#pwrlabl = np.array([(r'$t^{-1}$'),(r'$t^{-1}$')])

if 1: # aluminum TPC
	data_folders = np.array(['20241002-170523','20241003-152154','20241009-132504','20241007-130705'])
	colorz = np.array(['gray','steelblue','olivedrab','goldenrod','firebrick','sienna'])
	labl = np.array(['Al TPC','Al, overfill','Al, remove electrodes','cold gas'])

dtt = 0.01
tt = np.arange(dtt,1.5,dtt)
fitdp = np.array([1.0e-4*tt**-1,1.2e-4*tt**-1])



af = np.zeros((5,data_folders.shape[0]))

for ii in range(3,4):#data_folders.shape[0]):
	data_dir = '/Users/peter/Public/data/'+data_folders[ii]+'/'
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

	a0 = np.sum( aa[0:32,np.arange(0,ei,4)] ,axis=0 )
	a1 = np.sum( aa[0:32,np.arange(1,ei,4)] ,axis=0 )
	a2 = np.sum( aa[0:32,np.arange(2,ei,4)] ,axis=0 )
	a3 = np.sum( aa[0:32,np.arange(3,ei,4)] ,axis=0 )

	s10 = np.sum( s1[0:32,np.arange(0,ei,4)] ,axis=0 )
	s2cf = s2cf[np.arange(0,ei,4)]

	s1bot = np.sum(s1[16:32,np.arange(0,ei,4)],axis=0)
	s1top = np.sum(s1[0:16,np.arange(0,ei,4)],axis=0)
	pl.figure(13);pl.clf();pl.plot((s1top+s1bot),(s1top-s1bot)/(s1top+s1bot),'k.')
	pl.xlim([-100,15000])
	pl.ylim([-1,1])

	coin = np.sum(ee[:,:,:]>(1/3),axis=0)
	ees = np.sum(ee[:,:,:],axis=0)
	ees = ees * (coin>13)
	ee1 = ees[np.arange(1,ei,4),:]
	ee2 = ees[np.arange(2,ei,4),:]
	ee3 = ees[np.arange(3,ei,4),:]
	
	###
	cut = (s10>5000) & (s10<15000)
# 	cut = (s10>500) & (s10<1500)

	nz=np.nonzero(cut); nz = nz[0]
	N = np.sum(cut)
	print('cut keeps %d events'%N)
	###

	pl.figure(11);pl.clf();
	db=100
	beans = np.arange(0,2.5e4,db)
	beanc = (beans[1:]+beans[0:-1])/2
	[cts,beans] = np.histogram(s10,beans); cts[0]=0
	pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='k.')
	se = 25

	### aggregate the data
	triggerS1 = np.sum(s10[cut])/N

	dpt = np.array([0.01, 0.075, 0.5+0.05, 1.0+0.05, 5.0+0.05]) # trigger cascade (set by hardware)
	# Rce = 0.0063 # small-s2 limit of ratio S2ce/S2 -- should triple check
	### define the number of detected photons. subtract the number of phd identified as single e-
	detp = np.array([ triggerS1 , np.sum(a0[cut])/N*2 , np.sum(a1[cut])/N , np.sum(a2[cut])/N , np.sum(a3[cut])/N ])
	af[:,ii] = detp#/triggerS1
	
	if 1:
		pbw = 0.1 # plot bin width, fixed to cascade event window!
		pl.figure(8);#pl.clf()
		pl.errorbar(dpt,detp/triggerS1,yerr=np.sqrt(detp*N)/N/triggerS1,xerr=np.array([0.0025,0.025,0.05,0.05,0.05]),fmt='o',color=colorz[ii],markersize=5,markerfacecolor='white',label=(labl[ii]))
# 		pl.errorbar(dpt,detp,yerr=np.sqrt(detp*N)/N,xerr=np.array([0.0025,0.025,0.1,0.1,0.1]),fmt='o',color=colorz[ii],markersize=5,markerfacecolor='white',label=(labl[ii]))

		# photons	
		pl.plot(tt,1.2e-4*tt**-1.3,'k-',linewidth=0.5) # PTFE
		pl.plot(tt,0.8e-4*tt**-1.3,'k:',linewidth=0.5) # Aluminum
		pl.plot(tt,73e-4*tt**-1.0,'k--',linewidth=0.5) # cold gas
		
# 		pl.plot(tt,fitdp[ii],'-',linewidth=0.5,color=colorz[ii],label=pwrlabl[ii])
				
# 		pl.plot(np.array([1e-3,1e-1]),np.array([1,1]),'k:',label='progenitor window')	
		pl.xlabel('time (ms)')
		pl.ylabel('photons / 0.1 ms')
# 		pl.ylim([1e-5,2])
		pl.xscale('log')
		pl.yscale('log')
		pl.legend()
		pl.title(data_dir[-16:-1])
		pl.ylim([1e-5,2])

# 		pl.figure(88);pl.clf()
# 		pl.plot(s10[cut],a0[cut],'b.')
# 		pl.xlim([1e3,16e3])
		
# 		pl.plot(np.arange(1e-2,1,1e-2),1e-4*np.arange(1e-2,1,1e-2)**-2)

	print('cascade times:')
	print(dpt)
	print('detected photons:')
	print(detp)
	print('yerr:')
	print(np.sqrt(detp*N)/N/triggerS1)
	
	if 0:
		Ef = np.array([0,1,2,3,4,5]) # electric field
		pl.figure(7);pl.clf()
		pl.plot(Ef,af[4,:],'ro')


# pl.figure(8);
# pl.plot(np.array([0.01,10]),1e-4*np.array([1,1]),'k:',label='accidental')	

