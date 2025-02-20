import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob


# data_dir = '/media/xaber/extradrive1/crystalize_data/data-202403/20240314/20240314-190415/'
# data_dir = '/Users/peter/Public/data/20240314-190415/'
# data_dir = '/Users/peter/Public/data/20240319-150343/' # after opening to air

data_dir = '/Users/peter/Public/data/20240411-103253/' # after 12 hr bake at 100C
dpt = np.array([0.0001, 0.5+0.05, 1.0+0.05, 5.0+0.05]) # trigger cascade (set by hardware)

data_dir = '/Users/peter/Public/data/20240411-095815/' # after 12 hr bake at 100C, 0.15, 0.30, 0.50 ms cascade
dpt = np.array([0.0001, 0.15+0.05, 0.30+0.05, 0.50+0.05]) # trigger cascade (set by hardware)
# 
Rce = 0.0063 # small-s2 limit of ratio S2ce/S2 -- should triple check

aa_file_list = glob.glob(data_dir+"./aa/*v1.npz")
print('found %d files'%len(aa_file_list))
h_file = np.load(data_dir+"/compressed_filtered_data/headers.npz")
h_array = h_file["arr_0"]
h_n_events = int(np.floor(h_array.size/8))

aa = np.zeros([32,h_n_events])
ee = np.zeros([32,h_n_events])
s1 = np.zeros([32,h_n_events])
is1 = np.zeros([h_n_events])
is2 = np.zeros([h_n_events])
s2f = np.zeros([h_n_events])
eec = np.zeros([h_n_events])
ss = np.zeros([50000,h_n_events])
aa_last = 0
for aa_file in aa_file_list:
	#print('load file %s'%aa_file)
#	np.savez(aa_file,aa,ee,s1,is1,is2,s2f)
	with np.load(aa_file) as data:
		a = data["arr_0"].shape[1]
		aa[:,aa_last:(aa_last+a)] = data["arr_0"]
		ee[:,aa_last:(aa_last+a)] = data["arr_1"]
		s1[:,aa_last:(aa_last+a)] = data["arr_2"]
#		is1[:,aa_last:(aa_last+data["arr_3"].shape[1])] = data["arr_3"]
# 		is2[:,aa_last:(aa_last+data["arr_4"].shape[1])] = data["arr_4"]
		s2f[aa_last:(aa_last+a)] = data["arr_5"]
# 		eec[aa_last:(aa_last+a)] = data["arr_6"]
		try:
			ss[:,aa_last:(aa_last+a)] = data["arr_7"]
		except:
			print("assign ss failed")		        
# 		a = data["arr_0"]
# 		e = data["arr_1"]
# 		g = data["arr_2"]
# 	aa[:,aa_last:(aa_last+a.shape[1])] = a
# 	ee[:,aa_last:(aa_last+a.shape[1])] = e
# 	s1[:,aa_last:(aa_last+a.shape[1])] = g
	aa_last = aa_last + a

asum = np.sum(aa,axis=0)
ei = (asum.shape[0]-1) # end index
ei = int(np.floor(ei/4)*4)
a0 = ( asum[np.arange(0,ei,4)] )
a1 = ( asum[np.arange(1,ei,4)] )
a2 = ( asum[np.arange(2,ei,4)] )
a3 = ( asum[np.arange(3,ei,4)] )
s2f = s2f[np.arange(0,ei,4)]

ee0 = np.sum(ee[:,np.arange(0,ei,4)],axis=0)
ee1 = np.sum(ee[:,np.arange(1,ei,4)],axis=0)
ee2 = np.sum(ee[:,np.arange(2,ei,4)],axis=0)
ee3 = np.sum(ee[:,np.arange(3,ei,4)],axis=0)

try:
	ss0 = ss[:,np.arange(0,ei,4)]
	ss1 = ss[:,np.arange(1,ei,4)]
	ss2 = ss[:,np.arange(2,ei,4)]
	ss3 = ss[:,np.arange(3,ei,4)]
except:
	print("ss variable does not exist or wrong size")
	
###
cut = (a0>1000) & (a0<2000) & (s2f>0.67)
# cut = (a0>1500) & (a0<2500) & (s2f>0.67)
N = np.sum(cut)
print('cut keeps %d events'%N)
###


if 1:
	db=7
	beans = np.arange(0,100,db)
	beanc = (beans[1:]+beans[0:-1])/2
	pl.figure(10);pl.clf();
	# [cts,beans] = np.histogram(np.concatenate([ee1,ee2,ee3]),beans); cts[0]=0
	# pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='b.')
	[cts,beans] = np.histogram(np.concatenate([ee1[cut],ee2[cut],ee3[cut]]),beans); cts[0]=0
	pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='k.')
	se = 29
	#se = np.average(beanc,axis=0,weights=cts) # over-estimate

### aggregate the data


### define the number of detected photons. subtract the number of phd identified as single e-
detp = np.array([ np.sum(a0[cut])/Rce , np.sum(a1[cut]-ee1[cut]) , np.sum(a2[cut]-ee2[cut]) , np.sum(a3[cut]-ee3[cut]) ]) 
#detp = np.array([ np.sum(a0[cut])/0.01 , np.sum(a1[cut]) , np.sum(a2[cut]) , np.sum(a3[cut]) ]) 
detp = detp/N # units are detected photons per 0.1 ms


### number of detected electrons in the cascades
dete = np.array([ detp[0]/se , np.sum(ee1>0)/N , np.sum(ee2>0)/N , np.sum(ee3>0)/N ]) # this is too simplistic, need to improve lower-level algorithm

dtt = 0.1
# dtt = 1
tt = np.arange(0.0001,1000,dtt)



if (data_dir[-7:-1]=='103253'):
	fit = 35/tt # photons
	elz = 1.85/tt # electrons, LZ based on 1% target
	elbl = 0.30/tt
if (data_dir[-7:-1]=='095815'):
	fit = 35/tt # photons
	elz = 1.85/tt # electrons, LZ based on 1% target
	elbl = 0.30/tt



if 1:
	pbw = 0.1 # plot bin width, fixed to cascade event window!
	pl.figure(1);#pl.clf()
	
	if (data_dir[-7:-1]=='095815'): # short-time cascade
		# photons
		pl.errorbar(dpt,detp,yerr=np.sqrt(detp*N)/N,fmt='o',color='steelblue',markersize=5,markerfacecolor='white',label=(r'photon signal'))
		texp = np.arange(0.01,0.99,0.001)
# 		pl.plot(texp,2.54e5*np.exp(-(texp)/0.03),'--',color='black',linewidth=0.5)
# 		pl.plot(texp,2.5e6*np.exp(-(texp)/0.02),'--',color='black',linewidth=0.5)
# 		pl.plot(texp,1.18e5*np.exp(-(texp)/0.12),'--',color='black',linewidth=0.5)
# 		pl.plot(texp,5.0e4*np.exp(-(texp-0.07)/0.05),'--',color='black',linewidth=0.5)
		pl.plot(texp,7.5e4*np.exp(-(texp-0.07)/0.02),'--',color='black',linewidth=0.5)
		# electrons
		pl.errorbar(dpt,dete,yerr=np.sqrt(dete*N)/N,fmt='o',color='darkorange',markersize=5,markerfacecolor='white',label=(r'electron signal'))

		ttss0 = np.arange(0.05,0.1-1e-6,2e-6) # based on what was defined in cascade-analysis.py
		ss0sum = np.sum(ss0[ttss0.shape[0]:,cut],axis=1)
		ss1sum = np.sum(ss1[:,cut],axis=1)
		ss2sum = np.sum(ss2[:,cut],axis=1)
		ss3sum = np.sum(ss3[:,cut],axis=1)

# 		detp0 = np.mean(ss0sum/N*25000) # single chunk of 50 us
		detp0 = np.sum(ss0sum)/N*2 # single chunk of 50 us
		pl.errorbar(0.075,detp0*0.9,yerr=np.sqrt(detp0*N)/N,fmt='o',linewidth=0.5,color='steelblue',markersize=5,markerfacecolor='white')

		if 1:
			# break the 50 us back half of the initial trigger waveform into chunks
			detp0 = np.array([ np.sum(ss0sum[0:12500]) , np.sum(ss0sum[12500:-15]) ])/N *4 # factor x4? bc this is 50 us/2 and I am plotting rates/100 us
			# plot the integrated photons signal *0.9 to account for SE hiding in the traces
			pl.errorbar(np.array([0.0625,0.0875]),detp0*0.9,yerr=np.sqrt(detp0*N)/N,fmt='o',linewidth=0.5,color='steelblue',markersize=5,markerfacecolor='white')

			# break the 100 us  of the 1st cascade waveform into chunks
			detp1 = np.array([ np.sum(ss1sum[0:25000]) , np.sum(ss1sum[25000:-15]) ])/N *2 # factor x2? bc this is 100 us/2 and I am plotting rates/100 us
			# plot the integrated photons signal *0.9 to account for SE hiding in the traces
			pl.errorbar(np.array([0.175,0.225]),detp1*0.9,yerr=np.sqrt(detp1*N)/N,fmt='o',linewidth=0.5,color='steelblue',markersize=5,markerfacecolor='white')

			# break the 100 us  of the 2nd cascade waveform into chunks
			detp2 = np.array([ np.sum(ss2sum[0:25000]) , np.sum(ss2sum[25000:-15]) ])/N *2 # factor x2? bc this is 100 us/2 and I am plotting rates/100 us
			# plot the integrated photons signal *0.9 to account for SE hiding in the traces
			pl.errorbar(np.array([0.325,0.375]),detp2*0.9,yerr=np.sqrt(detp2*N)/N,fmt='o',linewidth=0.5,color='steelblue',markersize=5,markerfacecolor='white')


		if 1:	# plot entire waveform snippets if available
			pl.plot(ttss0[:-15],ss0sum[:-15]/N*25000*2,'-',color='steelblue')

			ttss1 = np.arange(0.15,0.25-1e-6,2e-6) # based on what was defined in cascade-analysis.py
			pl.plot(ttss1[15:-15],ss1sum[15:-15]/N*50000,'-',color='steelblue')

			ttss2 = np.arange(0.30,0.40-1e-6,2e-6) # based on what was defined in cascade-analysis.py
			pl.plot(ttss2[:],ss2sum[:]/N*50000,'-',color='steelblue')

			ttss3 = np.arange(0.50,0.60-1e-6,2e-6) # based on what was defined in cascade-analysis.py
			pl.plot(ttss3[:],ss3sum[:]/N*50000,'-',color='steelblue')
	

	if (data_dir[-7:-1]=='103253'):	# longer-times cascade
		pl.errorbar(dpt,detp,yerr=np.sqrt(detp*N)/N,fmt='o',color='steelblue',markersize=5,markerfacecolor='white')
		pl.errorbar(dpt,dete,yerr=np.sqrt(dete*N)/N,fmt='o',color='darkorange',markersize=5,markerfacecolor='white')
		# photons -- plot 1/t	
		pl.step(tt,fit,'-',linewidth=0.5,color='steelblue',label=(r'$1/t$ ($\gamma$)'))
		delayed_p_3ms = np.sum(fit[tt>3])*dtt/pbw
		print('delayed photons in 3-1000 ms window: %1.2f%% (%1.0f total)'% (delayed_p_3ms/(detp[0])*100,delayed_p_3ms) )
#	 	pl.step(tt,50*tt**-0.5,'-',linewidth=0.5,color='steelblue',label=(r'$1/t$ ($\gamma$)'))

		# electrons -- plot 1/t
		pl.step(tt,elbl,'-',linewidth=0.5,color='darkorange',label=(r'$1/t$ (e-)'))
		lbl_delayed_e_3ms = np.sum(elbl[tt>3])*dtt/pbw
		print('(LBL) delayed electrons in 3-1000 ms window: %1.2f%% (%1.0f total)'% (lbl_delayed_e_3ms/dete[0]*100,lbl_delayed_e_3ms) )
# 	 	pl.plot(3*np.array([1,1]),np.array([1e-3,fit[tt>3][0]]),'k:')


	
	pl.xlabel('time (ms)')
	pl.ylabel('(e or p) / 0.1 ms')
	pl.ylim([1e-3,1e6])
	pl.xscale('log')
	pl.yscale('log')
	pl.legend()
	pl.title(data_dir[-16:-1])

if 0:
	pl.figure(2);pl.clf()
	pl.plot(ttss0[:-15],ss0sum[:-15]/N*25000,'-',color='steelblue')
	pl.plot(ttss1[:],ss1sum[:]/N*50000,'-',color='steelblue')
	pl.plot(ttss2[:],ss2sum[:]/N*50000,'-',color='steelblue')
	pl.plot(ttss3[:],ss3sum[:]/N*50000,'-',color='steelblue')

