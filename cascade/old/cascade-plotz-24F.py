import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob


# Run 24F
data_dir = '/Users/peter/Public/data/20240529-165141/' 
data_dir = '/Users/peter/Public/data/20240529-171140/' # LED
# data_dir = '/Users/peter/Public/data/20240529-174148/' # LED +power

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
ss = np.zeros([25000,h_n_events])
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


###
cut = (a0>1e3) & (a0<3e3)# & (s2f>0.65)
N = np.sum(cut)
print('cut keeps %d events'%N)
###



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
dpt = np.array([0.0001, 0.5, 1.0, 5.0]) # trigger cascade (set by hardware)
Rce = 0.0063 # small-s2 limit of ratio S2ce/S2 -- should triple check
#Rce=0.01
### define the number of detected photons. subtract the number of phd identified as single e-
detp = np.array([ np.sum(a0[cut])/Rce , np.sum(a1[cut]-ee1[cut]) , np.sum(a2[cut]-ee2[cut]) , np.sum(a3[cut]-ee3[cut]) ]) 
#detp = np.array([ np.sum(a0[cut])/0.01 , np.sum(a1[cut]) , np.sum(a2[cut]) , np.sum(a3[cut]) ]) 
detp = detp/N

### number of detected electrons in the cascades
dete = np.array([ detp[0]/se , np.sum(ee1>0)/N , np.sum(ee2>0)/N , np.sum(ee3>0)/N ]) # this is too simplistic, need to improve lower-level algorithm

dtt = 0.1
# dtt = 1
tt = np.arange(0.0001,1000,dtt)



if (data_dir[-12:-1]=='0529-165141'):
	fit = 50/tt # photons
	elbl = 3.0/tt
if (data_dir[-12:-1]=='0529-171140'):
	fit = 500/tt # photons
	elbl = 4.0/tt
if (data_dir[-12:-1]=='0529-174148'):
	fit = 500/tt # photons
	elbl = 4.0/tt



if 1:
	pbw = 0.1 # plot bin width, fixed to cascade event window!
	pl.figure(1);pl.clf()
	pl.errorbar(dpt,detp,yerr=np.sqrt(detp*N)/N,fmt='o',color='steelblue',markersize=5,markerfacecolor='white',label=(r'photon signal'))
	#tau = 3.0; pl.plot(tt,120*np.exp(-tt/tau),'r-',label=(r'$\tau=%1.1d$ ms'%tau))

	# photons	
	pl.step(tt,fit,'-',linewidth=0.5,color='steelblue',label=(r'$1/t$ ($\gamma$)'))
	delayed_p_3ms = np.sum(fit[tt>3])*dtt/pbw
	print('delayed photons in 3-1000 ms window: %1.2f%% (%1.0f total)'% (delayed_p_3ms/(detp[0])*100,delayed_p_3ms) )

	# electrons
	pl.errorbar(dpt,dete,yerr=np.sqrt(dete*N)/N,fmt='o',color='darkorange',markersize=5,markerfacecolor='white',label=(r'electron signal'))
	pl.step(tt,elbl,'-',linewidth=0.5,color='darkorange',label=(r'$1/t$ (e-)'))
	lbl_delayed_e_3ms = np.sum(elbl[tt>3])*dtt/pbw
	print('(LBL) delayed electrons in 3-1000 ms window: %1.2f%% (%1.0f total)'% (lbl_delayed_e_3ms/dete[0]*100,lbl_delayed_e_3ms) )

	pl.plot(3*np.array([1,1]),np.array([1e-3,fit[tt>3][0]]),'k:')
	
	pl.xlabel('time (ms)')
	pl.ylabel('(e or p) / 0.1 ms')
	pl.ylim([1e-3,1e6])
	pl.xscale('log')
	pl.yscale('log')
	pl.legend()
	pl.title(data_dir[-16:-1])




