import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob


# 24B
data_dir = '/Users/peter/Public/data/20240314-190415/' # first good cascade data, 0.5,1.0,5.0 ms

# 24C
#data_dir = '/Users/peter/Public/data/20240319-150343/' # after open to air, 0.5,1.0,5.0 ms

# 24D
#data_dir = '/Users/peter/Public/data/20240411-103253/' # after 100C bake overnight 0.5,1.0,5.0 ms

# 24G
#data_dir = '/Users/peter/Public/data/20240604-163805/' # regular cascade
#data_dir = '/Users/peter/Public/data/20240604-181209/' # regular cascade, S1 trigger attempt

# 24H
data_dir = '/Users/peter/Public/data/20240611-173543/' # 133Ba cascade
#data_dir = '/Users/peter/Public/data/20240611-180620/' # 133Ba cascade
#data_dir = '/Users/peter/Public/data/20240612-173658/' # 133Ba+LED0.5us cascade	
data_dir = '/Users/peter/Public/data/20240613-155550/' # LED0.5us cascade	

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
a0 = ( asum[np.arange(0,ei,4)] )
a1 = ( asum[np.arange(1,ei,4)] )
a2 = ( asum[np.arange(2,ei,4)] )
a3 = ( asum[np.arange(3,ei,4)] )
s2cf = s2cf[np.arange(0,ei,4)]

if 0: # orignial method
	ee0 = np.sum(ee[:,np.arange(0,ei,4),:],axis=0)
	ee1 = np.sum(ee[:,np.arange(1,ei,4),:],axis=0)
	ee2 = np.sum(ee[:,np.arange(2,ei,4),:],axis=0)
	ee3 = np.sum(ee[:,np.arange(3,ei,4),:],axis=0)
else: # 17 June 2100 hrs method
	coin = np.sum(ee[:,:,:]>(1/3),axis=0)
	ees = np.sum(ee[:,:,:],axis=0)
	ees = ees * (coin>13)
	ee1 = ees[np.arange(1,ei,4),:]
	ee2 = ees[np.arange(2,ei,4),:]
	ee3 = ees[np.arange(3,ei,4),:]
	

dtt = 0.1
tt = np.arange(0.001,1000,dtt)

###
losat=1
if losat:
	cut = (a0>1e4) & (a0<1.5e5)#& (np.sum(ee1[:,:],axis=1)<450)
	fit = 50/tt # photons
else:
	cut = (a0>1e4) & (a0<2.5e5)
	fit = 150/tt # photons

if (data_dir[-12:-1]=='0604-181209'):
	cut = (a0>100) & (a0<1e4)
	fit = 2.5/tt # photons

nz=np.nonzero(cut); nz = nz[0]
N = np.sum(cut)
print('cut keeps %d events'%N)
###


if (data_dir[-12:-1]=='0314-190415'):
	if losat:
		pwr=1.5;elbl = 0.025*tt**-pwr
	else:
		pwr=1.5;elbl = 0.04*tt**-pwr
if (data_dir[-12:-1]=='0319-150343'):
	if losat:
		pwr=1.5;elbl = 0.04*tt**-pwr
	else:
		pwr=1;elbl = 0.1*tt**-pwr
if (data_dir[-12:-1]=='0411-103253'):
	if losat:
		pwr=1.5;elbl = 0.04*tt**-pwr
	else:
		pwr=0.5;elbl = 4.0*tt**-pwr	
if (data_dir[-12:-1]=='0604-181209'):
	pwr=1.5;elbl = 0.01*tt**-pwr
if (data_dir[-12:-1]=='0604-163805'):
	pwr=1.;elbl = 0.25*tt**-pwr
if (data_dir[-12:-1]=='0611-173543'):
	pwr=0.5;elbl = 4*tt**-pwr
if (data_dir[-12:-1]=='0612-093213'):
	pwr=0.5;elbl = 4*tt**-pwr
if (data_dir[-12:-1]=='0612-145819'):
	pwr=0.5;elbl = 4*tt**-pwr
if (data_dir[-12:-1]=='0612-173658'):
	pwr=0.5;elbl = 4*tt**-pwr
if (data_dir[-12:-1]=='0613-155550'):
	pwr=0.5;elbl = 4*tt**-pwr
#	pwr=0.5;elbl = 0.08*tt**-pwr
if (data_dir[-12:-1]=='0613-154900'):
	pwr=0.5;elbl = 4*tt**-pwr



db=3
beans = np.arange(0,100,db)
beanc = (beans[1:]+beans[0:-1])/2
pl.figure(10);pl.clf();
[cts,beans] = np.histogram(np.concatenate([ee1[cut,:],ee2[cut,:],ee3[cut,:]]),beans); cts[0]=0
pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='k.')
se = 25
#se = np.average(beanc,axis=0,weights=cts) # over-estimate

### aggregate the data
dpt = np.array([0.01, 0.5, 1.0, 5.0]) # trigger cascade (set by hardware)
# Rce = 0.0063 # small-s2 limit of ratio S2ce/S2 -- should triple check
### define the number of detected photons. subtract the number of phd identified as single e-
bfrac = 0.33
detp = np.array([ np.sum(a0[cut])/N/bfrac , np.sum(a1[cut]-np.sum(ee1[cut,:],axis=1))/N , np.sum(a2[cut]-np.sum(ee2[cut,:],axis=1))/N , np.sum(a3[cut]-np.sum(ee3[cut,:],axis=1))/N ])
#if (data_dir[-12:-1]=='0613-155550'): # adjust for LED-only case
#	detp[0] = detp[0] * 6.9 # 344820

### number of detected electrons in the cascades
dete = np.array([ detp[0] , np.sum(ee1[cut,:])/N , np.sum(ee2[cut,:])/N , np.sum(ee3[cut,:])/N ]) / se 

ptrainf = detp[1]/detp[0]*100
print("photon train fraction in 100 us window at 500 us delayt is %1.4f"%ptrainf)



if 1:
	pbw = 0.1 # plot bin width, fixed to cascade event window!
	pl.figure(8);pl.clf()
	pl.errorbar(dpt,detp/detp[0],yerr=np.sqrt(detp*N)/N/detp[0],fmt='o',color='steelblue',markersize=5,markerfacecolor='white',label=(r'photon signal'))
	#tau = 3.0; pl.plot(tt,120*np.exp(-tt/tau),'r-',label=(r'$\tau=%1.1d$ ms'%tau))

	# photons	
	fitn = 2e-4/tt
	pl.plot(tt,fitn,'-',linewidth=0.5,color='steelblue',label=(r'$t^{-1}$ ($\gamma$)'))
	if (data_dir[-12:-1]=='0612-173658'):
		pl.plot(tt,5*2e-4/tt,'--',linewidth=0.5,color='steelblue')#,label=(r'$t^{-1}$ ($\gamma$)'))
	if (data_dir[-12:-1]=='0613-155550'):
		pl.plot(tt,35*2e-4/tt,'--',linewidth=0.5,color='steelblue')#,label=(r'$t^{-1}$ ($\gamma$)'))
		
# 	pl.step(tt,fit,'-',linewidth=0.5,color='steelblue',label=(r'$1/t$ ($\gamma$)'))
#	delayed_p_3ms = np.sum(fit[tt>3])*dtt/pbw
#	print('delayed photons in 3-1000 ms window: %1.2f%% (%1.0f total)'% (delayed_p_3ms/(detp[0])*100,delayed_p_3ms) )

	# electrons
	pl.errorbar(dpt,dete/detp[0],yerr=np.sqrt(dete*N)/N/detp[0],fmt='o',color='darkorange',markersize=5,markerfacecolor='white',label=(r'electron signal'))
	if (data_dir[-12:-1]=='0612-173658'):
		pl.plot(tt,3*elbl/detp[0],'--',linewidth=0.5,color='darkorange')#,label=(r'$t^{-%1.1f}$ (e-)'%pwr))
	if (data_dir[-12:-1]=='0613-155550'):
		#pl.plot(tt,1/50*elbl/detp[0],'--',linewidth=0.5,color='darkorange')#,label=(r'$t^{-%1.1f}$ (e-)'%pwr))
		print('skip curve plot for this dataset')
	else:
		pl.plot(tt,elbl/detp[0],'-',linewidth=0.5,color='darkorange',label=(r'$t^{-%1.1f}$ (e-)'%pwr))

# 	pl.step(tt,elbl,'-',linewidth=0.5,color='darkorange',label=(r'$1/t$ (e-)'))
#	lbl_delayed_e_3ms = np.sum(elbl[tt>3])*dtt/pbw
#	print('delayed electrons in 3-1000 ms window: %1.2f%% (%1.0f total)'% (lbl_delayed_e_3ms/dete[0]*100,lbl_delayed_e_3ms) )

# 	pl.plot(3*np.array([1,1]),np.array([1e-3,fit[tt>3][0]]),'k:')

	pl.plot(np.array([1e-3,1e-1]),np.array([1,1]),'k:',label='progenitor window')
	
	pl.xlabel('time (ms)')
	pl.ylabel('(e- or $\gamma$) / 0.1 ms')
	pl.ylim([1e-9,2])
	pl.xscale('log')
	pl.yscale('log')
	pl.legend()
	pl.title(data_dir[-16:-1])




#pl.figure(90);pl.plot(np.sum(ee[:,:,0],axis=0),np.sum(ee[:,:,0]>0.3,axis=0),'k.')

if 0:
	pbw = 0.1 # plot bin width, fixed to cascade event window!
	pl.figure(1);pl.clf()
	pl.errorbar(dpt,detp,yerr=np.sqrt(detp*N)/N,fmt='o',color='steelblue',markersize=5,markerfacecolor='white',label=(r'photon signal'))
	#tau = 3.0; pl.plot(tt,120*np.exp(-tt/tau),'r-',label=(r'$\tau=%1.1d$ ms'%tau))

	# photons	
	pl.plot(tt,fit,'-',linewidth=0.5,color='steelblue',label=(r'$t^{-1}$ ($\gamma$)'))
# 	pl.step(tt,fit,'-',linewidth=0.5,color='steelblue',label=(r'$1/t$ ($\gamma$)'))
	delayed_p_3ms = np.sum(fit[tt>3])*dtt/pbw
	print('delayed photons in 3-1000 ms window: %1.2f%% (%1.0f total)'% (delayed_p_3ms/(detp[0])*100,delayed_p_3ms) )

	# electrons
	pl.errorbar(dpt,dete,yerr=np.sqrt(dete*N)/N,fmt='o',color='darkorange',markersize=5,markerfacecolor='white',label=(r'electron signal'))
	pl.plot(tt,elbl,'-',linewidth=0.5,color='darkorange',label=(r'$t^{-%1.1f}$ (e-)'%pwr))
# 	pl.step(tt,elbl,'-',linewidth=0.5,color='darkorange',label=(r'$1/t$ (e-)'))
	lbl_delayed_e_3ms = np.sum(elbl[tt>3])*dtt/pbw
	print('delayed electrons in 3-1000 ms window: %1.2f%% (%1.0f total)'% (lbl_delayed_e_3ms/dete[0]*100,lbl_delayed_e_3ms) )

# 	pl.plot(3*np.array([1,1]),np.array([1e-3,fit[tt>3][0]]),'k:')

	pl.plot(np.array([1e-3,1e-1]),np.array([1,1]),'k:',label='progenitor window')
	
	pl.xlabel('time (ms)')
	pl.ylabel('(e- or $\gamma$) / 0.1 ms')
	pl.ylim([1e-3,1e6])
	pl.xscale('log')
	pl.yscale('log')
	pl.legend()
	pl.title(data_dir[-16:-1])



