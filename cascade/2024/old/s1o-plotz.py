import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob
from matplotlib import rc
rc('font', **{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
rc('xtick', labelsize=14)
rc('ytick', labelsize=14)


data_dir_base = '/Users/peter/Public/data/'
data_dirs = np.array(['20240604-181209/','20240604-192140/','20240604-202159/','20240604-212218/','20240604-222237/','20240604-232256/'])

ky0=0
# it=0
N=0

for dd in data_dirs:
	data_dir = data_dir_base+dd
	aa_file_list = glob.glob(data_dir+"./aa-s1o/*v1.npz")
	print('found %d files'%len(aa_file_list))
	h_file = np.load(data_dir+"/compressed_filtered_data/headers.npz")
	h_array = h_file["arr_0"]
	h_n_events = int(np.floor(h_array.size/8))

	aa = np.zeros([32,h_n_events])
	ee = np.zeros([32,h_n_events])
	s1 = np.zeros([32,h_n_events])
	eec = np.zeros([h_n_events])
	ss = np.zeros([50000,h_n_events])
	tt = np.zeros([h_n_events,10])
	se = np.zeros([h_n_events,10])
	y0s = np.zeros([50000,h_n_events])

	aa_last = 0
	for aa_file in aa_file_list:
		#print('load file %s'%aa_file)
	#	np.savez(aa_file,aa,ee,s1,is1,is2,s2f)
		with np.load(aa_file) as data:
			a = data["arr_0"].shape[1]
			aa[:,aa_last:(aa_last+a)] = data["arr_0"]
			s1[:,aa_last:(aa_last+a)] = data["arr_1"]
			ee[:,aa_last:(aa_last+a)] = data["arr_2"]
			eec[aa_last:(aa_last+a)] = data["arr_3"]
			ss[:,aa_last:(aa_last+a)] = data["arr_4"]
			tt[aa_last:(aa_last+a),:] = data["arr_5"]
			se[aa_last:(aa_last+a),:] = data["arr_6"]
			y0s[:,aa_last:(aa_last+a)] = data["arr_7"][:,0:a] # kludge
		aa_last = aa_last + a


	asum = np.sum(aa,axis=0)
# 	ei = (asum.shape[0]-1) # end index
# 	ei = int(np.floor(ei/4)*4)
# 	a0 = ( asum[np.arange(0,ei,4)] )
# 	a1 = ( asum[np.arange(1,ei,4)] )
# 	a2 = ( asum[np.arange(2,ei,4)] )
# 	a3 = ( asum[np.arange(3,ei,4)] )
# # 
# 	ee0 = np.sum(ee[:,np.arange(0,ei,4)],axis=0)
# 	ee1 = np.sum(ee[:,np.arange(1,ei,4)],axis=0)
# 	ee2 = np.sum(ee[:,np.arange(2,ei,4)],axis=0)
# 	ee3 = np.sum(ee[:,np.arange(3,ei,4)],axis=0)
	eesum = np.sum(ee,axis=0)
	s1sum = np.sum(s1,axis=0)
	###
# 	cut = (tt[:,0]>0) & (np.sum(ee,axis=0)<100)#(s1sum>1000) & (s1sum<1e4) #(s1sum>1000)
# 	n = np.sum(cut)
# 	print('cut keeps %d events'%n)
# 	nz=np.nonzero(cut); nz=nz[0]

	###
# 	N = N+n
# 	ky0 = ky0 + np.sum(y0s[:,cut],axis=1) # stacked waveform
# 	if it==0:
# 		ktt = tt[cut]
# 		kse = se[cut]
# 	else:
# 		ktt = np.vstack((ktt,tt[cut])) #np.append(tt[cut],ktt)
# 		kse = np.vstack((kse,se[cut])) #np.append(se[cut],kse)
# 	it+=1
	
# dtt = 0.1
dtt = 0.002
t = np.arange(0.0001,100,dtt)


cut = (s1sum>5000) & (s1sum<10000) & (se[:,0]<100) & (s1sum/np.sum(ss,axis=0)>0.75) & (eec<5)# & (np.sum(y0s,axis=0)<2e5)
N = np.sum(cut)
print('cut keeps %d events'%N)
nz=np.nonzero(cut); nz=nz[0]


db=3
beans = np.arange(0,100,db)
beanc = (beans[1:]+beans[0:-1])/2
pl.figure(10);pl.clf();
[cts,beans] = np.histogram(np.concatenate([eesum]),beans); cts[0]=0
pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='k.')

pl.figure(1);pl.clf();
pl.semilogy(t,8*np.sum(y0s[:,cut],axis=1),label='stacked S1o') # ad-hoc scaling to match data point
#pl.plot(2/1e3*tt[cut,:],se[cut,:],'o',markeredgecolor='darkorange',markerfacecolor='None')

pl.plot(np.array([1,1])*2e-3*4930,np.array([0.1,1e5]),'k:')
pl.plot(np.array([1,1])*2e-3*(4930+3200),np.array([0.1,1e5]),'k:')
pl.plot(np.array([1,1])*2e-3*(4930+3200*2),np.array([0.1,1e5]),'k:')
pl.plot(np.array([1,1])*2e-3*(4930+3200*3),np.array([0.1,1e5]),'k:')
pl.plot(np.array([1,1])*2e-3*(4930+3200*4),np.array([0.1,1e5]),'k:')
pl.minorticks_on()

db=5
beans = np.arange(25,105,db)
beanc = (beans[1:]+beans[0:-1])/2
[cts,beans] = np.histogram(se[cut,0],beans); cts[0]=0; cts[1]=0; cts[2]=0
pl.errorbar(beanc,100*cts/db/N,yerr=np.sqrt(cts)/N,xerr=db/2,fmt='.',color='darkorange',label=('electrons/100 us'))

tl = np.arange(0.0001,5000,dtt)
pl.plot(tl[tl>70],900/tl[tl>70],'k-',linewidth=0.5,label='1/t')
pl.plot(tl[tl>70],9/tl[tl>70],'k-',linewidth=0.5,label='1/t')
pl.plot(tl[tl>70],170*np.exp(-(tl-70)/50)[tl>70],'k-',linewidth=0.5,label='1/t')


### cascade
ndp = np.array([np.sum(np.sum(aa[:,nz+1],axis=0)),np.sum(np.sum(aa[:,nz+2],axis=0)),np.sum(np.sum(aa[:,nz+3],axis=0))])
nde = np.array([np.sum(eec[nz+1]),np.sum(eec[nz+2]),np.sum(eec[nz+3])])

lastbit = np.sum(np.sum(y0s[40000:,cut],axis=1))*2/25 / N /(1/5) # approx hack with the *2/25 since I did not save y0spe variable
pl.errorbar(87.5,lastbit,xerr=12.5,yerr=np.sqrt(lastbit),fmt='o',color='darkblue')
pl.errorbar(np.array([550,1050,5050]),ndp/N,yerr=np.sqrt(ndp)/N,xerr=50,fmt='o',color='darkblue')
pl.errorbar(np.array([550,1050,5050]),nde/N,yerr=np.sqrt(nde)/N,xerr=50,fmt='o')

pl.xscale('log')
pl.legend(loc='best')
pl.xlabel('time ($\mu$s)',fontsize=14)
pl.ylabel('see legend',fontsize=14)
pl.xlim([1,6e3])
pl.ylim([1e-3,1e6])

# QC
# pl.figure(1);pl.hist(se[cut,0])
# pl.figure(2);pl.hist( s1sum[cut]/np.sum(ss[:,cut],axis=0) ,100)
# pl.figure(3);pl.hist(s1sum,np.arange(100,50e3,1e3))
# pl.figure(3);pl.hist(np.sum(y0s,axis=0),np.arange(1e3,5e5,1e4))
