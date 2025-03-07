import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os
import glob

csvfilez = glob.glob("/Users/peter/Desktop/*.csv")
nf = len(csvfilez)
# nn=0
# for ii in range(0,nf):
# 	vv=pd.read_csv(csvfilez[ii]) 
# 	nn+=vv.shape[0]
	
#vv=pd.read_csv('./20231208T2327.csv')


mk = np.array(['--',':'])
pl.figure(1);pl.clf()

for ii in range(0,nf):
# 	vv = pd.read_csv(csvfilez[nf-1-ii]) # read list backwards to get in chronological order
	vv = pd.read_csv(csvfilez[ii]) # read list 
	nn=vv.shape[0]
	lastn=0
	lastt=0
	t=np.zeros(nn)
	p=np.zeros(nn)
	T0=np.zeros(nn)
	T5=np.zeros(nn)
	T6=np.zeros(nn)
	SLM=np.zeros(nn)
	wb=np.zeros(nn)

	MM = csvfilez[0][6:8]
	DD = csvfilez[0][8:10]
	HH = csvfilez[0][11:13]
	mm = csvfilez[0][13:15]
	for n in range(0,vv.shape[0]):
		t[lastn+n]=vv.iat[n,2]+lastt;
		T0[lastn+n]=vv.iat[n,3];
		T5[lastn+n]=vv.iat[n,4];
		T6[lastn+n]=vv.iat[n,5];
		SLM[lastn+n]=vv.iat[n,7];
		p[lastn+n]=vv.iat[n,8];
		wb[lastn+n]=vv.iat[n,9];
	lastt=t[n]



	tt = t/3600
# 	ra = np.arange(0,1500,1)
# 	pl.plot(tt[ra],T0[ra],mk[ii],markersize=3,markerfacecolor='None')
	pl.plot(tt,T0,mk[ii],markersize=3,markerfacecolor='None',label='T0 [C]')
	pl.plot(tt,T5,mk[ii],markersize=3,markerfacecolor='None',label='T5 [C]')
	pl.plot(tt,T6,mk[ii],markersize=3,markerfacecolor='None',label='T6 [C]')
	pl.plot(tt,p,mk[ii],markersize=3,markerfacecolor='None',label='presure [Bar]')
	pl.plot(tt,SLM,mk[ii],markersize=3,markerfacecolor='None',label='flow [SLM]')
	pl.plot(tt,wb,mk[ii],markersize=3,markerfacecolor='None',label='top heater [W]')
	#pl.xlabel('hours since %s-%s @%s%s'%(MM,DD,HH,mm))
	pl.minorticks_on()
	pl.grid(True)
	pl.xlabel('hours')
	pl.ylabel('see legend')
	pl.legend()
	
if 0:	
	ws = 20 # window samples
	ww=np.lib.stride_tricks.sliding_window_view(wb,ws)
	ma_ww = ww.mean(axis=1)
	pl.plot(tt[0:-ws+1],ma_ww,'b.')

	
	pp=np.lib.stride_tricks.sliding_window_view(p,ws)
	ma_pp = pp.mean(axis=1)
# 	pl.plot(tt,p,'ko')
	pl.plot(tt[0:-ws+1],ma_pp,'k.')

	setpp = 1.17
	pl.plot(tt,np.ones(t.size)*(setpp*0.99),'k-',linewidth=1)
	pl.plot(tt,np.ones(t.size)*setpp,'k--',linewidth=1)
	pl.plot(tt,np.ones(t.size)*(setpp*1.01),'k-',linewidth=1)

#	pl.plot(tt,np.ones(t.size)*0.82,'r--')
	#pl.plot(tt,np.ones(t.size)*0.85,'r-')

	#pl.plot(np.array([7.5,31.5]),np.array([1,1])*0.805,'k--')

	pl.ylabel('pressure / bar')
	#pl.axis([-0.5,tt[-1]+0.5,0,1])
	#pl.figure(2);pl.clf()
	#pl.plot(tt,T0,'r-')
