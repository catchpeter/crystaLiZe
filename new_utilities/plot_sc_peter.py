import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os
import glob

csvfilez = glob.glob("./*.csv")
nf = len(csvfilez)
nn=0
for ii in range(0,nf):
	vv=pd.read_csv(csvfilez[ii]) 
	nn+=vv.shape[0]
	
#vv=pd.read_csv('./20231208T2327.csv')

t=np.zeros(nn)
p=np.zeros(nn)
T0=np.zeros(nn)
wb=np.zeros(nn)
lastn=0
lastt=0
for ii in range(0,nf):
	vv=pd.read_csv(csvfilez[nf-1-ii]) # read list backwards to get in chronological order
	MM = csvfilez[0][6:8]
	DD = csvfilez[0][8:10]
	HH = csvfilez[0][11:13]
	mm = csvfilez[0][13:15]
	for n in range(0,vv.shape[0]):
		t[lastn+n]=vv.iat[n,2]+lastt;
		T0[lastn+n]=vv.iat[n,3];
		p[lastn+n]=vv.iat[n,7];
		wb[lastn+n]=vv.iat[n,9];
	#lastn=n+timebetweenfilestartendplushowfaritran
	lastt=t[n]

#pl.plot(t,p,'ko')

if 1:
	if 1:
		tt = t/3600
		pl.xlabel('hours since %s-%s @%s%s'%(MM,DD,HH,mm))
	else:
		tt = t/3600 + (48+7.67)
		pl.xlabel('hours since Mon Dec 4 @0950')

	
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
