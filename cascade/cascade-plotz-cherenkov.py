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

	data_folders = np.array(['20250311-174524','20250311-220055','20250312-072556','20250312-111635']) # 165 K no Xe, sipms at 50 V, cerenkov trigger

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


af = np.zeros((5,data_folders.shape[0]))

ii=0
data_dir = '/Users/peter/Public/data/'+data_folders[ii]+'/'
aa_file_list = glob.glob(data_dir+"./aa/*v2.npz")
print('looking in: %s'%data_folders[ii])
# print('found %d files'%len(aa_file_list))
print('loading header file')
h_file = np.load(data_dir+"/compressed_filtered_data/headers.npz")
h_array = h_file["arr_0"]
h_n_events = int(np.floor(h_array.size/8))

print('warning: invented a big number for h_n_events, may need to check')
h_n_events = 30000

# check data_size
check_data = np.load(aa_file_list[0])
nch = check_data["arr_2"].shape[0]
nss = check_data["arr_3"].shape[2]
aa = np.zeros([nch,h_n_events])
ss = np.zeros([nch,h_n_events,nss])
s1 = np.zeros([nch,h_n_events])
s1ap = np.zeros([nch,h_n_events])
aa_last = 0

for ii in range(0,data_folders.shape[0]):
	data_dir = '/Users/peter/Public/data/'+data_folders[ii]+'/'
	aa_file_list = glob.glob(data_dir+"./aa/*v2.npz")
	print('looking in: %s'%data_folders[ii])
	print('found %d files'%len(aa_file_list))
	for aa_file in aa_file_list:
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
# 	if    (data_folders[0][-15:-7]=='20250227') \
# 		| (data_folders[0][-15:-7]=='20250305') \
# 		| (data_folders[0][-15:-7]=='20250306') \
# 		| (data_folders[0][-15:-7]=='20250307'): # then looking at no-xenon data
if (int(data_folders[0][-15:-7])>20250227) :
	print('identified as no Xe data')
	xe=0
else:
	s1bot = s1bot*(5/3)
	print('\n*** Assuming alphas in Xe (ADC saturated): accounting for PMT ADC saturation factor 5/3, obtained from afterpulse size')
	xe=1

s1bot_ap = (s1ap[pmt_ch,np.arange(0,ei,4)])
s1top = np.sum(s1[0:pmt_ch,np.arange(0,ei,4)],axis=0)
s1top_ap = np.sum(s1ap[0:pmt_ch,np.arange(0,ei,4)],axis=0)
s1top = s1top + s1top_ap
s10 = s1top + s1bot
s1tba = (s1top-s1bot)/(s1top+s1bot)

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

	if 0:
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
cut = (s10>300) & (s10<500) \
	& (s1bot_ap<10)

				

nz=np.nonzero(cut); nz = nz[0]
N = np.sum(cut)
print('\n*** cut keeps %d events ***\n'%N)
###

		



### aggregate the data
psipm = 0

# 	ct = np.array([0.01, 0.075, 0.5+0.05, 1.0+0.05, 5.0+0.05]) # trigger cascade (set by hardware) -- used for ever all earlier data
ct = np.array([0.01, 0.075, 0.2+0.05, 0.5+0.05, 1.0+0.05]) # trigger cascade (set by hardware) -- afternoon Feb 20 +
if (data_folders[0][-6:]=='094512'):
	ct = np.array([0.01, 0.075, 0.3+0.05, 0.5+0.05, 1.0+0.05]) # trigger cascade (set by hardware) -- briefly used before noon on Feb 20

### define the number of detected photons. subtract the number of phd identified as single e-
dpt = np.array([ np.sum(s1top[cut])/N , np.sum(a0t[cut])*2/N , np.sum(a1t[cut])/N , np.sum(a2t[cut])/N , np.sum(a3t[cut])/N ])
dpb = np.array([ np.sum(s1bot[cut])/N , np.sum(a0b[cut])*2/N , np.sum(a1b[cut])/N , np.sum(a2b[cut])/N , np.sum(a3b[cut])/N ])
# 	af[:,ii] = detp#/triggerS1

colr='blue'
if 1:
	pbw = 0.1 # plot bin width, fixed to cascade event window!
	pl.figure(8);pl.clf(); ax = pl.gca()
	pl.errorbar(ct[1:],dpb[1:],yerr=np.sqrt(dpb[1:]*N)/N,xerr=np.array([0.025,0.05,0.05,0.05]),fmt='o',color=colr,markersize=5,markerfacecolor='white',lw=0.5,label='R8778 PMT')
	if psipm:
		pl.errorbar(ct[1:],dpt[1:],yerr=np.sqrt(dpt[1:]*N)/N,xerr=np.array([0.025,0.05,0.05,0.05]),fmt='s',color='grey',markersize=5,markerfacecolor='white',lw=0.5,label='S13371 SiPM')
	
	pl.plot(ct[0],dpb[0],'o',color=colr,markersize=7,markerfacecolor='white')
	if psipm:
		pl.plot(ct[0],dpt[0],'s',color='grey',markersize=7,markerfacecolor='white')

	pl.text(0.012,dpb[0],(r'Cherenkov pulse $\bar{a}$'),color='k',verticalalignment='center',size=14)

	(a, b, sigma_a, sigma_b) = linfit(np.log10(ct[1:4]),np.log10(dpb[1:4]))
	pl.plot(tt,(10**a)*tt**b,'--',color='blue',linewidth=0.5)#,label='$at^{-b}$') 
	print(b)
	print(1/(10**a/dpb[0]))

	# UCLA:
	#pl.plot(tt,(dpb[0]/5313)*tt**-1.097,'--',color=colr,linewidth=1,label='Xe expectation') # fitting last 3 points, AP<300
	#pl.text(0.012,dpb[0]/70,('delayed photons $d_p(t)$'),rotation=-22,color='k',verticalalignment='center')
	
	# paper:
	if 0:
		pl.plot(tt[0:],(dpb[0]/avals[1])*tt[0:]**bvals[1],'--',color='red',linewidth=0.5)
		pl.plot(tt[0:],(dpb[0]/avals[2])*tt[0:]**bvals[2],'--',color='orange',linewidth=0.5)
		pl.plot(tt[0:],(dpb[0]/avals[3])*tt[0:]**bvals[3],'--',color='gold',linewidth=0.5)
		pl.plot(tt[0:],(dpb[0]/avals[4])*tt[0:]**bvals[4],'--',color='green',linewidth=0.5)
		pl.plot(tt[0:],(dpb[0]/avals[5])*tt[0:]**bvals[5],'--',color='blue',linewidth=0.5)
		pl.plot(tt[0:],(dpb[0]/avals[6])*tt[0:]**bvals[6],'--',color='indigo',linewidth=0.5)
	if 0:
		sysc='k'
		pl.plot(tt[0:],(dpb[0]/avals[1])*tt[0:]**bvals[1],':',color=sysc,linewidth=0.5)
		pl.plot(tt[0:],(dpb[0]/avals[2])*tt[0:]**bvals[2],':',color=sysc,linewidth=0.5)
		pl.plot(tt[0:],(dpb[0]/avals[3])*tt[0:]**bvals[3],':',color=sysc,linewidth=0.5)
		pl.plot(tt[0:],(dpb[0]/avals[4])*tt[0:]**bvals[4],':',color=sysc,linewidth=0.5)
		pl.plot(tt[0:],(dpb[0]/avals[5])*tt[0:]**bvals[5],':',color=sysc,linewidth=0.5)
# 		pl.plot(tt[0:],(dpb[0]/avals[6])*tt[0:]**bvals[6],':',color=sysc,linewidth=0.5)
	
	pl.text(0.012,dpb[0]/50,('$d_p(t)=0.09t^{-1.03}$'),rotation=-22,color='k',verticalalignment='center',size=14)

	if 1:
		pl.text(0.012,0.3,('random photon background'),color='k',verticalalignment='center',size=14)
		sig_t = 0.07; mu_t = 0.43
		sig_b = 0.06; mu_b = 0.19
		pl.plot(np.array([1e-2,2]),np.ones(2)*mu_b,':',color='grey',linewidth=1)
		ax.add_patch(Rectangle((1e-2, mu_b-sig_b), 2, sig_b*2,facecolor=colr,alpha=0.25, edgecolor='None'))
		if psipm:
			pl.plot(np.array([1e-2,2]),np.ones(2)*mu_t,':',color='grey',linewidth=1)
			ax.add_patch(Rectangle((1e-2, mu_t-sig_t), 2, sig_t*2,color='grey',alpha=0.25,edgecolor='None'))
		


	if 0: # top sipm array
		(a, b, sigma_a, sigma_b) = linfit(np.log10(ct[1:4]),np.log10(dpb[1:4]))
		print(b)
		print(1/(10**a/dpb[0]))
		pl.plot(tt,(10**a)*tt**b,'-',color='grey',linewidth=1)#,label='fit to sipms') 

	pl.xlabel('time (ms)',size=14)
	pl.ylabel('photon counts',size=14)
# 		pl.ylim([1e-5,2])
	pl.xscale('log')
	pl.yscale('log')
	pl.legend(fontsize=14)
# 	pl.title(data_dir[-16:-1])
	pl.ylim([0.1,1e4])

pl.savefig('fig2.png',dpi=300)

