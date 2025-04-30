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



if 1: # PTFE, Xe liquid


# 	data_folders = np.array(['20250424-091853','20250424-092351']) # same but lower stats
	data_folders = np.array(['20250424-101744','20250424-094903'])
	labl = np.array(['S13371 SiPM','S13370 SiPM','','']) # 50 V and 47.5 V
	mrk = np.array(['s','d'])
	colorz = np.array(['gray','green','olivedrab','steelblue','sienna'])

	data_folders = np.array(['20250424-094149']); # line triggers

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

def dp(a,b,t,A,tp,toff):
	# Tyler Anderson's function from his thesis
	# a: amplitude coefficient
	# b: power law exponent
	# t: time base (assumed ms here)
	# A: progenitor pulse counts (photons detected)
	# tp: time of the pulse
	# toff: offset time in which the function is not calculated
	wf = a*A*((t-tp)/toff)**b
	wf[t<(tp+toff)] = 0	
	return (wf)

dtt = 0.01
tt = np.arange(0.012,1.5,dtt)
# tt = np.arange(0.012,100,dtt)
fitdp = np.array([1.0e-4*tt**-1,1.2e-4*tt**-1])



af = np.zeros((5,data_folders.shape[0]))
print('NOTE: code assumes windowless sipms are channels 0,5,10,15 (physical channels +1)')
pl.figure(8);pl.clf();ax=pl.gca()

for ii in range(0,data_folders.shape[0]):
#for ii in range(0,1):
	if ii==0:
		ww=8
	elif ii==1:
		ww=5
	
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
# 		print(aa_file)
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
		
		lr = np.ones(nch)*6  # lower fit range for finding gain
		lr[5] = 10
		lr[10] = 10
		ur = 50 # upper fit range for finding gain
		if (any(cts0123[int(lr[ch]):ur])):
			gains[ch] = np.average(beanc[int(lr[ch]):ur],axis=0,weights=cts0123[int(lr[ch]):ur]) #
		else:
			print('not enough counts for average')
		gains[gains==0]=1e-6
		
		if ((ch==5) | (ch==10) | (ch==8)): 
			pl.figure(6);pl.clf();
			pl.errorbar(beanc,cts0,yerr=np.sqrt(cts0),xerr=db/2,fmt='k.')
			pl.errorbar(beanc,cts1,yerr=np.sqrt(cts1),xerr=db/2,fmt='r.')
			pl.errorbar(beanc,cts2,yerr=np.sqrt(cts2),xerr=db/2,fmt='b.')
			pl.errorbar(beanc,cts3,yerr=np.sqrt(cts3),xerr=db/2,fmt='g.')

			pl.errorbar(beanc,cts0123,yerr=np.sqrt(cts0123),xerr=db/2,fmt='.',color='gray')
	
			pl.plot(np.ones(2)*gains[ch],np.array([1,max(cts0123)]),'-',color='powderblue')
			pl.plot(np.ones(2)*beanc[int(lr[ch])],np.array([1,max(cts0123)]),':',color='powderblue')
			pl.plot(np.ones(2)*beanc[ur],np.array([1,max(cts0123)]),':',color='powderblue')
				

			pl.yscale('log')
			pl.title('ch %d, gain %2.1f'%(ch,gains[ch]))
			pl.show();pl.pause(0.1)	
			input('press any key')
			
		s1[ch,:] = s1[ch,:] / gains[ch]
# 		s1ap[ch,:] = s1ap[ch,:] / gains[ch]
# 		aa[ch,:] = aa[ch,:] / gains[ch]

#		input('paused...')

	s10 = s1[:,np.arange(0,ei,4)]
# 	sst = ss*(ss>lr) # threshold
	# sum all the pulses, or put a threshold to only sum spe
	ss0 = np.sum( ss[:,np.arange(0,ei,4),:] ,axis=2 )
	ss1 = np.sum( ss[:,np.arange(1,ei,4),:] ,axis=2 )
	ss2 = np.sum( ss[:,np.arange(2,ei,4),:] ,axis=2 )
	ss3 = np.sum( ss[:,np.arange(3,ei,4),:] ,axis=2 )


	if 0: # check spectrum
		pl.figure(10);#pl.clf();
		db=50
		beans = np.arange(0,5e3,db)
		beanc = (beans[1:]+beans[0:-1])/2
		[cts,beans] = np.histogram(s10[ww,:],beans); cts[0]=0
		pl.errorbar(beanc,cts,yerr=np.sqrt(cts),xerr=db/2,fmt='c+',label='sipm top')
# 		pl.title(data_folders[ii])
		pl.xlabel('phd')
# 		pl.legend()
		pl.minorticks_on()

	###
	cut = (s10[ww,:]>2000) & (s10[ww,:]<4000) \
		& (ss1[ww,:]<s10[ww,:]) \
		& (ss2[ww,:]<s10[ww,:]) \
		& (ss3[ww,:]<s10[ww,:])
# 	print('*** hey - add a cut for pulses >> 1-2 spe, within a delayed window?')	

	if ((data_folders[0][-6:])=='094149'): # the line trigger dataset
		ww=5
		cut = (s10[ww,:]>0) & (s10[ww,:]<100) \
			& (ss0[ww,:]<100) \
			& (ss1[ww,:]<100) \
			& (ss2[ww,:]<100) \
			& (ss3[ww,:]<100)

	nz=np.nonzero(cut); nz = nz[0]
	N = np.sum(cut)
	print('\n*** cut keeps %d events ***\n'%N)

	### aggregate the data
	ct = np.array([0.01, 0.075, 0.2+0.05, 0.5+0.05, 1.0+0.05]) # trigger cascade (set by hardware) -- afternoon Feb 20 +

	### define the number of detected photons. subtract the number of phd identified as single e-
	dpb = np.array([ np.sum(s10[ww,cut])/N , np.sum(ss0[ww,cut]/gains[ww])*2/N , np.sum(ss1[ww,cut]/gains[ww])/N , np.sum(ss2[ww,cut]/gains[ww])/N , np.sum(ss3[ww,cut]/gains[ww])/N ])
	print(dpb)
	if 1:
		pbw = 0.1 # plot bin width, fixed to cascade event window!
		pl.figure(8)#;pl.clf();ax=pl.gca()
		pl.plot(ct[0],dpb[0],mrk[ii],color=colorz[ii],markersize=7,markerfacecolor='None')
		pl.errorbar(ct[1:],dpb[1:],yerr=np.sqrt(dpb[1:]*N)/N,xerr=np.array([0.025,0.05,0.05,0.05]),fmt=mrk[ii],color=colorz[ii],markersize=5,markerfacecolor='white',lw=0.5,label=labl[ii])
		
		if ii==1:
			pl.text(0.012,dpb[0],(r'Xe scintillation pulse $\bar{a}$'),color='k',verticalalignment='center',fontsize=14)
			pl.text(0.012,dpb[0]/30,('delayed photons $d_p(t)$'),rotation=-27,color='k',verticalalignment='center',fontsize=14)

		
		(a, b, sigma_a, sigma_b) = linfit(np.log10(ct[1:4]),np.log10(dpb[1:4]))
		print(b)
		tmp = 1/(10**a/dpb[0])
		print(1/(10**a/dpb[0]))
# 		pl.plot(tt[0:],(10**a)*tt[0:]**b,'-',color=colorz[ii],linewidth=1)
		pl.plot(tt[0:],dpb[0]/tmp*tt[0:]**b,':',color=colorz[ii],linewidth=1)


		pl.text(0.012,0.4,('random photon background'),color='green',verticalalignment='center',fontsize=14)		
		sig = 0.04; mu_b = 0.65;
		pl.plot(np.array([1e-2,2]),np.ones(2)*mu_b,':',color='grey',linewidth=1)
		ax.add_patch(Rectangle((1e-2, mu_b-sig), 2, sig*2,facecolor=colorz[ii],alpha=0.25,edgecolor='None'))
# 		ax.add_patch(Rectangle((1e-2, mu_t-sig), 2, sig*2,facecolor=colorz[ii],alpha=0.5,edgecolor='None'))
# 		pl.plot(np.array([1e-2,2]),np.ones(2)*mu_t,':',color='k',linewidth=1) 
		

# 		pl.plot(tt,1e-4*dp(0.42/40,-1.3,tt,2700,0.01,40),'r:') # is there a missing factor x1/40 in the normalization?
		
		pl.xlabel('time (ms)',fontsize=14)
		pl.ylabel('photon counts',fontsize=14)
# 		pl.ylim([1e-5,2])
		pl.xscale('log')
		pl.yscale('log')
		pl.legend(fontsize=14)
		pl.ylim([0.1,1e4])

pl.savefig('fig3.png',dpi=300)




ttt = np.arange(1e-3,5e3,1e-3) # ms
ws = np.zeros(len(ttt))
ws2 = np.zeros(len(ttt))
pgs = np.array([7e5,5e5,3e5,8e5])
tps = np.array([0,1000,1500,2700])
labz = np.array(['photon pulse (counts)','','',''])
pl.figure(9);pl.clf();ax=pl.gca()
for i in range(len(pgs)):
	pl.plot(tps[i],pgs[i],'ro',label=labz[i])
	ws = ws + dp(0.42,-1.3,ttt,pgs[i],tps[i],40)
	#ws2 = ws2 + dp(1.6,-1.3,ttt,pgs[i],tps[i],15)
	ws2 = ws2 + dp(0.42/2.56,-0.15,ttt,pgs[i],tps[i],40)
pl.plot(ttt,ws+3400,'b-')
pl.plot(ttt,ws2+3400,'k--',lw=0.5)

pl.yscale('log')
pl.ylim([3e3,1e6])
pl.xlabel('time (ms)')
pl.ylabel('photon rate (Hz)')
pl.legend()
#w0 = dp(0.2,-1.3,ttt,1e6,0,40)



# wf0 = np.zeros(len(ttt))
# trig=0.01; sz = 1000*1e3
# wf0[ttt>=(trig+40000e-3)] = 10*sz/6086*(ttt[ttt>=(trig+40000e-3)])**-1.298
# cp = 1000000
# wf1 = np.append(np.zeros(cp),wf0[0:-cp])
# cp = 2000000
# wf2 = np.append(np.zeros(cp),wf0[0:-cp])
# pl.semilogy(ttt,(wf0+wf1+wf2)*1e4 + 3410,'k.',markersize=3)
