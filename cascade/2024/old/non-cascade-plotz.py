import numpy as np
import matplotlib.pyplot as pl
import pandas as pd
import os, socket
import glob


data_dir = '/media/xaber/extradrive1/crystalize_data/data-202403/20240315/20240315-150207'
data_dir = '/Users/peter/Public/data/20240315-150207/'
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
s2 = np.zeros([32,h_n_events])
ge = np.zeros([32,h_n_events])

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
		s2[:,aa_last:(aa_last+a)] = data["arr_6"]
		ge[:,aa_last:(aa_last+a)] = data["arr_7"]		
	aa_last = aa_last + a

asum = np.sum(aa,axis=0)
es1 = np.sum(s1,axis=0)
es2 = np.sum(s2,axis=0)
ce = np.sum(aa,axis=0)
rr = ce/es2

###
cut = (es1>10) & (es1<700) & (es2<2e4)
N = np.sum(cut)
print('cut keeps %d events'%N)
###

# pl.figure(11);pl.clf();
# pl.plot(es1[cut],es2[cut],'ko',markersize=3)

pl.figure(12);pl.clf();
pl.plot(es2[cut],s2f[cut],'ko',markersize=3)

pl.figure(13);pl.clf();
pl.plot(es2[cut],ce[cut],'ko',markersize=3)


beans = np.arange(-4,-1,0.1)
binc = (beans[0:-1] + beans[1:])/2
[cts,beans] = np.histogram(np.log10(rr[cut]),beans)
pl.figure(10);pl.clf();
pl.plot(binc,cts,'ko-',linewidth=0.5)
pl.xlabel('log10 (cathodeEcho/s2)')
pl.ylabel('cathodeEcho / s2')
#pl.ylim([1e-3,1e6])
#pl.xscale('lin')
#pl.yscale('log')
# pl.legend()
pl.title(data_dir[-16:-1])

pl.figure(9);pl.clf();
pl.plot(es2[cut],rr[cut],'ko',markersize=3)
pl.xlabel('s2')
pl.ylabel('cts')
#pl.ylim([1e-3,1e6])
#pl.xscale('lin')
pl.yscale('log')
# pl.legend()
pl.title(data_dir[-16:-1])
