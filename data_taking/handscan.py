import numpy as np
import sys, os
sys.path.append(os.path.join(os.path.dirname(sys.path[0]), "pulse_finding_scripts"))
from rq_generate_32ch import make_rq


with open(sys.path[0]+"/handscan_path.txt", 'r') as path:
    data_dir = path.readline().strip() #only read first line of path.txt

#read RQ
listrq = np.load(data_dir+'rq.npz')

n_events = listrq['n_events'][()]
n_pulses = listrq['n_pulses']
n_s1 = listrq['n_s1']
n_s2 = listrq['n_s2']
s1_before_s2 = listrq['s1_before_s2']
p_area = listrq['p_area']
p_class = listrq['p_class']
drift_Time = listrq['drift_Time']
drift_Time_AS = listrq['drift_Time_AS']
p_max_height = listrq['p_max_height']
p_min_height = listrq['p_min_height']
p_width = listrq['p_width']
p_afs_2l = listrq['p_afs_2l']
p_afs_50 = listrq['p_afs_50']
p_area_ch = listrq['p_area_ch']
p_area_ch_frac = listrq['p_area_ch_frac']
p_area_top = listrq['p_area_top']
p_area_bottom = listrq['p_area_bottom']
p_tba = listrq['p_tba']
p_start = listrq['p_start']
p_end = listrq['p_end']
sum_s1_area = listrq['sum_s1_area']
sum_s2_area = listrq['sum_s2_area']
center_top_x = listrq['center_top_x']
center_top_y = listrq['center_top_y']
center_bot_x = listrq['center_bot_x']
center_bot_y = listrq['center_bot_y']

d_between_SiPM_center_x = 1.23 # cm
d_between_SiPM_center_y = 1.14 # cm
center_bot_x_d = center_bot_x * d_between_SiPM_center_x/2
center_bot_y_d = center_bot_y * d_between_SiPM_center_y/2
center_top_x_d = center_top_x * d_between_SiPM_center_x/2
center_top_y_d = center_top_y * d_between_SiPM_center_y/2

listrq.close()
#end of RQ read

handscan = True#np.any(p_area_ch<-10, axis = (1,2))*np.any((p_class == 1) + (p_class == 2), axis = 1)

#print(p_area_ch[handscan])
make_rq(data_dir, handscan = handscan)

