import numpy as np


#data_dir = "/home/xaber/Data/20220208/darkCounts_cold_52V/"
#save_dir = "/home/xaber/Analysis/testing/"
#save_file_name = "baselines.txt"

data_dir = "/home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq/"


n_boards = 3
n_sipms = [16,8,8]
n_all_ch = int(np.sum(n_sipms))
wsize = 3000+8 # 8 = size of header
load_dtype = "int16"


t_start = 161 #141


trig_lines_0 = np.arange(t_start,t_start+60+4,4,dtype=int)
trig_lines_1 = np.arange(t_start+68,t_start+96+4,4,dtype=int)
trig_lines_2 = np.arange(t_start+104,t_start+132+4,4,dtype=int)
trig_lines_all = np.concatenate((trig_lines_0,trig_lines_1,trig_lines_2), dtype=int)




# Copy over previous file
cf_name = "/home/xaber/Data/config_25us_2V.txt"
prev = open(cf_name, "r")
prev_lines = prev.readlines()
prev.close()




i=0
for bd in range(n_boards):
    for ch in range(n_sipms[bd]):
    
        # Load data
        ch_data = np.fromfile(data_dir + "waveforms_"+str(bd)+"_"+str(ch)+".dat", dtype=load_dtype)
        
        # Reshape 
        n_events = int(ch_data.size/wsize)
        ch_data = np.reshape(ch_data, (n_events, wsize))

        # Calculate baseline
        baseline = int(np.mean(ch_data[:,8:8+100]))   
        
        # Change trigger threshold
        prev_lines[trig_lines_all[i]] = "TriggerThreshold "+str(baseline+100)+"\n"

        i += 1
        

with open(cf_name, "w") as f:
    for line in prev_lines:
        f.write(line)

#for line in prev_lines:
#    print(line)

# Initialize baseline file
#open(save_dir+save_file_name, "w")

"""
# Loop over all channels
for bd in range(n_boards):
    for ch in range(n_sipms[bd]):

        # Load data
        ch_data = np.fromfile(data_dir + "waveforms_"+str(bd)+"_"+str(ch)+".dat", dtype=load_dtype)

        # Reshape 
        n_events = int(ch_data.size/wsize)
        ch_data = np.reshape(ch_data, (n_events, wsize))

        # Calculate baseline
        baseline = np.mean(ch_data[:,8:8+100])   #ch_data[8:8+100]
        #print(baseline)

        # Write value to file
        with open(save_dir+save_file_name,"a") as f:
            f.write(str(baseline)+"\n")
"""

print("done")
