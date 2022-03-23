import numpy as np
import sys

#data_dir = "/home/xaber/Data/20220208/darkCounts_cold_52V/"
#save_dir = "/home/xaber/Analysis/testing/"
#save_file_name = "baselines.txt"
def configFile(ch_num):

    if ch_num > 31: return

    data_dir = "/home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq"


    n_boards = 3
    n_sipms = [16,8,8]
    n_all_ch = int(np.sum(n_sipms))
    wsize = 3000+8 # 8 = size of header
    load_dtype = "int16"


    # Lines in config file to change
    trig_lines_0 = np.arange(141,201+4,4,dtype=int)
    trig_lines_1 = np.arange(209,237+4,4,dtype=int)
    trig_lines_2 = np.arange(245,273+4,4,dtype=int)
    trig_lines_all = np.concatenate((trig_lines_0,trig_lines_1,trig_lines_2), dtype=int)

    trig_enable_0 = np.arange(140,200+4,4,dtype=int)
    trig_enable_1 = np.arange(208,236+4,4,dtype=int)
    trig_enable_2 = np.arange(244,272+4,4,dtype=int)
    trig_enable_all = np.concatenate((trig_enable_0,trig_enable_1,trig_enable_2))

    # Copy over previous file
    cf_name = "/home/xaber/Data/config_SPE.txt"
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
            prev_lines[trig_lines_all[i]] = "TriggerThreshold "+str(baseline+20)+"\n"

            if i==ch_num:
                prev_lines[trig_enable_all[i]] = "EnableInput 1\n"
            else:
                prev_lines[trig_enable_all[i]] = "EnableInput 0\n"

            i += 1
            

    with open(cf_name, "w") as f:
        for line in prev_lines:
            f.write(line)



    return

def main():
    args = sys.argv[1:]
    configFile(int(args[0])) 

    

if __name__ == "__main__":
    main()
