import numpy as np
import sys


data_dir = "/home/xaber/Analysis/solid_xenon_tpc/data_taking/auto_daq/"

n_boards = 3
n_sipms = [16,8,8]
n_all_ch = int(np.sum(n_sipms))
wsize = 3000+8 # 8 = size of header
load_dtype = "int16"



def getBoard(i):
    if i < 16: return 0
    elif i < 24 and i > 15: return 1
    elif i < 32 and i > 23: return 2
    else: return -1


def getBase(bd):

    baselines = []
    for ch in range(n_sipms[bd]):
        
        # Load data
        ch_data = np.fromfile(data_dir + "waveforms_"+str(bd)+"_"+str(ch)+".dat", dtype=load_dtype)
            
        # Reshape 
        n_events = int(ch_data.size/wsize)
        ch_data = np.reshape(ch_data, (n_events, wsize))

        # Calculate baseline
        baselines.append( int(np.mean(ch_data[:,8:8+100])) )

    baselines = np.asarray(baselines)

    return baselines





def configFile(ch_num):

    if ch_num > 31: return

    


    # Lines in config file to change
    t_start = 161 #141


    trig_lines_0 = np.arange(t_start,t_start+60+4,4,dtype=int)
    trig_lines_1 = np.arange(t_start+68,t_start+96+4,4,dtype=int)
    trig_lines_2 = np.arange(t_start+104,t_start+132+4,4,dtype=int)
    trig_lines_all = np.concatenate((trig_lines_0,trig_lines_1,trig_lines_2), dtype=int)

    trig_enable_0 = np.arange(t_start-1,t_start-1+60+4,4,dtype=int)
    trig_enable_1 = np.arange(t_start-1+68,t_start-1+96+4,4,dtype=int)
    trig_enable_2 = np.arange(t_start-1+104,t_start-1+132+4,4,dtype=int)
    trig_enable_all = np.concatenate((trig_enable_0,trig_enable_1,trig_enable_2))

    # Copy over previous file
    cf_name = "/home/xaber/Data/config_SPE.txt"
    prev = open(cf_name, "r")
    prev_lines = prev.readlines()
    prev.close()




    this_bd = getBoard(ch_num)

    if this_bd == 0:
        prev_lines[23] = "[BOARD 0] Open USB 0 76540000\n"

        baselines = getBase(this_bd) 
        for i in range(n_sipms[this_bd]):
            prev_lines[trig_lines_0[i]] = "TriggerThreshold "+str(baselines[i]+44)+"\n"
            if i == ch_num:
                prev_lines[trig_enable_0[i]] = "EnableInput 1\n"
            else:
                prev_lines[trig_enable_0[i]] = "EnableInput 0\n"



    elif this_bd == 1:
        prev_lines[23] = "[BOARD 0] Open USB 0 32100000\n"

        baselines = getBase(this_bd) 
        for i in range(n_sipms[this_bd]):
            prev_lines[trig_lines_0[i]] = "TriggerThreshold "+str(baselines[i]+44)+"\n"
            if i == ch_num - 16:
                prev_lines[trig_enable_0[i]] = "EnableInput 1\n"
            else:
                prev_lines[trig_enable_0[i]] = "EnableInput 0\n"



    elif this_bd == 2:
        prev_lines[23] = "[BOARD 0] Open USB 0 42100000\n"

        baselines = getBase(this_bd) 
        for i in range(n_sipms[this_bd]):
            prev_lines[trig_lines_0[i]] = "TriggerThreshold "+str(baselines[i]+44)+"\n"
            if i == ch_num - 24:
                prev_lines[trig_enable_0[i]] = "EnableInput 1\n"
            else:
                prev_lines[trig_enable_0[i]] = "EnableInput 0\n"
        



    with open(cf_name, "w") as f:
        for line in prev_lines:
            f.write(line)



    return





def main():
    args = sys.argv[1:]
    configFile(int(args[0])) 

    

if __name__ == "__main__":
    main()
