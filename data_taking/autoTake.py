import numpy as np
import os
import time
import datetime

from read_pressure import read_pressure
from read_temp import read_temp
from read_hv import read_cathode, read_gate


"""Takes data automatically with wavedumbMB
Input the settings you want below
"""


# High-level location to save data. The full directory is created later
data_dir_high = "/home/xaber/Data/"
#data_dir_high = "/media/xaber/gpeter/data/"

# Run settings you need to input
dynamic_range = 1 # 0 = 2Vpp, 1 = 0.5Vpp
event_window_us = 15 # us
pre_trigger = 0.5 # Percentage of event window
trigger_threshold_mV = 6 # Per channel in mV
run_time_s = 10 # sec

# Run conditions you need to input
anode_v = 0 # kV
sipm_bias = 54 # V
extra = "" # any other info you want to include in dir

# Run conditions that are automatically read
cathode_v = read_cathode() # kV
gate_v = read_gate() # kV
icv_pressure = read_pressure() # bar
icv_bot_temperature = read_temp() # deg C



""" Do not edit below here if just taking data
"""


# Other globals
if dynamic_range == 0:
    trigger_threshold_ADCC = int(trigger_threshold_mV/(2000.0/16384.0) )
    dr = "2"
elif dynamic_range == 1:
    trigger_threshold_ADCC = int(trigger_threshold_mV/(500.0/16384.0) )
    dr = "0.5"
else:
    # Default to 2Vpp
    dynamic_range = 0
    trigger_threshold_ADCC = int(trigger_threshold_mV/(2000.0/16384.0) )
    dr = "2"

baseline_data_dir = "/home/xaber/Data/temp/"



def takeBaseData():

    """Takes short data set for baseline calibration
    """

    print("===========\nStarted baseline data collection")

    os.chdir(baseline_data_dir)
    os.system("wavedumpMB /home/xaber/Data/oneBaseConfig.txt 5 1 "+str(int(dynamic_range))+" >/dev/null 2>&1")

    print("===========\nFinished baseline data collection")
    time.sleep(1)
    
    return



def mkConfig():

    """Calculates baseline, edits other config file settings
    """

    print("===========\nStarted writing config file")

    n_boards = 3
    n_sipms = [16,8,8]
    wsize = 3000+8 # 8 = size of header
    load_dtype = "int16"

    # Config file trigger lines
    t_start = 161 #141
    trig_lines_0 = np.arange(t_start,t_start+60+4,4,dtype=int)
    trig_lines_1 = np.arange(t_start+68,t_start+96+4,4,dtype=int)
    trig_lines_2 = np.arange(t_start+104,t_start+132+4,4,dtype=int)
    trig_lines_all = np.concatenate((trig_lines_0,trig_lines_1,trig_lines_2), dtype=int)

    # Copy over previous file
    cf_name = "/home/xaber/Data/config_normal.txt"
    prev = open(cf_name, "r")
    prev_lines = prev.readlines()
    prev.close()

    # Change event window and pre-trigger
    evWinSamp = int(event_window_us*1000/2)
    preTrigSamp = int(evWinSamp*pre_trigger)
    prev_lines[32] = "RecordLength                  "+str(evWinSamp)+"\n"
    prev_lines[33] = "PreTrigger                    "+str(preTrigSamp)+"\n"

    # Loop over channels and change trigger threshold                                     
    i=0
    for bd in range(n_boards):
        for ch in range(n_sipms[bd]):
        
            # Load data
            ch_data = np.fromfile(baseline_data_dir + "waveforms_"+str(bd)+"_"+str(ch)+".dat", dtype=load_dtype)
            
            # Reshape 
            n_events = int(ch_data.size/wsize)
            ch_data = np.reshape(ch_data, (n_events, wsize))

            # Calculate baseline
            baseline = int(np.mean(ch_data[:,8:8+100]))   
            
            # Change trigger threshold
            prev_lines[trig_lines_all[i]] = "TriggerThreshold "+str(baseline+trigger_threshold_ADCC)+"\n"

            i += 1
            

    with open(cf_name, "w") as f:
        for line in prev_lines:
            f.write(line)

    print("===========\nFinished writing config file")


    return



def makeDataDir():

    """Creates data directory

    Returns: 
     data_dir - str of directory if successful, 1 if unsuccessful
    """

    # Get the day and time and format into string
    # There must be an easier way to do this...
    dt_now = datetime.datetime.now()
    year = str(dt_now.year)
    month = str(dt_now.month) if dt_now.month > 9 else "0"+str(dt_now.month)
    day = str(dt_now.day) if dt_now.day > 9 else "0"+str(dt_now.day)
    hour = str(dt_now.hour) if dt_now.hour > 9 else "0"+str(dt_now.hour)
    minute = str(dt_now.minute) if dt_now.minute > 9 else "0"+str(dt_now.minute)
    ym = year+month
    ymd = ym+day 
    hm = hour+minute

    # Format data_dir
    data_dir_dt = "data-"+ym+"/"+ymd+"/"+ymd+"-"+hm+"_"
    data_dir_daq = dr+"DR_"+str(trigger_threshold_mV)+"mVtrig_"+str(cathode_v)
    data_dir_v = "C_"+str(gate_v)+"G_"+str(anode_v)+"A_"+str(sipm_bias)+"SiPM_"
    data_dir_tp = str(icv_pressure)+"bar_"+str(icv_bot_temperature)+"ICVbot_"+extra+"/"
    data_dir = data_dir_high + data_dir_dt + data_dir_daq + data_dir_v + data_dir_tp

    mkdirCommand = "mkdir "+data_dir+" -p"
    ret = os.system(mkdirCommand)

    if ret == 0:
        print("===========\nCreated directory:\n    "+data_dir)
        return data_dir
    else:
        print("===========\nError in creating directory:\n    "+data_dir)
        return 1
        


def takeData(data_dir):

    """Take data
    """

    print("===========\nStarted data collection")

    os.chdir(data_dir)

    dataCommand = "wavedumpMB /home/xaber/Data/config_normal.txt "+str(run_time_s)+" 0 "+str(dynamic_range) 
    ret = os.system(dataCommand)
    if ret != 0:
        print("===========\nError in running wavedumb")
        return


    print("===========\nFinished data collection")

    return



def main():

    # Take baseline data
    takeBaseData()

    # Change config file
    mkConfig()

    # Create data directory
    data_dir = makeDataDir()
    #if data_dir: return

    # Take data
    takeData(data_dir)
    
    return



if __name__ == "__main__":
    main()
