import numpy as np
import os
import time
import datetime

from read_pressure import read_pressure
from read_temp import read_temp
from read_hv import read_cathode, read_gate
from calcBaseSPE import configFile


"""Takes data automatically with wavedumbMB
Input the settings you want below
"""


# High-level location to save data. The full directory is created later
data_dir_high = "/home/xaber/Data/"
#data_dir_high = "/media/xaber/gpeter/data/"


dynamic_range = 1 # 0 = 2Vpp, 1 = 0.5Vpp
event_window_us = 2 #15 # us
pre_trigger = 0.5 # Percentage of event window


trigger_threshold_mV = "NA" # Per channel in mV
run_time_s = 60 #120 # sec, per channel

# Run conditions you need to input
anode_v = 500 # V
sipm_bias = 54 # V
extra = "SPE" # any other info you want to include in dir

# Run conditions that are automatically read
cathode_v = read_cathode() # V
gate_v = read_gate() # V
icv_pressure = read_pressure() # bar
icv_bot_temperature = read_temp() # deg C


# Other globals
if dynamic_range == 0:
    #trigger_threshold_ADCC = int(trigger_threshold_mV/(2000.0/16384.0) )
    dr = "2"
elif dynamic_range == 1:
    #trigger_threshold_ADCC = int(trigger_threshold_mV/(500.0/16384.0) )
    dr = "0.5"
else:
    # Default to 2Vpp
    dynamic_range = 0
    #trigger_threshold_ADCC = int(trigger_threshold_mV/(2000.0/16384.0) )
    dr = "2"


""" Do not edit below here if just taking data
"""

baseline_data_dir = "/home/xaber/Data/temp/"


def SPE_takeBaseData():

    """Takes short data set for baseline calibration
    """

    print("===========\nStarted baseline data collection")

    os.chdir(baseline_data_dir)
    os.system("wavedumpMB /home/xaber/Data/oneBaseConfig.txt 5 1 "+str(int(dynamic_range))+" >/dev/null 2>&1")

    print("===========\nFinished baseline data collection")
    time.sleep(1)
    
    return



def SPE_makeDataDir():

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
    data_dir_daq = dr+"DR_"+str(trigger_threshold_mV)+"mVtrig_"+str(event_window_us)+"us_"
    data_dir_v = str(cathode_v)+"C_"+str(gate_v)+"G_"+str(anode_v)+"A_"+str(sipm_bias)+"SiPM_"
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


def SPE_takeData(data_dir):

    """Take data
    """

    print("===========\nStarted data collection")

    os.chdir(data_dir)

    dataCommand = "wavedumpMB /home/xaber/Data/config_SPE.txt "+str(run_time_s)+" 0 "+str(dynamic_range) 
    ret = os.system(dataCommand)
    if ret != 0:
        print("===========\nError in running wavedumb")
        return


    print("===========\nFinished data collection")

    return


def main():


    n_boards = 3
    n_sipms = [16,8,8]
    n_all_ch = int(np.sum(n_sipms))
    N = 0

    data_dir = SPE_makeDataDir()
    print(data_dir)

    for bd in range(n_boards):

        # Change directory for each board
        bdDir = data_dir+"/b"+str(bd)
        bdmkdirCommand = "mkdir "+bdDir
        os.system(bdmkdirCommand)
        os.chdir(bdDir)

        for ch in range(n_sipms[bd]):
            
            # Take baseline data
            SPE_takeBaseData()

            # Change config file
            configFile(N)

            # Take data
            SPE_takeData(bdDir)

            N += 1
        

    # Create transfer flag file
    os.system("touch "+data_dir+"readyToTransfer")

    
    return



if __name__ == "__main__":
    main()