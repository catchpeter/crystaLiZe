from subprocess import run
import subprocess
import numpy as np
import os, sys
import time
import datetime

from read_pressure import read_pressure
from read_temp import read_temp
from read_hv import read_cathode, read_gate
from grid_voltage import grid_voltage


"""Takes data automatically with wavedumbMB
Input the settings you want below
"""
t_delay = 0 #3500, will start data t_delay seconds later.
process_flag = False # will compress and process data if this is true.

while t_delay>0:
    print("Start taking data in {} mins.".format(t_delay/60))
    time.sleep(100)
    t_delay -= 100
    print("Start taking data in {} mins.".format(t_delay/60))
# High-level location to save data. The full directory is created later
data_dir_high = "/home/xaber/Data/"
#data_dir_high = "/media/xaber/gpeter/data/"

# Run settings you need to input
event_window_us = 20 #20 #15 # us
pre_trigger = 0.5 # Percentage of event window
trigger_threshold_mV = 10 # Per channel in mV
run_time_s = 30*60 # sec

# Run conditions you need to input
phase = "liquid" 
anode_v = 500 # V
sipm_bias = 54 # V
source = "Co57"
extra = "test" # any other info you want to include in dir

# Run conditions that are automatically read
cathode_v = read_cathode() # V
gate_v = read_gate() # V
icv_pressure = read_pressure() # bar
icv_bot_temperature = read_temp() # deg C
run_time_min = "{:n}".format(run_time_s/60)

                                        
""" Do not edit below here if just taking data
"""
# TO CHANGE DYNAMIC RANGE YOU NEED TO EDIT WAVEDUMP C CODE AND RECOMPILE
# YOU STILL NEED TO CHANGE THIS VALUE SO THE ANALYSIS WORKS
dynamic_range = 0 # 0 = 2Vpp, 1 = 0.5Vpp

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
            try:
                baseline = int(np.mean(ch_data[:,8:8+100]))
            except:
                print(bd, ch)
                baseline = 99999
            
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
    # Update: wow I actually wrote this shitty code
    dt_now = datetime.datetime.now()
    year = str(dt_now.year)
    month = str(dt_now.month) if dt_now.month > 9 else "0"+str(dt_now.month)
    day = str(dt_now.day) if dt_now.day > 9 else "0"+str(dt_now.day)
    hour = str(dt_now.hour) if dt_now.hour > 9 else "0"+str(dt_now.hour)
    minute = str(dt_now.minute) if dt_now.minute > 9 else "0"+str(dt_now.minute)
    second = str(dt_now.second) if dt_now.second > 9 else "0"+str(dt_now.second)
    ym = year+month
    ymd = ym+day 
    hms = hour+minute+second

    # Format data_dir
    data_dir = data_dir_high + "data-"+ym+"/"+ymd+"/"+ymd+"-"+hms+"/"

    # Old data dir naming
    #data_dir_dt = "data-"+ym+"/"+ymd+"/"+ymd+"-"+hm+"_"
    #data_dir_daq = dr+"DR_"+str(trigger_threshold_mV)+"mVtrig_"+str(event_window_us)+"us_"
    #data_dir_v = str(cathode_v)+"C_"+str(gate_v)+"G_"+str(anode_v)+"A_"+str(sipm_bias)+"SiPM_"
    #data_dir_tp = str(icv_pressure)+"bar_"+str(icv_bot_temperature)+"ICVbot_"+phase+"_source"+source+"_"+extra+"/"
    #data_dir = data_dir_high + data_dir_dt + data_dir_daq + data_dir_v + data_dir_tp

    
    mkdirCommand = "mkdir "+data_dir+" -p"
    ret = os.system(mkdirCommand)

    if ret == 0:
        print("===========\nCreated directory:\n    "+data_dir)
        return data_dir
    else:
        print("===========\nError in creating directory:\n    "+data_dir)
        return 1
        

def writeConditions(data_dir):


    fields = []
    fields.append(f"Date time")
    fields.append(f"Run time min")

    fields.append(f"Dynamic range")
    fields.append(f"Trigger threshold mV")
    fields.append(f"Event window us")

    fields.append(f"Cathode voltage")
    fields.append(f"Gate voltage")
    fields.append(f"Anode voltage")
    fields.append(f"SiPM bias")

    fields.append(f"ICV pressure bar")
    fields.append(f"ICV bottom temp deg C")
    fields.append(f"Phase")
    fields.append(f"Source")
    fields.append(f"Misc")



    cond_list = []
    cond_list.append(f"{data_dir[-16:-1]}")
    cond_list.append(f"{run_time_min}")

    cond_list.append(f"{dr}")
    cond_list.append(f"{trigger_threshold_mV}")
    cond_list.append(f"{event_window_us}")

    cond_list.append(f"{cathode_v}")
    cond_list.append(f"{gate_v}")
    cond_list.append(f"{anode_v}")
    cond_list.append(f"{sipm_bias}")

    cond_list.append(f"{icv_pressure}")
    cond_list.append(f"{icv_bot_temperature}")
    cond_list.append(f"{phase}")
    cond_list.append(f"{source}")
    cond_list.append(f"{extra}")


    rows = [fields, cond_list]
  


    np.savetxt(data_dir+"/conditions.csv", rows, delimiter=",", fmt="%s")






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

"""
def process_data(data_dir):
    #Place the path of this data to path.txt ready to be processed
    process_path = "/home/xaber/Analysis/solid_xenon_tpc/pulse_finding_scripts/"
    sys.path.append(process_path)

    from rq_generate_32ch import make_rq
    from cut_plot import make_plots
    from compression import compression

    os.chdir(process_path)
    compression(data_dir)
    make_rq(data_dir)
    make_plots(data_dir)
"""
    
def main():

    #Ramp up voltage
    # ramp_flag = grid_voltage(gate = 3000, cathode = 3200)
    # time.sleep(300)

    # Take baseline data
    takeBaseData()

    # Change config file
    mkConfig()

    # Create data directory
    data_dir = makeDataDir()
    if data_dir == 1: return

    # Create conditions log
    writeConditions(data_dir)

    # Take data
    takeData(data_dir)

    # Create transfer flag file
    os.system("touch "+data_dir+"readyToTransfer")

    # Ramp down voltage
    #ramp_flag = grid_voltage(gate = 0, cathode = 0)

    # Process the data
    #if process_flag: process_data(data_dir)

    # upload data to google drive
    #os.system("rclone -v copy " + data_dir +" gdrive:crystallize/data/"+ data_dir[data_dir.find("data-"):])



if __name__ == "__main__":
    main()