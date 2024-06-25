#!/usr/bin/env python3
"""
Script to plot slow control values. Rewritten from older version.
- RMG March 2024
"""

import os 
import sys 
import glob
import datetime
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl


# Set range of times to plot
t_plot_start = 0 # in hours
t_plot_end = -1 # in hours -1 means plotting all data points

download_flag = True  # set if need to download log file
no_of_download = 1 # number of files to download
no_of_files = 1 # number of files to plot


def where_to_save():

    # Add your local directory here to avoid typing it out

    if os.path.isdir("C:/Users/maque") or os.path.isdir("/Users/maque/"):
        if os.name == "posix": 
            local_path = "/Users/maque/OneDrive/work/Sxe/log_plot/"
        if os.name == "nt": 
            local_path = "C:/Users/maque/OneDrive/work/Sxe/log_plot/"
    elif os.path.isdir("C:/Users/rmg"):
        local_path = "C:/Users/rmg/Documents/Analysis/logs/"
    else:
        print("Input the local directory to save log. Leave blank for current directory.")
        local_path = input()
        if local_path == "":
            local_path = os.getcwd() + "/"

    return local_path


def download_data(local_path, download_flag):

    if download_flag == False: return

    server = "xaber@xena.dhcp.lbl.gov"
    command = "ls ~/ttlogs"
    process_g = subprocess.Popen(['ssh', server, command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err_g = process_g.communicate()

    all_info = (out.decode(sys.stdout.encoding))
    info_list = all_info.split("\n")
    #print(all_info)
    for i in range(no_of_download):
        t_file = (info_list[-3-i])

        subprocess.call(["scp", server+":"+"/home/xaber/ttlogs/"+t_file, local_path])
    #    print(t_file+" downloaded")

    return


def get_time_offset(file_name):

    t_str = file_name[-17:-4]
    year = int(t_str[:4])
    month = int(t_str[4:6])
    day = int(t_str[6:8])
    hour = int(t_str[9:11])
    minute = int(t_str[11:13])

    return datetime.datetime(year,month,day,hour,minute)

    
def t_p():

    local_path = where_to_save() # where to local logs are

    download_data(local_path, download_flag) # download data

    # Getting local files to plot
    all_files = glob.glob(local_path+"*.csv")
    all_files.sort()
    list_files = all_files[-1*no_of_files:]
    print("Will plot data from the following files:\n"+"\n".join(list_files))

    # Initialize arrays
    #raw_ymd = np.array([], dtype=str)
    #raw_hms = np.array([], dtype=str)
    elapsed_time_s = np.array([])
    t0 = np.array([])
    t5 = np.array([]) 
    t6 = np.array([]) 
    t7 = np.array([])
    t4 = np.array([]) 
    pressure= np.array([])
    top_power = np.array([])
    bot_power = np.array([])

    # Loop over files 
    for i, data_file in enumerate(list_files):

        if i == 0: t0_offset = get_time_offset(data_file)

        # Load each data file and add to arrays
        df = np.loadtxt(data_file, delimiter=",", dtype=str)
        #raw_ymd = np.append(raw_ymd, df[:,0].astype(str))
        #raw_hms = np.append(raw_hms, df[:,1].astype(str))
        try:
            temp_elapsed_time_s = df[:,2].astype(float)
            if i > 0: 
                temp_elapsed_time_s += (get_time_offset(data_file) - t0_offset).total_seconds()
            elapsed_time_s = np.append(elapsed_time_s, temp_elapsed_time_s)
            t0 = np.append(t0, df[:,3].astype(float))
            t5 = np.append(t5, df[:,4].astype(float))
            t6 = np.append(t6, df[:,5].astype(float))
            t7 = np.append(t7, df[:,6].astype(float))
            t4 = np.append(t4, df[:,7].astype(float))
            pressure = np.append(pressure, df[:,8].astype(float))
            top_power = np.append(top_power, df[:,9].astype(float))
            bot_power = np.append(bot_power, df[:,10].astype(float))
        except:
            temp_elapsed_time_s = df[2].astype(float)
            if i > 0: 
                temp_elapsed_time_s += (get_time_offset(data_file) - t0_offset).total_seconds()
            elapsed_time_s = np.append(elapsed_time_s, temp_elapsed_time_s)
            t0 = np.append(t0, df[3].astype(float))
            t5 = np.append(t5, df[4].astype(float))
            t6 = np.append(t6, df[5].astype(float))
            t7 = np.append(t7, df[6].astype(float))
            t4 = np.append(t4, df[7].astype(float))
            pressure = np.append(pressure, df[8].astype(float))
            top_power = np.append(top_power, df[9].astype(float))
            bot_power = np.append(bot_power, df[10].astype(float))


    elapsed_time_hrs = elapsed_time_s / 3600
    top_power = top_power*0.815+0.056 # correction

    # 1. Plot all four temperatures
    fig1, ax1 = pl.subplots()
    ax2=ax1.twinx()
    ax1.plot(elapsed_time_hrs, t0, label="T0, top ICV", color="blue")
    ax1.plot(elapsed_time_hrs, t5, label="T5, copper rod", color="orange")
    ax1.plot(elapsed_time_hrs, t6, label="T6, bot ICV", color="green")
    ax1.plot(elapsed_time_hrs, t7, label="T7, HV cable", color="red")
    ax2.plot(elapsed_time_hrs, t4, label="T4, Room Temperature", color="purple")
    ax1.grid("major", linewidth=0.5, linestyle="dotted")
    ax1.minorticks_on()
    ax1.legend(loc="upper left")
    ax2.legend(loc="upper right")
    ax2.tick_params(axis='y', colors="purple")
    ax1.set_xlabel("Elapsed time [hrs]")
    ax1.set_ylabel("Temperature [deg C]")
    fig1.tight_layout()

    # 2. Plot pressure and heaters
    fig2, ax2 = pl.subplots()
    ax2.plot(elapsed_time_hrs, pressure, color="black")
    ax2.set_ylabel("ICV pressure [bar]")
    ax2.yaxis.label.set_color("black")
    ax2.tick_params(axis="y", labelcolor="black")
    ax2a = ax2.twinx()
    ax2a.plot(elapsed_time_hrs, top_power, label="Top heater", color="darkgreen")
    ax2a.plot(elapsed_time_hrs, bot_power, label="Bot heater", color="darkgreen", linestyle="--")
    ax2a.set_ylabel("Heater power [W]")
    ax2a.yaxis.label.set_color("darkgreen")
    ax2a.tick_params(axis="y", labelcolor="darkgreen")
    ax2.minorticks_on()
    ax2a.minorticks_on()
    #ax2.grid("major", linewidth=0.5, linestyle="dotted")
    #ax2a.grid("major", linewidth=0.5, linestyle="dotted")
    ax2a.legend()
    ax2.set_xlabel("Elapsed time [hrs]")
    fig2.tight_layout()

    # 3 Plot pressure and heaters
    fig3, ax3 = pl.subplots()
    ax3.plot(elapsed_time_hrs, t6, color="black")
    ax3.set_ylabel("ICV bot temperature [deg C]")
    ax3.yaxis.label.set_color("black")
    ax3.tick_params(axis="y", labelcolor="black")
    ax3a = ax3.twinx()
    ax3a.plot(elapsed_time_hrs, top_power, label="Top heater", color="darkgreen")
    ax3a.plot(elapsed_time_hrs, bot_power, label="Bot heater", color="darkgreen", linestyle="--")
    ax3a.set_ylabel("Heater power [W]")
    ax3a.yaxis.label.set_color("darkgreen")
    ax3a.tick_params(axis="y", labelcolor="darkgreen")
    ax3.minorticks_on()
    ax3a.minorticks_on()
    #ax3.grid("major", linewidth=0.5, linestyle="dotted")
    #ax3a.grid("major", linewidth=0.5, linestyle="dotted")
    ax3a.legend()
    ax3.set_xlabel("Elapsed time [hrs]")
    fig3.tight_layout()

    pl.show()


    return


def main():

    t_p()

    return


if __name__ == "__main__":
    main()
