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
no_of_download = 2 # number of files to download
no_of_files = 2 # number of files to plot


def where_to_save():

    # Add your local directory here to avoid typing it out

    if os.path.isdir("C:/Users/maque") or os.path.isdir("/Users/maque/"):
        if os.name == "posix": 
            local_path = "/Users/maque/OneDrive/work/Sxe/log_plot/"
        if os.name == "nt": 
            local_path = "C:/Users/maque/OneDrive/work/Sxe/log_plot/"
    elif os.path.isdir("C:/Users/ryanm"):
        local_path = "C:/Users/ryanm/Documents/Research/logs/"
    elif os.path.isdir("C:/Users/rmg"):
        local_path = "C:/Users/rmg/Documents/Analysis/logs/"
    else:
        print("Please input the local directory to save log: ")
        local_path = input()

    return local_path


def download_data(local_path, download_flag):

    if download_flag == False: return

    server = "xaber@xena.dhcp.lbl.gov"
    command = "ls ~/ttlogs"
    process_g = subprocess.Popen(['ssh', server, command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err_g = process_g.communicate()

    all_info = (out.decode(sys.stdout.encoding))
    info_list = all_info.split("\n")
    print(all_info)
    for i in range(no_of_download):
        t_file = (info_list[-3-i])

        subprocess.call(["scp", server+":"+"/home/xaber/ttlogs/"+t_file, local_path])
        print(t_file+" downloaded")

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
    t4 = np.array([])
    t5 = np.array([]) 
    t6 = np.array([]) 
    t7 = np.array([]) 
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
            t4 = np.append(t4, df[:,3].astype(float))
            t5 = np.append(t5, df[:,4].astype(float))
            t6 = np.append(t6, df[:,5].astype(float))
            t7 = np.append(t7, df[:,6].astype(float))
            pressure = np.append(pressure, df[:,7].astype(float))
            top_power = np.append(top_power, df[:,8].astype(float))
            bot_power = np.append(bot_power, df[:,9].astype(float))
        except:
            temp_elapsed_time_s = df[2].astype(float)
            if i > 0: 
                temp_elapsed_time_s += (get_time_offset(data_file) - t0_offset).total_seconds()
            elapsed_time_s = np.append(elapsed_time_s, temp_elapsed_time_s)
            t4 = np.append(t4, df[3].astype(float))
            t5 = np.append(t5, df[4].astype(float))
            t6 = np.append(t6, df[5].astype(float))
            t7 = np.append(t7, df[6].astype(float))
            pressure = np.append(pressure, df[7].astype(float))
            top_power = np.append(top_power, df[8].astype(float))
            bot_power = np.append(bot_power, df[9].astype(float))


    elapsed_time_hrs = elapsed_time_s / 3600
    top_power = top_power*0.815+0.056 # correction

    # 1. Plot all four temperatures
    fig1, ax1 = pl.subplots()
    ax1.plot(elapsed_time_hrs, t4, label="T4, top flange", color="blue")
    ax1.plot(elapsed_time_hrs, t5, label="T5, copper rod", color="orange")
    ax1.plot(elapsed_time_hrs, t6, label="T6, bot ICV", color="green")
    ax1.plot(elapsed_time_hrs, t7, label="T7, HV cable", color="red")
    ax1.grid("major", linewidth=0.5, linestyle="dotted")
    ax1.minorticks_on()
    ax1.legend()
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
    ax2a.plot(elapsed_time_hrs, bot_power, label="Top heater", color="darkgreen", linestyle="--")
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
    ax3a.plot(elapsed_time_hrs, bot_power, label="Top heater", color="darkgreen", linestyle="--")
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


# Old version if ya need it
"""
#!/usr/bin/env python
# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import pandas
from sys import argv
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from glob import glob
from datetime import datetime
import subprocess
import sys
import os

if os.path.isdir("C:/Users/maque") or os.path.isdir("/Users/maque/"):
    if os.name == "posix": local_path = "/Users/maque/OneDrive/work/Sxe/log_plot/"
    if os.name == "nt": local_path = "C:/Users/maque/OneDrive/work/Sxe/log_plot/"
elif os.path.isdir("C:/Users/ryanm"):
    local_path = "C:/Users/ryanm/Documents/Research/logs/"
else:
    local_path = "C:/Users/rmg/Documents/Analysis/logs/"
    #print("Please input where you want to save log: ")
    #local_path = input()


 
#set the range to plot
t_plot_start = 0 #in hours
t_plot_end = -1# in hours -1 means plotting all data points
download_flag = True  #set if need to download log file
no_of_download = 1
no_of_files = 1   #number of files to plot


#download the latest log file
if download_flag:
    server = 'xaber@xena.dhcp.lbl.gov'

    command = "ls ~/ttlogs"


    process_g = subprocess.Popen(['ssh', server, command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err_g = process_g.communicate()

    all_info = (out.decode(sys.stdout.encoding))
    info_list = all_info.split("\n")
    print(all_info)
    for i in range(no_of_download):
        t_file = (info_list[-3-i])

        subprocess.call(["scp", server+":"+"/home/xaber/ttlogs/"+t_file, local_path])
        print(t_file+" downloaded")

all_files = glob(local_path+"*.csv")
all_files.sort()
list_files = all_files[-1*no_of_files:]
#list_files = ["C:/Users/maque/OneDrive/work/Sxe/log_plot/20210729T1238.csv"]
print("Will plot data from the following files:\n"+"\n".join(list_files))

pressure= np.array([])
temps= np.array([])
t5 = np.array([]) 
t6 = np.array([]) 
t7 = np.array([]) 
top_power = np.array([])
bot_power = np.array([])
#temps=t5
raw_times= np.array([])
dates= np.array([])
clock = np.array([])
times = np.array([])

for datafile in list_files:
    df = pandas.read_csv(datafile,names=["date","clock_time","time","temp","top_power","pressure","t5","t6","t7", "bot_power"])


    #df.as_matrix(columns=df.columns)
    #df.loc[[:],['temp']]
    pressure=np.append(pressure, df.get('pressure').values)
    temps=np.append(temps, df.get('temp').values)
    t5 = np.append(t5, df.get('t5').values)
    t6 = np.append(t6, df.get('t6').values)
    t7 = np.append(t7, df.get('t7').values)
    top_power = np.append(top_power, df.get('top_power').values)
    bot_power = np.append(bot_power, df.get('bot_power').values)
    #temps=t5
    clock = df.get('clock_time').values
    for mo in clock: 
        hour = np.floor(mo/10000)
        mins = np.floor((mo-10000*hour)/100)
        seconds = mo-10000*hour-100*mins
        raw_times = np.append(raw_times, hour+mins/60.+seconds/3600.)
    date_format = "%Y%m%d"
    days = []
    ref_date = datetime.strptime("20201024", date_format)
    for i in df.get("date").values:
        a = datetime.strptime(str(i), date_format)
        days.append((a-ref_date).days)
    dates=np.append(dates, days)

t6 = (t6+1.54)/0.97
top_power = top_power*0.815+0.056  #correct the wrong voltgage reading
times = raw_times +24.*dates
times = times - times[0]
temps_smooth = temps #savgol_filter(temps, 51, 3)
t6_smooth = t6 #savgol_filter(t6, 81, 3)
pressure_smooth = pressure #savgol_filter(pressure, 81, 3)
t5_smooth = t5 #savgol_filter(t5, 51, 3)
t7_smooth = t7 #savgol_filter(t7, 51, 3)

save_plot = True

t_fit_start = 0 #np.where(times>t_plot_start)[0][0]
t_fit_end = -1  #if t_plot_end<0 else np.where(times>t_plot_end)[0][0]
print("The range to plot is from {} to {} in hours.".format(t_fit_start, t_fit_end))
date_1 = df.get("date").values[t_fit_end]
clock_1 = df.get("clock_time").values[t_fit_end]

time_1 = datetime.strptime(str(date_1)+str(clock_1), "%Y%m%d%H%M%S").strftime("%Y-%m-%d, %H:%M:%S")
print("The most recent data is taken at "+time_1)

fig, ax1 = plt.subplots()
color1 = 'navy'
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('ICV bottom temperature (C)', color=color1)


ax1.plot(times[t_fit_start:t_fit_end],t6[t_fit_start:t_fit_end],label="ICV bottom", color='orange')

ax1.plot(times[t_fit_start:t_fit_end],t6_smooth[t_fit_start:t_fit_end], label = "ICV bottom smooth", color=color1)
ax1.tick_params(axis='y', labelcolor=color1)
plt.minorticks_on()

#print(fit_function(70))

color2 = 'darkgreen'
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Top heater power (W)', color=color2)  # we already handled the x-label with ax1
ax2.plot(times[t_fit_start:t_fit_end],top_power[t_fit_start:t_fit_end], color=color2, label = "Top heater")
ax2.plot(times[t_fit_start:t_fit_end],bot_power[t_fit_start:t_fit_end], color=color2, linestyle='dashed', label = "Bot heater")
ax2.tick_params(axis='y', labelcolor=color2)
plt.legend(loc="center right", bbox_to_anchor = (1.35, 0.5))
plt.minorticks_on()

#fig.tight_layout()  # otherwise the right y-label is slightly clipped
#if save_plot: plt.savefig(data_dir+'pressure.png', dpi = 150)
plt.suptitle("ICV bottom temperature")
fig.show()

fig, ax1 = plt.subplots()
color1 = 'navy'
ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('Pressure (bar)', color=color1)


ax1.plot(times[t_fit_start:t_fit_end],pressure[t_fit_start:t_fit_end], color='orange')

ax1.plot(times[t_fit_start:t_fit_end],pressure_smooth[t_fit_start:t_fit_end], color=color1)
ax1.tick_params(axis='y', labelcolor=color1)
plt.minorticks_on()

#print(fit_function(70))

color2 = 'darkgreen'
ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

ax2.set_ylabel('Top heater power (W)', color=color2)  # we already handled the x-label with ax1
ax2.plot(times[t_fit_start:t_fit_end],top_power[t_fit_start:t_fit_end], color=color2, label = "Top heater")
ax2.plot(times[t_fit_start:t_fit_end],bot_power[t_fit_start:t_fit_end], color=color2, linestyle='dashed', label = "Bot heater")
ax2.tick_params(axis='y', labelcolor=color2)
plt.legend(loc="center right", bbox_to_anchor = (1.35,0.5))
plt.minorticks_on()
fig.set_size_inches(8, 4.5)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
#if save_plot: plt.savefig(data_dir+'pressure.png', dpi = 150)
plt.suptitle("Pressure", y = 1.05, x = 0.42, fontsize = 18)
fig.show()


#lin_params = scipy.optimize.curve_fit(lin_func,times,temps,p0=[1,20])[0]
#print(lin_params)
#print(-lin_params[1]/lin_params[0]) 

fig, ax = plt.subplots()
#plt.plot(times[t_fit_start:t_fit_end],temps[t_fit_start:t_fit_end],label="Top Flange")
ax.plot(times[t_fit_start:t_fit_end],temps_smooth[t_fit_start:t_fit_end],label="Top Flange smooth")

#ax.plot(times[t_fit_start:t_fit_end],t5[t_fit_start:t_fit_end],label="copper rod")
ax.plot(times[t_fit_start:t_fit_end],t5_smooth[t_fit_start:t_fit_end],label="copper rod smooth")

#ax.plot(times[t_fit_start:t_fit_end],t6[t_fit_start:t_fit_end],label="ICV bottom")
ax.plot(times[t_fit_start:t_fit_end],t6_smooth[t_fit_start:t_fit_end], label = "ICV bottom smooth")

#ax.plot(times[t_fit_start:t_fit_end],t7[t_fit_start:t_fit_end],label="Xe space")
ax.plot(times[t_fit_start:t_fit_end],t7_smooth[t_fit_start:t_fit_end],label="Xe space smooth")

#ax.plot(times,exp_func(times,*exp_params),label="fit")
times_extrap=times*2
#ax.plot(times_extrap,exp_func(times_extrap,*exp_params),label="fit")
#ax.plot(times_extrap,lin_func(times_extrap,*lin_params),label="fit")
ax.legend()
ax.set_xlabel("Time (hr)")
ax.set_ylabel("RTD Temp (C)")
ax.grid()
#ax.yscale('log')
#ax.xlim(0,2)
#ax.ylim(-130,-50)
plt.show()
"""