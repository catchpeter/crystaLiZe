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
    print("Please input where you want to save log: ")
    local_path = input()



#set the range to plot
t_plot_start = 0 #in hours
t_plot_end = -1# in hours -1 means plotting all data points
download_flag = True  #set if need to download log file
no_of_download = 1
no_of_files = 1   #number of files to plot


#download the latest log file
if download_flag:
    server = 'xaber@128.3.183.210'

    command = "ls ~/ttlogs"


    process_g = subprocess.Popen(['ssh', server, command], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err_g = process_g.communicate()

    all_info = (out.decode(sys.stdout.encoding))
    info_list = all_info.split("\n")

    for i in range(no_of_download):
        t_file = (info_list[-4-i])

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
temps_smooth = savgol_filter(temps, 51, 3)
t6_smooth = savgol_filter(t6, 81, 3)
pressure_smooth = savgol_filter(pressure, 81, 3)
t5_smooth = savgol_filter(t5, 51, 3)
t7_smooth = savgol_filter(t7, 51, 3)

save_plot = True

t_fit_start = np.where(times>t_plot_start)[0][0]
t_fit_end = -1 if t_plot_end<0 else np.where(times>t_plot_end)[0][0]
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