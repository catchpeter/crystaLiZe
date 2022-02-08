#!/usr/bin/env ipython
# 2017-03-02 pfs - cobbling together from existing scripts written by James and Kate

import sys
import os
import time
import cPickle
import subprocess
import MeasurementComputingUSB as mcUSB
import serial
import numpy as np
import collections
from scipy.signal import savgol_filter


from shared_functions import getT, getV, setV	

def PID(Kp, Ki, Kd, MV_bar=0.6, beta=1, gamma=0):
    # initialize stored data
    eD_prev = 0
    t_prev = -100
    P = 0
    I = 0
    D = 0
    
    # initial control
    MV = MV_bar
    
    while True:
        # yield MV, wait for new t, SP, PV, TR
        data = yield MV
        
        # see if a tracking data is being supplied
        if len(data) < 4:
            t, PV, SP = data
        else:
            t, PV, SP, TR = data
            I = TR - MV_bar - P - D
        
        # PID calculations
        P = Kp*(beta*SP - PV)
        I = I + Ki*(SP - PV)*(t - t_prev)
        eD = gamma*SP - PV
        D = Kd*(eD - eD_prev)/(t - t_prev)
        MV = MV_bar + P + I + D
        print(SP, PV)
        print(MV_bar, P, I, D, t, t_prev)
        # Constrain MV to range 0 to 100 for anti-reset windup
        MV = 0 if MV < 0 else 100 if MV > 100 else MV
        #I = MV - MV_bar - P - D
        
        # update stored data for next iteration
        eD_prev = eD
        t_prev = t

def time_to_second(time_str):
	return int(time_str[:8])*24*3600+int(time_str[9:11])*3600+int(time_str[11:13])*60+int(time_str[13:])

def queue_action(deque_list, value): #append value to the end, and pop the first elememnt
	deque_list.popleft()
	deque_list.append(value)

#main
bottome_t_set_point = -105.9
controller = PID(0.08, 1.11e-5, 0, beta=0)   # create pid control
controller.send(None)                 # initialize
#save the most recent 2000 of data points
no_recent = 100
t_recent = collections.deque(np.zeros(no_recent))
bottom_t_recent = collections.deque(np.zeros(no_recent))

fstr = time.strftime("%Y%m%dT%H%M%S",time.localtime())
t_intitial = time_to_second(fstr)

i=0
start_sec = time.time()
TC=-1
while i<1e5: # reads every ~15 sec, so 1e5 => at least 4 days
	timestr = time.strftime("%Y%m%dT%H%M%S",time.localtime())
	read_mcc_success = 0
	while (read_mcc_success != 1):
		try:
			(TC,PB,TC5,TC6,TC7) = getT()
			read_mcc_success = 1
		except:
			#print ("TC = %3.2f" % TC)
			print "read failed (!!) waiting 15 secs and trying again"
			time.sleep(15)
			print "OK lets try again..."
	print ("%s : T4=%3.3f C, T5=%3.3f C, T6=%3.3f C, T7=%3.3f C, P=%1.3f Bar, iteration=%d" % (timestr,TC,TC5,TC6,TC7,PB,i))

	#save recent data for PID controll
	queue_action(t_recent, time_to_second(timestr) - t_intitial)
	queue_action(bottom_t_recent, TC6)
	bottom_t_recent_smooth = savgol_filter(bottom_t_recent, int(no_recent/4)*2+1, 3)  #window width needs to be odd


	(v_return,i_return) = getV()
	(v_bot_return, i_bot_return) = getV(2)

	#feed to PID
	current_bottom_power = float(v_bot_return)*float(i_bot_return)
	index_used = int(no_recent/2)
	next_bottom_power = controller.send([t_recent[index_used], bottom_t_recent[index_used], bottome_t_set_point]) 
	#next_bottom_power = controller.send([t_recent[index_used], bottom_t_recent[index_used], bottome_t_set_point, current_bottom_power]) 

	print("t_recent is {}, index is {}.".format(t_recent[index_used], index_used))
	print("Bottome heater is set at {:.3f}, PID suggest: {:.3f}".format(current_bottom_power, next_bottom_power))

	if (i>no_recent) and (i%80 == 0):
		next_bottom_power = 0.65 if next_bottom_power>0.65 else 0.55 if next_bottom_power<0.55 else next_bottom_power
		v_bot = np.sqrt(25.*next_bottom_power)
		(v_set, i_set) = setV(v_bot, 2)
		print("change made")


	 # default is Ch1
	if 1: # write to file
		# YYYYMMDD,HHMMSS,elapsed_seconds,degreeC,watts,pressureBar
		fid = open(("/home/xaber/ttlogs/%s.csv" % fstr),"a+")
		info = []
		info.append(time.strftime("%Y%m%d",time.localtime()))
		info.append(time.strftime("%H%M%S",time.localtime()))
		info.append(str(time.time()-start_sec))
		info.append("{:3.3f}".format(TC))
		info.append("{:1.3f}".format(float(v_return)*float(i_return)))
		info.append("{:1.3f}".format(PB))
		info.append("{:3.3f}".format(TC5))
		info.append("{:3.3f}".format(TC6))
		info.append("{:3.3f}".format(TC7))
		info.append("{:1.3f}\n".format(float(v_bot_return)*float(i_bot_return)))
		fid.write(", ".join(info))
		fid.close()
	i+=1		







