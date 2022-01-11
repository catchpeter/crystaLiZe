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

from shared_functions import getT, getV, setV	

#setting the melting parameter
bot_initial_power = 16
top_initial_power = 16
bot_v = np.sqrt(25.*bot_initial_power)
top_v = np.sqrt(25.*top_initial_power)
time_high_power = 1500 #in seconds
bot_power_stable = 1.6
top_power_stable = 1.7

bot_v_stable = np.sqrt(25.*bot_power_stable)
top_v_stable = np.sqrt(25.*top_power_stable)

fstr = time.strftime("%Y%m%dT%H%M",time.localtime())

i=0
start_sec = time.time()

TC=-1
lower_power = True
#initialize the heaters
(v_set_i, i_set_i) = setV(top_v)
(v_set, i_set) = setV(bot_v, 2)
melting_log = open(("/home/xaber/freezelog/%s.txt" % fstr),"a+")
melting_log.write("{}. Start melting, the top heater power is now set as {:.2f} V, the bot heater is set as {:.2f}.\n\n".format(time.strftime("%Y%m%dT%H%M",time.localtime()), top_v, bot_v))
melting_log.close()

while i<1e5: # reads every ~15 sec, so 1e5 => at least 4 days
	#change top heater power
	timestr = time.strftime("%Y%m%dT%H%M",time.localtime())
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
	(v_return,i_return) = getV() # default is Ch1
	(v_bot_return, i_bot_return) = getV(2)
	if 1: # write to file
		# YYYYMMDD,HHMMSS,elapsed_seconds,degreeC,watts,pressureBar
		fid = open(("/home/xaber/ttlogs/%s.csv" % fstr),"a+")
		fid.write("%s,%s,%f,%3.3f,%1.3f,%1.3f,%3.3f,%3.3f,%3.3f, %1.3f\n" % (time.strftime("%Y%m%d",time.localtime()),time.strftime("%H%M%S",time.localtime()),time.time()-start_sec,TC,float(v_return)*float(i_return),PB,TC5,TC6,TC7, float(v_bot_return)*float(i_bot_return)))
		fid.close()
	i+=1		

	#lower the power to stable setting
	if (((time.time()-start_sec)>time_high_power) or PB>1.5) and lower_power:
		(v_top_set, i_top_set) = setV(top_v_stable)
		(v_bot_set, i_bot_set) = setV(bot_v_stable, 2)
		melting_log = open(("/home/xaber/freezelog/%s.txt" % fstr),"a+")
		melting_log.write("{} change the top power to {:.2f}, bot power to {:.2f}\n\n ".format(time.strftime("%Y%m%dT%H%M",time.localtime()), top_v_stable, bot_v_stable))
		melting_log.close()
		lower_power = False
	

	







