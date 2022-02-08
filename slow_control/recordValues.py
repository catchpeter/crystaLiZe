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

#print "here we go!"

maxV = 16
scaling = 2 # tried 1, was not enough slope

fstr = time.strftime("%Y%m%dT%H%M",time.localtime())

i=0
start_sec = time.time()
TC=-1

current_top_power = -1
current_bottom_power = -1

while i<1e5: # reads every ~15 sec, so 1e5 => at least 4 days
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
	(v_return,i_return) = getV()
	(v_bot_return, i_bot_return) = getV(2)
	 # default is Ch1
	if 1: # write to file
		# YYYYMMDD,HHMMSS,elapsed_seconds,degreeC,watts,pressureBar
		fid = open(("/home/xaber/ttlogs/%s.csv" % fstr),"a+")
		fid.write("%s,%s,%f,%3.3f,%1.3f,%1.3f,%3.3f,%3.3f,%3.3f,%1.3f\n" % (time.strftime("%Y%m%d",time.localtime()),time.strftime("%H%M%S",time.localtime()),time.time()-start_sec,TC,float(v_return)*float(i_return),PB,TC5,TC6,TC7, float(v_bot_return)*float(i_bot_return)))
		fid.close()
		# also right it to a tmp file, it is used to read out remotely
		with open("/home/xaber/ttlogs/current.csv", "w") as f:
			f.write("%s,%s,%f,%3.3f,%1.3f,%1.3f,%3.3f,%3.3f,%3.3f,%1.3f\n" % (time.strftime("%Y%m%d",time.localtime()),time.strftime("%H%M%S",time.localtime()),time.time()-start_sec,TC,float(v_return)*float(i_return),PB,TC5,TC6,TC7, float(v_bot_return)*float(i_bot_return)))

	i+=1		

	#read the heater setting
	with open("heaters_settting.txt") as f:
	    f.readline()
	    power_line = f.readline()
	    powers = power_line.split(",")
	top_heater_power = float(powers[0])
	bottom_heater_power = float(powers[1])

	if abs(top_heater_power-current_top_power)>0.005:
		print("the current top heater power is {:.2f}W, will change to {:.2f}W".format(current_top_power, top_heater_power))
		top_heater_v = np.sqrt(25.*top_heater_power)
		(v_top_now, i_top_now) = setV(top_heater_v)
		current_top_power = top_heater_power
		
	if abs(bottom_heater_power-current_bottom_power)>0.005:
		print("the current bottome heater power is {:.2f}W, will change to {:.2f}W".format(current_bottom_power, bottom_heater_power))
		bot_heater_v = np.sqrt(25.*bottom_heater_power)
		(v_bottom_now, i_bottom_now) = setV(bot_heater_v, 2)
		current_bottom_power= bottom_heater_power
		







