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

#delay start in seconds. 
#time.sleep(1200)

#Set up freezing parameters, here we assume the bottom heater is off, the top heater is on and its power is changed every once a while.
quick_drop = np.linspace(3.20, 2.7, 29)
slow_drop = np.linspace(2.7, 1.85, 61)
top_power_list = np.append(quick_drop, slow_drop)
initial_top_power = top_power_list[0]
initial_bot_power = 0.
# top_power_limit = 1.3  #the top power will not go lower than this. 
# change_step = 0.08
interval = 3600  # in seconds
bottom_heater_stable_power = 0.7
top_heater_stable_power = 1.05
frozen_pressure = 0.72

#when approching fully frozen, determined by slow_pressure, increase the power so the last part is slower. 
# slow_top_power = 1.5
# slow_bot_power = 0.
# slow_pressure = 0.82  

bot_heater_stable_v = np.sqrt(25.*bottom_heater_stable_power)
top_heater_stable_v = np.sqrt(25.*top_heater_stable_power)

fstr = time.strftime("%Y%m%dT%H%M",time.localtime())

i=0
j=0
start_sec = time.time()
freeze_refer = start_sec # as a starting referrence to count when to change heater power
continue_freeze = True   # should be true to start freezing
#change_to_slow_freezing = True

TC=-1

#initialize the heaters
v_top_i = np.sqrt(25.*initial_top_power)
v_bot_i = np.sqrt(25.*initial_bot_power)
(v_set_i, i_set_i) = setV(v_top_i)
(v_set, i_set) = setV(v_bot_i, 2)
print("{}: Start freezing, the top heater power is now set as {:.2f} Watts, {:.2f} V\n".format(time.strftime("%Y%m%dT%H%M",time.localtime()), initial_top_power, float(v_set_i)/1.23))
freeze_log = open(("/home/xaber/freezelog/%s.txt" % fstr),"a+")
freeze_log.write("Start freezing, the top heater power is now set as {:.2f} Watts, {:.2f} V\n\n".format(initial_top_power, float(v_set_i)/1.23))
freeze_log.close()

flag_check_pressure = True   # should be true to start freezing

while i<1e6: # reads every ~15 sec, so 1e5 => at least 4 days
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
		# also right it to a tmp file, it is used to read out remotely
		with open("/home/xaber/ttlogs/current.csv", "w") as f:
			f.write("%s,%s,%f,%3.3f,%1.3f,%1.3f,%3.3f,%3.3f,%3.3f,%1.3f\n" % (time.strftime("%Y%m%d",time.localtime()),time.strftime("%H%M%S",time.localtime()),time.time()-start_sec,TC,float(v_return)*float(i_return),PB,TC5,TC6,TC7, float(v_bot_return)*float(i_bot_return)))

	i+=1		

	#Check the pressure, if below the  triple point, stop freezing, change heaters power to stable values
	if PB<frozen_pressure and flag_check_pressure:
		flag_check_pressure = False
		continue_freeze = False
		(v1_set, i1_set) = setV(top_heater_stable_v)
		(v2_set, i2_set) = setV(bot_heater_stable_v, 2)
		print("{}: It is fully frozen, top heater is set as {:.2f} V, bottom heater to {} V\n".format(time.strftime("%Y%m%dT%H%M",time.localtime()), top_heater_stable_v, bot_heater_stable_v))
		freeze_log = open(("/home/xaber/freezelog/%s.txt" % fstr),"a+")
		freeze_log.write("{}: It is almost fully frozen, top heater is set as {:.2f} V, bottom heater to {} V\n".format(time.strftime("%Y%m%dT%H%M",time.localtime()), top_heater_stable_v, bot_heater_stable_v))
		freeze_log.close()

	# if PB<slow_pressure and change_to_slow_freezing:
	# 	change_to_slow_freezing = False
	# 	continue_freeze = False
	# 	slow_top_v = np.sqrt(25.*slow_top_power)
	# 	slow_bot_v = np.sqrt(25.*slow_bot_power)
	# 	(v1_set, i1_set) = setV(slow_top_v)
	# 	(v2_set, i2_set) = setV(slow_bot_v, 2)
	# 	print("{}: It is close to being fully frozen, top heater is set as {:.2f} V, bottom heater to {} V\n".format(time.strftime("%Y%m%dT%H%M",time.localtime()), slow_top_v, slow_bot_v))
	# 	freeze_log = open(("/home/xaber/freezelog/%s.txt" % fstr),"a+")
	# 	freeze_log.write("{}: It is close to being fully frozen, top heater is set as {:.2f} V, bottom heater to {} V\n".format(time.strftime("%Y%m%dT%H%M",time.localtime()), slow_top_v, slow_bot_v))
	# 	freeze_log.close()

	#Reduce the top heater power once a while
	if (time.time()-freeze_refer) >interval and (j<(np.size(top_power_list)-1)) and continue_freeze:
		j = j+1
		freeze_refer = time.time()
		new_top_power = top_power_list[j]
		new_heater_voltage = np.sqrt(new_top_power*25.) #heater resistance is 25 ohms.
		(v_set, i_set) = setV(new_heater_voltage)
		print("{}: {:.2f} hours has past, the top heater power is now set as {:.2f} Watts, {:.2f} V\n".format(time.strftime("%Y%m%dT%H%M",time.localtime()), (freeze_refer-start_sec)/3600., new_top_power, float(v_set)/1.23))
		freeze_log = open(("/home/xaber/freezelog/%s.txt" % fstr),"a+")
		freeze_log.write("%s,%s,%f,%3.3f,%1.3f,%1.3f,%3.3f,%3.3f,%3.3f, %1.3f\n" % (time.strftime("%Y%m%d",time.localtime()),time.strftime("%H%M%S",time.localtime()),time.time()-start_sec,TC,float(v_return)*float(i_return),PB,TC5,TC6,TC7, float(v_bot_return)*float(i_bot_return)))
		freeze_log.write("{}: {:.2f} hours has past, the top heater power is now set as {:.2f} Watts, {:.2f} V\n\n".format(time.strftime("%Y%m%dT%H%M",time.localtime()), (freeze_refer-start_sec)/3600., new_top_power, float(v_set)/1.23))
		freeze_log.close()

	

	







