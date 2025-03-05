#!/usr/bin/env python3
# 2023-12-06 pfs/rmg 
# 2023-12-08 pfs add new cases since pressure is rising beyond where we want to keep it
# 2023-12-10 pfs further tuning on acceptable pressure range
# 2024-03-06 pfs/rmg updating for new MeasurementComputing ADC and new slowcontrol box (xena)
# 2025-03-05 pfs update to keep stable at -100 C (not using pressure)

import sys
import os
import time
import subprocess
import serial
import numpy as np
import csv


def wriiit(Wt,Wb):
	f=open(filename2,"w")
	f.write("top, bot\n%1.2f, %1.2f\n" % (Wt,Wb))
	f.close()

filename1='/home/xaber/ttlogs/current.csv'
filename2='/home/xaber/crystaLiZe/new_utilities/heaters_setting.txt'

# read heaters_setting.txt ONCE at the start to get the baseline settings
with open(filename2) as file:
	csvFile = csv.reader(file)
	for lines in csvFile: # ignores the first line
		tmp=lines[0]
	Wt = float(lines[0]) # starting read of top heater power (Watts) -- does not get modified
	Wb = float(lines[1]) # starting read of bot heater power (Watts) -- does not get modified

dp = 0.1
# 0.818 Bar triple point	
# set point p
#spp = 0.822 # crystal triple phase setting
# spp = 1.17 # liquid/vapor (note: approx stable with 2.15,2.15)
spp = -103

#for i in range(1,2):
while True:
	try:
		with open(filename1) as file:
			csvFile = csv.reader(file)
			for lines in csvFile:
				T0 = float(lines[3]) # latest read of pressure (Bar)
# 				print(lines)
# 				PB = float(lines[7]) # latest read of pressure (Bar)

		with open(filename2) as file:
			csvFile = csv.reader(file)
			for lines in csvFile: # ignores the first line
				tmp=lines[0]
			Wtnow = float(lines[0]) # starting read of top heater power (Watts) -- does not get modified
			Wbnow = float(lines[1]) # starting read of bot heater power (Watts) -- does not get modified

		print ("baseline goal p=%1.3f with Wt=%1.3f and Wb=%1.3f" % (spp,Wt,Wb))
		print (".......... latest: p=%1.3f and heater powers: %1.3f W (top) %1.3f W (bot)" % (PB,Wtnow,Wbnow) )

		if ( (T0<=(spp+dp)) and (T0>=(spp-dp)) ): # in range - do nothing
			print("T0=%1.3f C - in range!" % T0)
			wriiit(Wt,Wb)							

		if (T0>(spp+dp)):
			print("T0=%1.3f C - above range!" % T0)
			wriiit(Wt,Wb*0.9)							
			if (T0>(spp+2*dp)):
				wriiit(Wt,Wb*0.8)
				if (T0>(spp+3*dp)):
					wriiit(Wt,Wb*0.7)
					if (T0>(spp+4*dp)):
						wriiit(Wt,Wb*0.6)
						if (T0>(spp+5*dp)):
							wriiit(Wt,Wb*0.5)

		if (T0<(spp-dp)):
			print("T0=%1.3f C - below range!" % T0)
			wriiit(Wt,Wb*1.1)
			if (T0<(spp-2*dp)):
				wriiit(Wt,Wb*1.2)
				if (T0<(spp-3*dp)):
					wriiit(Wt,Wb*1.3)
					if (T0<(spp-4*dp)):
						wriiit(Wt,Wb*1.4)
						if (T0<(spp-5*dp)):
							wriiit(Wt,Wb*1.5)
		time.sleep(5)
		
	except KeyboardInterrupt: # does not work...
		stored_exception=sys.exc_info()
	







