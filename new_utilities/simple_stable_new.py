#!/usr/bin/env python3
# 2023-12-06 pfs/rmg 
# 2023-12-08 pfs add new cases since pressure is rising beyond where we want to keep it
# 2023-12-10 pfs further tuning on acceptable pressure range
# 2024-03-06 pfs/rmg updating for new MeasurementComputing ADC and new slowcontrol box (xena)

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
	Wt = float(lines[0]) # latest read of top heater power (Watts)
	Wb = float(lines[1]) # latest read of bot heater power (Watts)

#for i in range(1,2):
while True:
	try:
		with open(filename1) as file:
			csvFile = csv.reader(file)
			for lines in csvFile:
				PB = float(lines[7]) # latest read of pressure (Bar)

		print ("pressure now %1.3f and heater powers: %1.3f W (top) %1.3f W (bot)" % (PB,Wt,Wb) )

		# 0.818 Bar triple point	
		# set point p
		spp = 0.84 # crystal triple phase setting
		spp = 1.17 # liquid/vapor (note: approx stable with 2.15,2.15)
		if ( (PB<=(spp+0.01)) and (PB>=(spp-0.01)) ): # do nothing
			print("P=%1.3f Bar - in range!" % PB)
			wriiit(Wt,Wb)							

		if (PB>(spp+0.01)):
			print("P=%1.3f Bar - above range!" % PB)
			wriiit(Wt,Wb*0.9)							
			if (PB>(spp+0.02)):
				wriiit(Wt,Wb*0.8)
				if (PB>(spp+0.03)):
					wriiit(Wt,Wb*0.7)
					if (PB>(spp+0.04)):
						wriiit(Wt,Wb*0.6)
						if (PB>(spp+0.05)):
							wriiit(Wt,Wb*0.5)

		if (PB<(spp-0.01)):
			print("P=%1.3f Bar - below range!" % PB)
			wriiit(Wt,Wb*1.1)
			if (PB<(spp-0.02)):
				wriiit(Wt,Wb*1.2)
				if (PB<(spp-0.03)):
					wriiit(Wt,Wb*1.3)
					if (PB<(spp-0.04)):
						wriiit(Wt,Wb*1.4)
						if (PB<(spp-0.05)):
							wriiit(Wt,Wb*1.5)
		time.sleep(5)
		
	except KeyboardInterrupt: # does not work...
		stored_exception=sys.exc_info()
	







