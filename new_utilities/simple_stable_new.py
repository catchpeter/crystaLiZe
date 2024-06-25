#!/usr/bin/env python3
# 2023-12-06 pfs/rmg 
# 2023-12-08 pfs add new cases since pressure is rising beyond where we want to keep it
# 2023-12-10 pfs further tuning on acceptable pressure range
# 2024-03-06 pfs/rmg updating for new MeasurementComputing ADC and new slowcontrol box (xena)
# 2024-04-22 sx optimize heater control algorithm

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

scale=2
dp = 0.005
# 0.818 Bar triple point	
# set point p
spp = 0.81 # crystal triple phase setting
# spp = 1.17 # liquid/vapor (note: approx stable with 2.15,2.15)

#for i in range(1,2):
while True:
	try:
		with open(filename1) as file:
			csvFile = csv.reader(file)
			for lines in csvFile:
				PB = float(lines[8]) # latest read of pressure (Bar)

		with open(filename2) as file:
			csvFile = csv.reader(file)
			for lines in csvFile: # ignores the first line
				tmp=lines[0]
			Wtnow = float(lines[0]) # starting read of top heater power (Watts) -- does not get modified
			Wbnow = float(lines[1]) # starting read of bot heater power (Watts) -- does not get modified

		print ("baseline goal p=%1.3f with Wt=%1.3f and Wb=%1.3f" % (spp,Wt,Wb))
		print (".......... latest: p=%1.3f and heater powers: %1.3f W (top) %1.3f W (bot)" % (PB,Wtnow,Wbnow) )

		if ( (PB<=(spp+dp)) and (PB>=(spp-dp)) ): # do nothing
			print("P=%1.3f Bar - in range!" % PB)
			if abs(Wb-Wbnow)>0.01:
				wriiit(Wt,Wb)			
			else:
				print ("same setting as previous step")				
		elif (PB>(spp+5*dp)):
			print("P=%1.3f Bar - above range!" % PB)
			if abs(Wb*0.3-Wbnow)>0.01:
				wriiit(Wt,Wb*0.3)
			else:
				print ("same setting as previous step")				
		elif (PB>(spp+4*dp)):
			print("P=%1.3f Bar - above range!" % PB)
			if abs(Wb*0.4-Wbnow)>0.01:	
				wriiit(Wt,Wb*0.4)
			else:
				print ("same setting as previous step")				
		elif (PB>(spp+3*dp)):
			print("P=%1.3f Bar - above range!" % PB)
			if abs(Wb*0.5-Wbnow)>0.01:	
				wriiit(Wt,Wb*0.5)
			else:
				print ("same setting as previous step")				
		elif (PB>(spp+2*dp)):
			print("P=%1.3f Bar - above range!" % PB)
			if abs(Wb*(1-scale*0.2)-Wbnow)>0.01:	
				wriiit(Wt,Wb*(1-scale*0.2))
			else:
				print ("same setting as previous step")				
		elif (PB>(spp+dp)):
			print("P=%1.3f Bar - above range!" % PB)
			if abs(Wb*(1-scale*0.1)-Wbnow)>0.01:	
				wriiit(Wt,Wb*(1-scale*0.1))							
			else:
				print ("same setting as previous step")				
		elif (PB<(spp-5*dp)):
			print("P=%1.3f Bar - below range!" % PB)
			if abs(Wb*1.8-Wbnow)>0.01:	
				wriiit(Wt*1.5,Wb*1.8)	
			else:
				print ("same setting as previous step")				
		elif (PB<(spp-4*dp)):
			print("P=%1.3f Bar - below range!" % PB)
			if abs(Wb*1.7-Wbnow)>0.01:	
				wriiit(Wt*1.4,Wb*1.7)
			else:
				print ("same setting as previous step")				
		elif (PB<(spp-3*dp)):
			print("P=%1.3f Bar - below range!" % PB)
			if abs(Wb*(1+scale*0.3)-Wbnow)>0.01:	
				wriiit(Wt*1.3,Wb*(1+scale*0.3))
			else:
				print ("same setting as previous step")				
		elif (PB<(spp-2*dp)):
			print("P=%1.3f Bar - below range!" % PB)
			if abs(Wb*(1+scale*0.2)-Wbnow)>0.01:	
				wriiit(Wt*1.2,Wb*(1+scale*0.2))
			else:
				print ("same setting as previous step")				
		elif (PB<(spp-dp)):
			print("P=%1.3f Bar - below range!" % PB)
			if abs(Wb*(1+scale*0.1)-Wbnow)>0.01:	
				wriiit(Wt,Wb*(1+scale*0.1))
			else:
				print ("same setting as previous step")				
		time.sleep(5)
		
	except KeyboardInterrupt: # does not work...
		stored_exception=sys.exc_info()
	







