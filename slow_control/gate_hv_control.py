#! /bin/user/env ipython
import serial
import time
import os
import sys
import datetime
from sys import argv
from sys import exit


ser = serial.Serial('/dev/ttyUSB2',9600,parity='N',timeout=2)
log = open('gate_ramp_log.txt', 'a')
print("\n\n--- SRS PS355 High Voltage Controller ---")
print("Change voltage at specified rate\n")

# User io
if len(argv) == 1:
	targetV = abs(float(raw_input("\nEnter voltage (V): ")))
	dvdt = abs(float(raw_input("Enter the voltage change in each step (V) ")))
	interval = abs(int(raw_input("Interval between each step (s): ")))	
elif len(argv) == 4:
	try:
		targetV = abs(float(argv[1]))
		dvdt = abs(float(argv[2]))
		interval = abs(float(argv[3]))
	except ValueError:
		print("Parameters given are not numbers, please correct")
		exit()
else:
	print("Please give parameters in the following order: \n target voltage (V), increase step (V), interval (s).\n")
	exit()

log.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+'\n')
log.write("Start to ramp gate voltage with parameters below\n")
log.write("Target voltage ramp to: {}\n".format(targetV))
log.write("Voltage step: {}\n".format(dvdt))
log.write("Time interval between voltage change: {}\n".format(interval))
# dvdt Check	
dvdtMax = 500
ans = "a"
if dvdt > dvdtMax:
	ans = raw_input("Rate is >500 V/s. Do you want to continue? (y/n) ")
	if ans == "y":
		time.sleep(1)
		print("I hope you don't break anything")
	else:
		exit()

# Reading current voltage
ser.write('VOUT?\n')
v0 = ser.readline()
v0 = v0[:-1]
v0 = abs(float(v0))

# Changing voltage
time.sleep(1)
if targetV  > v0:
	ww = 1
else:
	ww = -1
teaTime = abs((targetV-v0))/dvdt
i = 1
while i <= teaTime:
	vi = v0+dvdt*i*ww
	vSet = 'VSET %f\n' %-vi
	ser.write(vSet)
	
	time.sleep(1)  #wait 1 seconds for the new Vset is applied. 

	#check trip status, if tripped, break out of the loop
	ser.write('*STB? 2\n')
	trip_status = ser.readline()
	if trip_status == '1\n':
		try:
			log.write("The gate is tripped, the trip voltage is: {}.\n\n".format(v_mon))		
			log.close()
			sys.exit('The gate is tripped! The tripped voltage is:{}'.format(v_mon))	
		except(NameError):
			pass
	#print the voltage output. 
	ser.write('VOUT?\n')
	v_mon = ser.readline()
	print(v_mon)

	time.sleep(interval-1)
	i = i+1

# Final confirmation
vSet = 'VSET %f\n' %(-targetV)
ser.write(vSet)
time.sleep(1)
ser.write('VOUT?\n')
vf = ser.readline()
print("Now set at",vf[:-1])
time.sleep(1)
log.write("The gate did not trip, its voltage is: {}\n\n".format(targetV))
#cont = str(raw_input("Set voltage again? (y/n) "))


ser.close()	
log.close()
