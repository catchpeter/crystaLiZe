#! /bin/user/env ipython
import serial
import time
import os
import sys

print('Check the status of HV power supply\n')

ps_selection = raw_input('Which HV power supply do you want to check? 1 for cathode or 2 for the gate.')
if ps_selection == '1':
	device = '/dev/ttyUSB0'
else:
	device = '/dev/ttyUSB2'

ser = serial.Serial(device, 9600,parity='N',timeout=2)


time.sleep(1)
ser.write('VOUT?\n') 
vmon = ser.readline()

print(vmon)
ser.close()


