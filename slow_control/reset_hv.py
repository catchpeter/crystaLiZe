#! /bin/user/env ipython
import serial
import time
import os
import sys

print('If the HV power supply had a current trip, this will reset the voltage limit to zero and power it on, after running this script, it will be ready for the ramping script to run.\n')

ps_selection = raw_input('Which HV power supply do you need to reset? 1 for cathode or 2 for the gate.')
if ps_selection == '1':
	device = '/dev/ttyUSB0'
else:
	device = '/dev/ttyUSB2'

ser = serial.Serial(device, 9600,parity='N',timeout=2)
ser.write('TCLR\n')
ser.write('VSET?\n')
vset = ser.readline()

print(vset)

time.sleep(1)
ser.write('VSET 0\n')
ser.write('HVON\n')

ser.close()


