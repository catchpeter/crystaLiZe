#!/usr/bin/env ipython
# 2020-08-16 swk - copied/edited from scripts by Peter, James, Kate

import sys
import os
import time
import cPickle
import subprocess
import MeasurementComputingUSB as mcUSB
import serial
import numpy as np

def getT():
	mcc=mcUSB.MCC_PMD_1208LS()
	time.sleep(1.0) # this is important
	mcc.blink_led()
	oldid=mcc.get_id()
	mcc.set_id( (oldid+1) & 255 )
	id_dict=cPickle.loads(mcc.read_user_memory())
	newid=mcc.get_id()
	id_dict['serial']=newid
	id_dict['last_used_date']=time.asctime()
	id_dict['last_used_time']=time.time()
	mcc.write_user_memory(data=cPickle.dumps(id_dict) )

	ch0 = mcc.analog_input(0, gain=mcc.GAIN1_SE)
	ch1 = mcc.analog_input(1, gain=mcc.GAIN1_SE)
	ch4 = mcc.analog_input(4, gain=mcc.GAIN1_SE)
	ch5 = mcc.analog_input(5, gain=mcc.GAIN1_SE)
	ch6 = mcc.analog_input(6, gain=mcc.GAIN1_SE)
	ch7 = mcc.analog_input(7, gain=mcc.GAIN1_SE)

	# conversions from James' script daqDaemon.py, have not checked
	#tempC = ch0*30-150
	#pressureBar = ch1*1.8
	#massflowSLPM = ch4

	TC = (ch0*30-149.3) # changed scale factor from James values but calibrated only at -102
	TC5 = (ch5*30-149.3)
	TC6 = (ch6*30-149.3)
	TC7 = (ch7*30-149.3) 
	PB = (ch1*1.760) # changed scale factor from James values but calibrated only at 1.60
	#print "temp C: %1.3f" % TC 
	#print "pressure: %1.3f" % pressureBar 
	#print "massflowSLPM: %1.3f" % massflowSLPM

	mcc.close()
	return (TC,PB,TC5,TC6,TC7)
	
def query_serial(ser,query):
	if ser.isOpen():
		ser.write("%s%s%s" % (query,serial.CR,serial.LF))
		output = ser.readline()
	else:
		ser.open()
		ser.write("%s%s%s" % (query,serial.CR,serial.LF))
		output = ser.readline()
	return output
	
def setV(voltage, channel=1):
	# Ch 1 is usually top heater
	# parameters:
	voltage_range = 'HIGH' # or 'LOW'
	ser = serial.Serial('/dev/ttyUSB1',9600,parity='N',bytesize=8,stopbits=2,timeout=2)
	output = query_serial(ser,"*IDN?")
	query_serial(ser,"SYSTem:REMote")
	query_serial(ser,"INSTrument:NSELect %s" % channel)
	query_serial(ser,"VOLTage:RANGe %s" %(voltage_range))
	query_serial(ser,"VOLTage %f" % (voltage))
	query_serial(ser,"OUTPut ON")
	v_return = query_serial(ser,"MEASure:VOLTage?")
	i_return = query_serial(ser,"MEASure:CURRent?")
	return (v_return,i_return)	

def getV(channel=1):
	# Ch 1 is usually top heater
	# parameters:
	ser = serial.Serial('/dev/ttyUSB1',9600,parity='N',bytesize=8,stopbits=2,timeout=2)
	output = query_serial(ser,"*IDN?")
	query_serial(ser,"SYSTem:REMote") # do we want this? or leave alone?
	query_serial(ser,"INSTrument:NSELect %s" % channel)
	v_return = query_serial(ser,"MEASure:VOLTage?")
	i_return = query_serial(ser,"MEASure:CURRent?")
	return (v_return,i_return)	
