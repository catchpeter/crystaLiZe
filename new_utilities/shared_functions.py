#!/usr/bin/env python3

"""
These are basic functions for reading devices for slow control.
Code rewritten in glorious python3 from the old slow control scripts.
- RMG Feb 2024
    
uldaq package install: follow instructions https://pypi.org/project/uldaq/
You will need to first get the C library, this should be straightforward.
"""


import numpy as np
import serial
import uldaq


# =========================================
# Slow control panel
# =========================================

def read_t_p(verbose=False):

    # Device communication setup (copying example script)
    all_devices = uldaq.get_daq_device_inventory(uldaq.InterfaceType.ANY) # get list of devices
    daq_device = uldaq.DaqDevice(all_devices[0]) # there should only be one device 
    ai_device = daq_device.get_ai_device()
    descriptor = daq_device.get_descriptor()
    daq_device.connect(connection_code=0)
    ai_info = ai_device.get_info()
    input_mode = uldaq.AiInputMode.SINGLE_ENDED
    ranges = ai_info.get_ranges(input_mode)
    range_index = 0

    # Read raw values from DAQ
    ch_to_read = np.array([0,1,4,5,6,7], dtype=int)
    n_ch = max(ch_to_read) + 1
    raw_data = np.zeros(n_ch, dtype=float)
    n_samples = 100 # want to reduce noise as much as possible
    for i in range(n_samples):
        for ch in ch_to_read:
            raw_data[ch] += ai_device.a_in(ch, input_mode, ranges[range_index], uldaq.AInFlag.DEFAULT)
    raw_data /= n_samples

    # Convert to actual values
    # Conversions and variables are legacy from James (idk who that even is)
    tc = raw_data[0]*30 - 149.3
    tc5 = raw_data[5]*30 - 149.3
    tc6 = raw_data[6]*30 - 149.3
    tc7 = raw_data[7]*30 - 149.3
    pb = raw_data[1]*1.760
    all_values = np.round(np.array([tc,tc5,tc6,tc7,pb]),3)
    
    if verbose:
        print(all_values)
    
    return all_values


# =========================================
# High voltage supplies
# =========================================

def read_hv(usb_n):

    ser = serial.Serial(f'/dev/ttyUSB{usb_n}', 9600, parity='N', timeout=2)

    cmd = 'VOUT?\n'
    ser.write(cmd.encode() + serial.CR + serial.LF)
    v_raw = ser.readline()
    v = abs(float(v_raw[:-1]))
    ser.close()

    return v


def change_hv(usb_n, value):

    ser = serial.Serial(f'/dev/ttyUSB{usb_n}', 9600, parity='N', timeout=2)

    cmd = f'VSET {-value}\n'
    ser.write(cmd.encode() + serial.CR + serial.LF)
    ser.close()

    return


def check_trip(usb_n):

    ser = serial.Serial(f'/dev/ttyUSB{usb_n}', 9600, parity='N', timeout=2)

    cmd = '*STB? 2\n'
    ser.write(cmd.encode() + serial.CR + serial.LF)
    trip_status = ser.readline()
    ser.close()

    if trip_status == '1\n':
        return 1
    else:
        return 0


# =========================================
# Heater power supply
# =========================================

def read_heater(usb_n, channel):

    ser = serial.Serial(f'/dev/ttyUSB{usb_n}',9600,parity='N',bytesize=8,stopbits=2,timeout=2)

    cmd_list = ["*IDN?", "SYSTem:REMote", f"INSTrument:NSELect {channel}", "MEASure:VOLTage?", "MEASure:CURRent?"]
    for i, cmd in enumerate(cmd_list):
        ser.write(cmd.encode() + serial.CR + serial.LF)
        output = ser.readline()
        if i == 3: v_return = output
        elif i == 4: i_return = output
    ser.close()

    return [float(v_return), float(i_return)]


def change_heater(usb_n, channel, voltage):

    ser = serial.Serial(f'/dev/ttyUSB{usb_n}',9600,parity='N',bytesize=8,stopbits=2,timeout=2)

    cmd_list = ["*IDN?","SYSTem:REMote",f"INSTrument:NSELect {channel}","VOLTage:RANGe HIGH",f"VOLTage {voltage}","OUTPut ON"]
    for i, cmd in enumerate(cmd_list):
        ser.write(cmd.encode() + serial.CR + serial.LF)
        output = ser.readline()
    ser.close()

    return read_heater(usb_n, channel)


# =========================================
# Misc
# =========================================

def get_usb_n(name):

    # Yeah this one isn't very pythonic

    usb_file = np.loadtxt("usb_ports", delimiter=",", dtype=str)

    if name == "cathode" or name == "Cathode":
        usb_n = int(usb_file[1,0])
    elif name == "gate" or name == "Gate":
        usb_n = int(usb_file[1,1])
    elif name == "heaters" or name == "heater" or name == "Heaters" or name == "Heater":
        usb_n = int(usb_file[1,2])
    else:
        print("Error! Name options are cathode, gate, heaters")
        usb_n = -1

    return usb_n


