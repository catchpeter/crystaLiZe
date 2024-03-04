# New utilities for slow control of LBNL crystaLiZe setup 
2024/03/04: initial write up - RMG

## Installation
* Requires Linux. Tested on Ubuntu 20.04
* Requires python3 with numpy, pyserial, and uldaq.
* To install uldaq, you need to first get the C library. This should be straightforward, follow instructions here: https://pypi.org/project/uldaq/

## Configure USB ports
* There are three USB ports you need to configure: cathode, gate, and heaters.
* The configuration file is "usb_ports", the options are 0,1,2 (maybe higher if other stuff gets plugged in)
* To check which port is which, the easiest way is guess-and-check with using "read_hv" and "read_heater" in shared_functions.py
* Alternatively, one can unplug and plug and run "sudo dmesg" to see the device that was just plugged in as e.g., "ttyUSB2" --> 2

## Running slow control
* The main script that records and changes heaters is "record_values.py". This should be run in the background in a screen session.
* To change the heaters, edit "heaters_setting.txt"
* The values recorded are kept in csv files located in /home/xaber/ttlogs/
* There is an email alarm system that fires whenever this script stops. To add your email, edit the function in record_values.

## Misc scripts
* To read or ramp the cathode/gate, use e.g., read_cathode.py and change_cathode.py. When changing HV, this is recorded in log files in this directory. 


