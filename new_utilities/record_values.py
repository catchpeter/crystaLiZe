#!/usr/bin/env python3

"""
Python script that records values of temps and pressure in the bkg.
I.e., nominal slow control script.
Re-written from the old version, now using the new shared_functions.py
- RMG March 2024
"""

import time
import numpy as np
import shared_functions as sf

read_period_sec = 20 # read about every 20 sec
n_reads = int(28*24*3600/read_period_sec) # four weeks of max running
n_tries = 10

nl_str = "\n=========================================\n"


def exit_email():
    print(f"{nl_str}Too many errors, exiting{nl_str}")
    return



def main():
    
    print(f"{nl_str}Starting to record values{nl_str}")

    for i in range(n_reads):

        timestr = time.strftime("%Y%m%dT%H%M",time.localtime())
        
        # In each iteration, try up to n_tries to read everything
        fail = False
        for j in range(n_tries):
            
            # Try reading temps and pressure
            try:
                t_p_arr = sf.read_t_p()
                fail = False
            except:
                print(f"{nl_str}Error in reading temperatures and pressure! Try {j+1}/{n_tries}{nl_str}")
                fail = True
                time.sleep(read_period_sec)
                continue

            # Try reading heater power supply
            try:
                v1, i1 = read_heater(heater_usb_n, 1)
                v2, i2 = read_heater(heater_usb_n, 2)
                fail = False
            except:
                print(f"{nl_str}Error in reading heaters! Try {j+1}/{n_tries}{nl_str}")
                fail = True
                time.sleep(read_period_sec)
                continue

            if not fail: 
                break
       
        # If it does fail, then email and exit
        if fail:
            exit_email()
            return

        
        # Record values in csv file



        # Change heaters settings


    

        time.sleep(read_period_sec)

    return


if __name__ == "__main__":
    main()
