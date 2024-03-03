#!/usr/bin/env python3

"""
Script to ramp the cathode, rewritten in python3 from old script.
- RMG March 2024

Two ways to run this script:
    1. No arguments, in which case you are prompted to input parameters 
    2. With arguments: targetV, dvdt, interal.

"""


import sys
import time
import datetime
import shared_functions as sf


def main():
    
    # Checking inputs for ramping
    n_args = len(sys.argv)
    if n_args == 1:
        targetV = abs(float(input("\nEnter voltage (V): ")))
        dvdt = abs(float(input("Enter the voltage change in each step (V) ")))
        interval = abs(int(input("Interval between each step (s): ")))
    elif n_args == 4:
        try:
            targetV = abs(float(sys.argv[1]))
            dvdt = abs(float(sys.argv[2]))
            interval = abs(float(sys.argv[3]))
        except ValueError:
            print("Parameters given are not numbers, please correct")
            sys.exit()
    else:
        print("Please give parameters in the following order: \n target voltage (V), increase step (V), interval (s).\n")
        sys.exit()

    # Initial read of HV supply
    usb_n = sf.get_usb_n(name="cathode")
    initial_value = sf.read_hv(usb_n)
    print(f"Cathode at -{initial_value} V")

    # Start log    
    log = open('cathode_ramp_log.txt', 'a')    
    log.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+'\n')
    log.write("Start to ramp cathode voltage with parameters below\n")
    log.write("Target voltage ramp to: {}\n".format(targetV))
    log.write("Voltage step: {}\n".format(dvdt))
    log.write("Time interval between voltage change: {}\n".format(interval))
    time.sleep(1)

    # Voltage ramp loop
    sign = 1
    if targetV < initial_value:
        sign = -1
    n_steps = int(abs(targetV - initial_value)/dvdt)
    for i in range(n_steps):
        
        # Change voltage
        sf.change_hv(usb_n, initial_value + i*dvdt*sign)
        time.sleep(1)

        # Check trip
        trip_status = sf.check_trip(usb_n)
        if trip_status == 1:
            log.write("The cathode is tripped, the trip voltage is: {}.\n\n".format(v_mon))
            log.close()
            sys.exit('The cathode is tripped! The tripped voltage is:{}'.format(v_mon))

        # Print voltage after each step
        current_value = sf.read_hv(usb_n)
        print(f"Cathode at -{current_value} V, {i+1}/{n_steps}")
        
        # zzzzzzzz
        time.sleep(interval-1)


    # Final check
    sf.change_hv(usb_n, targetV)
    time.sleep(1)
    final_value = sf.read_hv(usb_n)
    print(f"Cathode finally at -{final_value} V")
    log.write("The cathode did not trip, its voltage is: {}\n\n".format(targetV))
    log.close()
    

    return


if __name__ == "__main__":
    main()
