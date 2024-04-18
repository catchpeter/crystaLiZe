#!/usr/bin/env python3

"""
Python script that records values of temps and pressure in the bkg.
I.e., nominal slow control script.
Re-written from the old version, now using the new shared_functions.py
- RMG March 2024
"""

import os
import time
import numpy as np
import shared_functions as sf

read_period_sec = 20 # read about every 20 sec
n_reads = int(28*24*3600/read_period_sec) # four weeks of max running
n_tries = 10

nl_str = "\n=========================================\n"

heater_usb_n = sf.get_usb_n(name="heater")


def exit_email(fail=False):

    if fail:
        email_str = "ERROR: record values has stopped!"
    else:
        email_str = "Record values has concluded without errors"

    print(f"{nl_str}{email_str}{nl_str}")
    save_dir = "/home/xaber/crystaLiZe/new_utilities/exit_email_record_values.txt"
    np.savetxt(save_dir, np.array(["To: rmg@lbl.gov","Subject: Record values has stopped","From: crystalize_slow_control@xena.dhcp.lbl.gov", "" ,email_str]), fmt="%s")
    os.system("sendmail -t < /home/xaber/crystaLiZe/new_utilities/exit_email_record_values.txt")

    return



def record_values():
    
    current_top_power = -1
    current_bottom_power = -1

    print(f"{nl_str}Starting to record values{nl_str}")
    fstr = time.strftime("%Y%m%dT%H%M",time.localtime())
    start_sec = time.time()

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
                v1, i1 = sf.read_heater(heater_usb_n, 1)
                v2, i2 = sf.read_heater(heater_usb_n, 2)
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
            exit_email(fail)
            return

        # Create list of values to save
        # Ordering is legacy...
        ymd_str = time.strftime("%Y%m%d",time.localtime())
        hms_str = time.strftime("%H%M%S",time.localtime())
        elapsed_time = time.time()-start_sec
        tc4 = t_p_arr[0]
        tc5 = t_p_arr[1]
        tc6 = t_p_arr[2]
        tc7 = t_p_arr[3]
        pb = t_p_arr[4]
        pow_top = v1*i1
        pow_bot = v2*i2
        save_list = (ymd_str, hms_str, elapsed_time, tc4, tc5, tc6, tc7, pb, pow_top, pow_bot)

        print ("%s : T4=%3.3f C, T5=%3.3f C, T6=%3.3f C, T7=%3.3f C, P=%1.3f Bar, iteration=%d" % (timestr,tc4,tc5,tc6,tc7,pb,i))
        
        # Save values 
        with open(("/home/xaber/ttlogs/%s.csv" % fstr),"a+") as fid:
            fid.write("%s,%s,%f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f\n" % save_list)
        # Also save to tmp file for remote read out
        with open("/home/xaber/ttlogs/current.csv", "w") as ftemp:
            ftemp.write("%s,%s,%f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f,%3.3f\n" % save_list)

        # Read heaters settings
        try:
            with open("/home/xaber/crystaLiZe/new_utilities/heaters_setting.txt") as f:
                f.readline()
                power_line = f.readline()
                powers = power_line.split(",")
            top_heater_power = float(powers[0])
            bottom_heater_power = float(powers[1])

#             # Change heaters settings
#             if abs(top_heater_power-current_top_power)>0.005:
#                 print("the current top heater power is {:.2f}W, will change to {:.2f}W".format(current_top_power, top_heater_power))
#                 top_heater_v = np.sqrt(25.*top_heater_power)
#                 (v_top_now, i_top_now) = sf.change_heater(heater_usb_n, 1, top_heater_v)
#                 current_top_power = top_heater_power
#             if abs(bottom_heater_power-current_bottom_power)>0.005:
#                 print("the current bottome heater power is {:.2f}W, will change to {:.2f}W".format(current_bottom_power, bottom_heater_power))
#                 bot_heater_v = np.sqrt(25.*bottom_heater_power)
#                 (v_bottom_now, i_bottom_now) = sf.change_heater(heater_usb_n, 2, bot_heater_v)
#                 current_bottom_power= bottom_heater_power
        except:
            print("error in reading/writing")

        # Sleep until next iteration
        time.sleep(read_period_sec)


    # End of line
    exit_email()

    return


def main():

    try:
        record_values()
    except:
        exit_email(fail=True)

    return



if __name__ == "__main__":
    main()
