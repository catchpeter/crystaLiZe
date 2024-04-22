#!/usr/bin/env python3
"""
Script to freeze the detector. 
Run this in the background in parallel with record_values.py.
Rewritten from old script
2023/03/20 - RMG
"""

import os
import csv
import time
import numpy as np


def exit_email(success):

    if success:
        email_str = "Freezing has concluded successfully!"
    else:
        email_str = "ERROR: freezing has stopped!"

    #print(f"{nl_str}{email_str}{nl_str}")
    save_dir = "/home/xaber/crystaLiZe/new_utilities/exit_email_freezing.txt"
    np.savetxt(save_dir, np.array(["To: rmg@lbl.gov","Subject: Freezing has stopped","From: crystalize_slow_control@xena.dhcp.lbl.gov", "" ,email_str]), fmt="%s")
    os.system("sendmail -t < /home/xaber/crystaLiZe/new_utilities/exit_email_record_values.txt")

    return


def write_heater_file(filename,top,bot):
    with open(filename,"w") as f: 
        f.write("top, bot\n%1.2f, %1.2f\n" % (top,bot))
        f.close()
    return


def freezing(top_heater_start):

    print("Starting to freeze")

    dt = 2*3600 # time interval for lowering heater
    dPdt = 0.01 # change in heater power every time interval
    top_heater_low = 1 # lower limit on heater power
    heater_file = "/home/xaber/crystaLiZe/new_utilities/heaters_setting.txt"
    n_tries = 1000 # max amount of attempts to write, shouldn't be an issue
    top_heater_end = 1.5 # stable heater power
    bot_heater_end = 0.5 # stable heater power

    n_steps_max = int( (top_heater_start - top_heater_low)/dPdt )
    for i in range(n_steps_max):

        # Check pressure. If we're good, set heaters to stable and exit
        current_csv = "/home/xaber/ttlogs/current.csv"
        with open(current_csv) as f:
            csvFile = csv.reader(f)
            for lines in csvFile:
                PB = float(lines[7]) # latest read of pressure (Bar)
        if PB < 0.78:
            print("Freezing has concluded! Exiting")
            write_heater_file(heater_file, top_heater_end, bot_heater_end)
            exit_email(success=True)
            return

        # Next step for the top heater
        top_heater_i = top_heater_start - (i+1)*dPdt
        print(f"Changing top heater to {top_heater_i}")

        # Write to heaters_settings
        error_flag = False
        for j in range(n_tries):
            try:
                write_heater_file(heater_file, top_heater_i, 0)
                error_flag = False
                break
            except:
                print(f"Error in trying to change heaters. Try {j+1}/{n_tries}")
                error_flag = True
                time.sleep(10)
        if error_flag == True:
            print("Too many errors. Exiting")
            exit_email(success=False)
            return

        time.sleep(dt)

    # If you made it all the way through without finishing
    exit_email(success=False)

    return



def main():

    top_heater_start = 3.7
    freezing(top_heater_start)

    return


if __name__ == "__main__":
    main()
