import numpy as np
import time
import glob
import datetime as dt
import os


def recordAcquisitions(start_time):

    # Tracking spreadsheet file name
    run_name = "24B"
    csv_file_name = f"/home/xaber/crystalize/Analysis/acquisitions/run_{run_name}.xlsx"

    # If it exists, get the final line, which is the start time to start recording new entries
    csv_exists = os.path.isfile(csv_file_name)
    if csv_exists:
        with open(csv_file_name, "r") as f:
            last_line = f.readlines()[-1]
            last_dt = last_line[:15]
            start_time = dt.datetime(int(last_dt[:4]), int(last_dt[4:6]), int(last_dt[6:8]), int(last_dt[9:11]), int(last_dt[11:13]), int(last_dt[13:15]))
    
    # Where to look for data
    data_dir_base = "/media/xaber/outSSD2/crystalize_data/"
    data_dir_base_2 = "/media/xaber/extradrive1/crystalize_data/"
    long_list = sorted(glob.glob(data_dir_base+"*/*/*/"))+sorted(glob.glob(data_dir_base_2+"*/*/*/"))

    # Loop over data dirs
    for i, data_dir in enumerate(long_list):

        # Again, does the tracking spreadsheet exist?
        csv_exists = os.path.isfile(csv_file_name)

        # Get the datetime
        dt_stamp = data_dir[-16:-1]
        try:
            dt_dt = dt.datetime(int(dt_stamp[:4]), int(dt_stamp[4:6]), int(dt_stamp[6:8]), int(dt_stamp[9:11]), int(dt_stamp[11:13]), int(dt_stamp[13:15]))
        except:
            continue

        # Should we record it?
        if (dt_dt - start_time).total_seconds() > 1:
            print(data_dir)

            # Open the conditions file
            try:
                conds = np.loadtxt(data_dir+"conditions.csv", delimiter=",", dtype=str)
            except:
                break
            fields = conds[0].tolist()
            values = conds[1].tolist()
            
            # Add the location on theo
            fields.append("Location")
            values.append(data_dir)

            # If it doesn't exist, record the fields as well
            if not csv_exists:
                with open(csv_file_name, "w") as f:
                    np.savetxt(f, [fields,values], delimiter=",", fmt="%s")
            # Otherwise, record just the values
            else:
                with open(csv_file_name, "a") as f:
                    np.savetxt(f, [values], delimiter=",", fmt="%s")


    # Send it to the cloud
    command = "rclone -v copy " + csv_file_name +" gdrive:crystallize/data/tracking/acquisitions/"
    csv_exists = os.path.isfile(csv_file_name)
    if csv_exists:
        try:
            os.system(command)
        except:
            print("Error in upload, probably reached the limit")


    
    return



def main():

    start_time = dt.datetime(2024,3,6,0,0,0)

    while True:

        recordAcquisitions(start_time)

        time.sleep(5*60)


    return


if __name__ == "__main__":
    main()