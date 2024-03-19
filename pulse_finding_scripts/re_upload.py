"""
Simple script to upload data to google drive.
Should be run in the background, seperate from autoAna to avoid infinite upload loops in processing.
Edit data_dir_list for where it should look for data to upload.
"""

import glob
import os
import time


def re_upload():

    data_dir_list = glob.glob("/media/xaber/extradrive2/crystalize_data/data-*/*/*/")

    for data_dir in data_dir_list:

        print(data_dir)

        command = "rclone -v copy " + data_dir +" gdrive:crystallize/data/"+ data_dir[data_dir.find("data-"):]
        try:
            os.system(command)
        except:
            print("Error in upload, probably reached the limit")

    print("Done looping through data")


def main():

    while True: # lmao

        re_upload()

        # Checks every so often
        time.sleep(10)

    return 


if __name__ == "__main__":
    main()