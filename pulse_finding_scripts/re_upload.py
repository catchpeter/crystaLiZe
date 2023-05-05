import glob
import os
import time


def re_upload():

    #data_dir_list = glob.glob("/media/xaber/f5d91b31-9a7d-3278-ac5b-4f9ae16edd60/crystalize_data/data-*/*/*/")
    data_dir_list = glob.glob("/media/xaber/G-Drive2/crystalize_data/data-*/*/*/")

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

        # Checks every minute (roughly)
        time.sleep(5*60)

    return 


if __name__ == "__main__":
    main()