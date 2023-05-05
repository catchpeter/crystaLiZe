import glob
import os
import time

from compression import compression
from rq_generate_32ch import make_rq
from cut_plot import make_plots




def autoAna(data_dir_list, upload=True):

    """
    Automatic data analysis.
    This should run in the background during normal TPC operation.
    For analyzing specific data sets, use re_ana.py
    """

    for data_dir in data_dir_list:

        print(data_dir)

        # Check if data has finished transfering
        if not os.path.exists(data_dir + "transferDone"):
            continue

        # Check to compress
        if not os.path.exists(data_dir + "compressed_filtered_data"):
            try:
                compression(data_dir)
            except:
                print("uh oh compression didn't work")

        # Check to generate rq's and plots
        if not os.path.exists(data_dir + "rq_filtered.npy"):
            try:
                make_rq(data_dir,phase="solid")
                #make_plots(dir)
            except:
                print("uh oh rq didn't work")

        # Send it to the cloud 
        if upload:
            command = "rclone -v copy " + data_dir +" gdrive:crystallize/data/"+ data_dir[data_dir.find("data-"):]
            os.system(command)

    return 


def main():

    # f5d91b31-9a7d-3278-ac5b-4f9ae16edd60
    

    while True: # lmao

        location = "/media/xaber/G-Drive2/crystalize_data/data-202305/*/*/" 
        data_dir_list = glob.glob(location)

        autoAna(data_dir_list, upload=False)

        # Checks every minute (roughly)
        time.sleep(5*60)

    return 


if __name__ == "__main__":
    main()
