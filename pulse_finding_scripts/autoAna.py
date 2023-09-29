import glob
import os
import time

from compression import compression
from rq_generate_alphas import make_rq
from read_settings import get_vscale, get_phase


def autoAna(data_dir_list, upload=False, use_old_liquid_spe=False):

    """
    Automatic data analysis.
    This should run in the background during normal TPC operation.
    For re-analyzing specific data sets, use re_ana.py
    """

    for data_dir in data_dir_list:

        # Check if data has finished transfering
        if not os.path.exists(data_dir + "transferDone"):
            continue


        # Check to compress
        if not os.path.exists(data_dir + "compressed_filtered_data/headers.npz"):
            try:
                threshold = 5
                vscale = get_vscale(data_dir)
                if vscale == (500.0/16384.0): threshold = 5
                compression(data_dir, threshold=threshold, save_everything=False)
            except:
                print("uh oh compression didn't work")


        # Check to generate rq's
        if not os.path.exists(data_dir + "rq_filtered_alphas_newSPE.npy") and not os.path.exists(data_dir + "rq_SPE_filtered_new.npy"):
            
            try:
                phase = get_phase(data_dir)
                if use_old_liquid_spe: phase = "old_liquid"
                make_rq(data_dir,phase=phase, dead=True)
            except:
                print("uh oh rq didn't work")
             

        # Send it to the cloud 
        if upload:
            command = "rclone -v copy " + data_dir +" gdrive:crystallize/data/"+ data_dir[data_dir.find("data-"):]
            os.system(command)

    return 


def main():

    while True: # lmao

        location = "/media/xaber/G-Drive2/crystalize_data/data-202309/202309*/*/" 
        data_dir_list = glob.glob(location)
        
        autoAna(data_dir_list, upload=False, use_old_liquid_spe=True)

        time.sleep(60)



if __name__ == "__main__":
    main()