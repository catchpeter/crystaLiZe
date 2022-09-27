import glob
import os
import time

from compression import compression
from rq_generate_32ch import make_rq
from cut_plot import make_plots




def autoAna():

    data_base_dir = "/media/xaber/extradrive1/crystalize_data/"
    all_dir = glob.glob(data_base_dir + "*/*/*/")

    # Not particularly efficient at finding data. Balanced by our beefy desktop
    for dir in all_dir:

        # Check if data has finished transfering
        if not os.path.exists(dir + "transferDone"):
            continue

        # Check to compress
        if not os.path.exists(dir + "compressed_data/compressed_0.npy"):
            try:
                compression(dir)
            except:
                print("uh oh compression didn't work")

        # Check to generate rq's and plots
        if not os.path.exists(dir + "rq.npz"):
            try:
                make_rq(dir)
                make_plots(dir)
            except:
                print("uh oh rq didn't work")

        # Send it to the cloud 
        command = "rclone -v copy " + dir +" gdrive:crystallize/data/"+ dir[dir.find("data-"):]
        os.system(command)

    return 


def main():

    while True: # lmao

        autoAna()

        # Checks every minute (roughly)
        time.sleep(60)

    return 


if __name__ == "__main__":
    main()
