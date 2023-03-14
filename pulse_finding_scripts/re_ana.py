import numpy as np 
import glob
import os

from compression import compression
from rq_generate_32ch import make_rq
from cut_plot import make_plots
from rq_generate_single_electrons import find_single_electrons


def re_ana(data_dir_list, compress=True, rq=True, plots=False, rq_se=True, upload=True):
    
    """
    Function for analyzing specific data sets, usually re-doing specific ones

    Inputs:
        data_dir_list: list of strings of data directories
        compress, rq, plots, rq_se, upload: bools for what you want to do 
    """

    for data_dir in data_dir_list:

        if not os.path.exists(data_dir + "transferDone"):
            print("\n\n Data not finished transfering in "+data_dir+"\n\n")
            continue

        if compress:
            try:
                print("\n\nCompressing in "+data_dir+"\n\n")
                compression(data_dir)
            except:
                print("Error in compressing "+data_dir)
        
        if rq:
            try:
                print("\n\nFinding rq's in "+data_dir+"\n\n")
                make_rq(data_dir)
            except:
                print("\n\nError finding rq's in "+data_dir+"\n\n")

        if plots:
            try:
                print("\n\nMaking plots in "+data_dir+"\n\n")
                make_plots(data_dir)
            except:
                print("\n\nError making plots in "+data_dir+"\n\n") 

        if rq_se:
            try:
                print("\n\nFinding SE's in "+data_dir+"\n\n")
                find_single_electrons(data_dir)
            except:
                print("\n\nError finding SE's in "+data_dir+"\n\n")

        if upload:
            try:
                command = "rclone -v copy " + data_dir +" gdrive:crystallize/data/"+ data_dir[data_dir.find("data-"):]
                os.system(command)
            except:
                print("\n\nError uploading "+data_dir+"\n\n")

    return



def main():

    location = "/media/xaber/f5d91b31-9a7d-3278-ac5b-4f9ae16edd60/crystalize_data/data-202303/20230308/*Co*/"
    #location = "/media/xaber/extradrive1/crystalize_data/data-202303/20230308/*Co*/"
    data_dir_list = glob.glob(location)

    re_ana(data_dir_list, compress=True, rq=True, plots=False, rq_se=True, upload=True)

    return


if __name__ == "__main__":
    main()