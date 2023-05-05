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
                make_rq(data_dir,handscan=False)
            except:
                print("\n\nError finding rq's in "+data_dir+"\n\n")

        if plots:
            try:
                print("\n\nMaking plots in "+data_dir+"\n\n")
                make_plots(data_dir)
            except:
                print("\n\nError making plots in "+data_dir+"\n\n") 

        if rq_se:
            #try:
            #    print("\n\nFinding SE's in "+data_dir+"\n\n")
            find_single_electrons(data_dir,phase="solid")
            #except:
            #    print("\n\nError finding SE's in "+data_dir+"\n\n")

        if upload:
            try:
                command = "rclone -v copy " + data_dir +" gdrive:crystallize/data/"+ data_dir[data_dir.find("data-"):]
                os.system(command)
            except:
                print("\n\nError uploading "+data_dir+"\n\n")

    return



def main():

    #location = "/media/xaber/extradrive2/crystalize_data/data-202303/20230321/*/"
    #location = "/media/xaber/f5d91b31-9a7d-3278-ac5b-4f9ae16edd60/crystalize_data/data-202211/202211*/*degraded*/"
    #location = "/media/xaber/G-Drive2/crystalize_data/data-202304/20230418/20230418-0945_2DR_10mVtrig_15us_3201.0C_3001.0G_500A_54SiPM_0.7bar_-118.5ICVbot_2fold_plainMesh_solid_252Cf_extraPb_60min/"
    #data_dir_list = ["/media/xaber/extradrive4/crystalize_data/data-202207/20220715/20220715-1214_2DR_6mVtrig_15us_3201.0C_3001.0G_500A_54SiPM_0.72bar_-115.48ICVbot_2fold_CoOCVtop_120min/"]
    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202304/20230425/20230425-0954_2DR_10mVtrig_15us_3202.0C_3001.0G_500A_54SiPM_0.6bar_-121.52ICVbot_2fold_plainMesh_rightAfterFreeze_CoOCVTop_20min/"]
    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202304/20230425/20230425-1533_2DR_10mVtrig_15us_3202.0C_3001.0G_500A_54SiPM_0.58bar_-121.52ICVbot_2fold_plainMesh_rightAfterFreeze_CoOCVTop_20min/"]
    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202304/20230425/20230425-1558_2DR_10mVtrig_50us_3202.0C_3001.0G_500A_54SiPM_0.6bar_-120.92ICVbot_2fold_plainMesh_rightAfterFreeze_CoOCVTop_20min/"]
    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202304/20230425/20230425-1619_2DR_10mVtrig_100us_3202.0C_3001.0G_500A_54SiPM_0.64bar_-120.92ICVbot_2fold_plainMesh_SE_1msdelay_CoOCVTop_20min/"]
    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202304/20230426/20230426-1055_2DR_10mVtrig_20us_3202.0C_3001.0G_500A_54SiPM_0.64bar_-121.22ICVbot_2fold_plainMesh_SE_1msdelay_CoOCVTop_20min/"]
    
    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202304/20230426/20230426-1801_2DR_10mVtrig_50us_3202.0C_3001.0G_500A_54SiPM_0.64bar_-120.01ICVbot_2fold_plainMesh_SE_1msdelay_CoOCVTop_20min/"]
    
    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202304/20230427/20230427-1052_2DR_10mVtrig_20us_3201.0C_3001.0G_500A_54SiPM_0.58bar_-121.52ICVbot_2fold_plainMesh_SE_1msdelay_CoOCVTop_20min/"]
    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202304/20230427/20230427-1508_2DR_10mVtrig_20us_3202.0C_3001.0G_1000A_54SiPM_0.62bar_-120.92ICVbot_2fold_plainMesh_SE_500usdelay_CoOCVTop_20min/"]

    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202304/20230428/20230428-1621_2DR_10mVtrig_50us_3202.0C_3001.0G_1000A_54SiPM_0.69bar_-119.71ICVbot_2fold_plainMesh_SE_1p5msdelay_CoOCVTop_20min/"]

    data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202304/20230428/20230428-1640_2DR_10mVtrig_50us_3202.0C_3001.0G_500A_54SiPM_0.65bar_-119.71ICVbot_2fold_plainMesh_SE_1p5msdelay_CoOCVTop_20min/"]
    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202304/20230428/20230428-1632_2DR_10mVtrig_50us_3202.0C_3001.0G_1000A_54SiPM_0.67bar_-119.71ICVbot_2fold_plainMesh_SE_1p5msdelay_CoOCVTop_20min/"]
    
    
    #data_dir_list = ["/media/xaber/G-Drive2/crystalize_data/data-202305/20230501/20230501-1139_2DR_10mVtrig_50us_3201.0C_3001.0G_500A_54SiPM_0.65bar_-119.71ICVbot_2fold_plainMesh_CoOCVTop_120min/"]

    #data_dir_list = glob.glob(location)

    re_ana(data_dir_list, compress=False, rq=False, plots=False, rq_se=True, upload=False)

    return


if __name__ == "__main__":
    main()