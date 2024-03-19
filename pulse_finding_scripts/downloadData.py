"""
Simple script to grab data from google drive
"""

import os

# Data location on google drive
base_command = "rclone -v copy gdrive:crystallize/data/data-202210/20221004/20221004-1538_2DR_20mVtrig_15us_3202.0C_3001.0G_500A_54SiPM_1.22bar_-97.36ICVbot_2fold_CoOCVtop_20min/"

# Local location
data_dir = "/media/xaber/G-Drive2/crystalize_data/data-202210/20221004/20221004-1538_2DR_20mVtrig_15us_3202.0C_3001.0G_500A_54SiPM_1.22bar_-97.36ICVbot_2fold_CoOCVtop_20min/"

command = base_command + " " + data_dir
os.system(command)
