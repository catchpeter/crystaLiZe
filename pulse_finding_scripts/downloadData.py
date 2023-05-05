import os


base_command = "rclone -v copy gdrive:crystallize/data/data-202207/20220715/20220715-1214_2DR_6mVtrig_15us_3201.0C_3001.0G_500A_54SiPM_0.72bar_-115.48ICVbot_2fold_CoOCVtop_120min/" 
#data_dir = "/media/xaber/extradrive4/crystalize_data/data-202211/20221105/20221105-0732_2DR_10mVtrig_30us_5202.0C_5002.0G_500A_54SiPM_1.36bar_-102.19ICVbot_2fold_degraded_afterFill_180min/"
data_dir = "/media/xaber/extradrive4/crystalize_data/data-202207/20220715/20220715-1214_2DR_6mVtrig_15us_3201.0C_3001.0G_500A_54SiPM_0.72bar_-115.48ICVbot_2fold_CoOCVtop_120min/"
n_channels = [16,8,8]
for bd in range(3):
    for ch in range(n_channels[bd]):
        command = base_command+"waveforms_"+str(bd)+"_"+str(ch)+".dat "+ data_dir

        os.system(command)

