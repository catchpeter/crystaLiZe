import glob
import os 
import time
import subprocess



def transfer():

    data_base_dir = "/home/xaber/Data/"
    all_dir = glob.glob(data_base_dir + "*/*/*/")

    for dir in all_dir:

        # Check if data aquisition is finished and hasn't already been transferred 
        if not os.path.exists(dir + "readyToTransfer") or os.path.exists(dir + "transferDone"):
            continue

        print(dir)

        # Create path to directory on remote machine
        m = dir[17:28]
        d = dir[29:38]
        create_month = "ssh xaber@theo.dhcp.lbl.gov mkdir -p /media/xaber/extradrive1/crystalize_data/"+m
        os.system(create_month)
        create_day = "ssh xaber@theo.dhcp.lbl.gov mkdir -p /media/xaber/extradrive1/crystalize_data/"+m+"/"+d
        os.system(create_day)



        # Transfer data
        m_and_d = dir[17:38]
        remote_dir = "/media/xaber/extradrive1/crystalize_data/" + m_and_d
        command = "scp -r " + dir + " xaber@theo.dhcp.lbl.gov:" +remote_dir
        os.system(command)

        # Create transfer flag
        create_flag_command = "touch "+ dir + "transferDone"
        os.system(create_flag_command)

        # Transfer transfer flag
        transfer_flag_command = "scp -r " + dir + "transferDone" + " xaber@theo.dhcp.lbl.gov:/media/xaber/extradrive1/crystalize_data/" + dir[17:]
        os.system(transfer_flag_command)






    return




def main():

    while True: # hehe hoho haha

        transfer()

        # Checks every minute (roughly)
        time.sleep(60)

    return 


if __name__ == "__main__":
    main()
