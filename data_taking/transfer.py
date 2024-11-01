import glob
import os 
import time
import subprocess



def transfer():

    # f5d91b31-9a7d-3278-ac5b-4f9ae16edd60

    remote_drive_name = "extradrive1"

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
        create_month = f"ssh xaber@theo.dhcp.lbl.gov mkdir -p /media/xaber/{remote_drive_name}/crystalize_data/"+m
        os.system(create_month)
        create_day = f"ssh xaber@theo.dhcp.lbl.gov mkdir -p /media/xaber/{remote_drive_name}/crystalize_data/"+m+"/"+d
        os.system(create_day)



        # Transfer data
        m_and_d = dir[17:38]
        remote_dir = f"/media/xaber/{remote_drive_name}/crystalize_data/" + m_and_d
        command = "scp -r " + dir + " xaber@theo.dhcp.lbl.gov:" +remote_dir
        os.system(command)

        # Create transfer flag
        create_flag_command = "touch "+ dir + "transferDone"
        os.system(create_flag_command)

        # Transfer transfer flag
        transfer_flag_command = "scp -r " + dir + "transferDone" + f" xaber@theo.dhcp.lbl.gov:/media/xaber/{remote_drive_name}/crystalize_data/" + dir[17:]
        os.system(transfer_flag_command)






    return




def main():

    while True: # hehe hoho haha

        transfer()

        # Checks every minute (roughly)
        time.sleep(5)

    return 


if __name__ == "__main__":
    main()
