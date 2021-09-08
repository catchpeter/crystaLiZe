import os
import subprocess
import time
from shutil import copy2

#folder to have a proper spe
spe_file = "/media/xaber/gpeter/data/20210810/20210810_1622_Co_OCVtop_1.38_spe/spe.txt"

source = "Po_Co_OCVtop"
pressure = "1.3bar"
voltage = "2.8g_3.0c"
length_mins = 100
trigger = "3mv"
window = "3us"
config_file = "/home/xaber/caen/wavedump-3.8.2/data/config_self_trig_25us.txt"

length = str(length_mins)+"min"

#folder to save data, need to end with "/"
data_dir = "/media/xaber/gpeter/data/"  


for i in range(1):
	date = time.strftime("%Y%m%d", time.localtime())
	clock = time.strftime("%H%M", time.localtime())

	os.chdir(data_dir)
	try:
		os.mkdir(date)
	except:
		pass
	os.chdir("./"+date)

	try:
		new_path = "_".join((date, clock, source, voltage, pressure, trigger, window, "circ", length))
		os.mkdir(new_path)
	except:
		pass
	os.chdir("./"+new_path)
	
	#subprocess.call(['wavedump', config_file, str(length_mins*60)]) 

	print("Done, good")
	full_path = data_dir+date+"/"+new_path+"/"
	#Place the path of this data to path.txt ready to be processed
	process_path = "/home/xaber/caen/wavedump-3.8.2/data/solid_xenon_tpc/pulse_finding_scripts/"
	copy2(spe_file, full_path)

	os.chdir(process_path)
	with open("path.txt", "w") as f:
		f.write(full_path)




