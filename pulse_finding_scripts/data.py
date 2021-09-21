import os
import subprocess
import time
from shutil import copy2
from rq_generate import make_rq
from cut_plot import make_plots
from glob import glob

#time.sleep(3600*3)
#folder to have a proper spe
spe_file = "/media/xaber/gpeter/data/20210915/20210915_1057_Po_2.4g_2.6c_1.24bar_spe/spe.txt"

#set up parameters
#source = "Po"
source = "Po"
pressure = "1.24bar"
voltage = "0g_0c"
length_mins = 20
trigger = "3mv"
window = "25us"
config_file = "/home/xaber/caen/wavedump-3.8.2/data/config_self_trig_"+window+".txt"

process_flag = True # Option to process the data right after taking it.
length = str(length_mins)+"min"

#folder to save data, need to end with "/"
data_dir = "/media/xaber/gpeter/data/"  

for i in range(1):
	#make a new folder to store the data
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

	# Take data
	os.chdir("./"+new_path)
	
	subprocess.call(['cp', config_file, './'])
	subprocess.call(['wavedump', config_file, str(length_mins*60)]) 

	# Process the data
	full_path = data_dir+date+"/"+new_path+"/"
	#Place the path of this data to path.txt ready to be processed
	process_path = "/home/xaber/caen/wavedump-3.8.2/data/solid_xenon_tpc/pulse_finding_scripts/"
	copy2(spe_file, full_path)

if process_flag:
	os.chdir(process_path)
	make_rq(full_path)
	make_plots(full_path)




