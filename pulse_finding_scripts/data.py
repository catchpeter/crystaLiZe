import os
import subprocess
import time
from shutil import copy2
from rq_generate import make_rq
from cut_plot import make_plots
from glob import glob
from baseline import baseline_set
from grid_voltage import grid_voltage
from read_pressure import read_pressure
#time.sleep(3600*3)

#set up parameters
#source = "Po"
data_plan = open("data_plan.txt", 'r')
source = data_plan.readline().strip()
#folder to have a proper spe
spe_file = data_plan.readline().strip()
info_help = data_plan.readline()

#read in the data_plan.txt, and make a list of tasks
voltage_sets = []

for line in data_plan: 
	try:
		line_int = list(map(int, line.strip().split(",")))
	except:
		break

	for i in range(line_int[3]):
		voltage_sets.append(line_int[:3])

pressure = "{}bar".format(read_pressure())
trigger = "3mv"
window = "25us"
config_file = "/home/xaber/caen/wavedump-3.8.2/data/config_self_trig_"+window+".txt"

process_flag = True # Option to process the data right after taking it.

#folder to save data, need to end with "/"
data_dir = "/media/xaber/gpeter/data/"  



for one_data in voltage_sets:
	#set voltage and time
	gate_voltage = one_data[0]
	cathode_voltage = one_data[1]
	length_mins = one_data[2]
	length = str(length_mins)+"min"
	voltage = "{:.1f}g_{:.1f}c".format(gate_voltage/1000., cathode_voltage/1000.)

	#Ramp up voltage
	ramp_flag = grid_voltage(gate = gate_voltage, cathode = cathode_voltage)
	time.sleep(120)

	#set up the baseline and triggers
	baseline_set()

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
	
	#record if the grid voltages are correct. 
	with open("note.txt", "w") as f:
		f.write("The grid voltage is set up: "+str(ramp_flag))

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




