import os
import subprocess
import numpy as np

def baseline_set():
	n_sipms = 8
	#threshold in mV
	threshold = np.zeros(n_sipms)
	threshold[0] = 3
	threshold[1] = 3
	threshold[2] = 10
	threshold[3] = 3
	threshold[4] = 3
	threshold[5] = 3
	threshold[6] = 3
	threshold[7] = 10



	event_window = 25
	config_file = "/home/xaber/caen/wavedump-3.8.2/data/config_self_trig_25us.txt"



	#folder to save data, need to end with "/"
	data_dir = "/media/xaber/gpeter/data/"  


	os.chdir(data_dir)
	try:
		os.mkdir("baseline")
	except:
		pass
	os.chdir("./baseline")

	subprocess.call(['wavedump', config_file, '6']) 

	# Calculate the baseline
	wsize = int(500 * event_window)  # samples per waveform # 12500 for 25 us
	vscale = (2000.0/16384.0) # = 0.122 mV/ADCC, vertical scale
	tscale = (8.0/4096.0)     # = 0.002 µs/sample, time scale

	save_avg_wfm = False # get the average waveform passing some cut and save to file

	post_trigger = 0.5 # Was 0.2 for data before 11/22/19
	trigger_time_us = event_window*(1-post_trigger)
	trigger_time = int(trigger_time_us/tscale)

	ch_data = []
	load_dtype = "int16"

	for ch_ind in range(n_sipms):
		ch_data.append(np.fromfile(data_dir + "./baseline/wave"+str(ch_ind)+".dat", dtype=load_dtype))

	baseline = np.zeros(n_sipms)
	for ch_ind in range(n_sipms):
	    V = ch_data[ch_ind]
	    V = V[:int(len(V) / wsize) * wsize]
	    V = V.reshape(int(V.size / wsize), wsize) # reshape to make each channel's matrix of events
	    V_average = np.average(V, axis = 0)   #average all events taken
	    baseline[ch_ind] = (np.average(V_average[:int(wsize*(1-post_trigger)/4)]))  #average the 1/4 of all pre-trigger sampels

	#calculate triggers
	triggers = (baseline+threshold/vscale)
	triggers = np.round(triggers)
	triggers = triggers.astype(int)

	#edit the config file
	with open(config_file, "r") as f:
	    a = f.readlines()

	#find lines needs to be changed
	index_to_change = []
	for i in range(len(a)):
	    if a[i].startswith("TRIGGER_THRESHOLD    "): index_to_change.append(i)

	for i in range(8):
	    a[index_to_change[i]] = 'TRIGGER_THRESHOLD      {}\n'.format(triggers[i])

	with open(config_file, "w") as f:
	    for text_line in a:
	        f.write(text_line)

	with open("/home/xaber/caen/wavedump-3.8.2/data/baseline.txt", "w") as f:
		for n in baseline:
			f.write("{}\n".format(n))

def main():
    baseline_set()

if __name__ == "__main__":
    main()