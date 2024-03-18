import re
import numpy as np

def get_event_window(data_dir, old=False):

    if old:
        event_window = int(re.findall(r'_(\d+)us',data_dir)[0])
    else:
        conds = np.loadtxt(data_dir+"/conditions.csv", delimiter=",", dtype=str)
        #fields = conds[0].tolist()
        values = conds[1].tolist()
        event_window = int(values[4])

    return event_window



def get_vscale(data_dir, old=False):

    if old:
        if data_dir.find("0.5DR") != -1:
            vscale = (500.0/16384.0)
        elif data_dir.find("2DR") != -1:
            vscale = (2000.0/16384.0) # = 0.122 mV/ADCC, vertical scale
        else:
            vscale = (2000.0/16384.0) # default to 2V
    
    else:
        conds = np.loadtxt(data_dir+"/conditions.csv", delimiter=",", dtype=str)
        #fields = conds[0].tolist()
        values = conds[1].tolist()
        dr = values[2]
        if dr == "2":
            vscale = (2000.0/16384.0) 
        elif dr == "0.5":
            vscale = (500.0/16384.0)

    return vscale



def get_sipm_bias(data_dir, old=False):

    if old:
        sipm_bias = int(re.findall(r'_(\d+)SiPM',data_dir)[0])
    else:
        conds = np.loadtxt(data_dir+"/conditions.csv", delimiter=",", dtype=str)
        #fields = conds[0].tolist()
        values = conds[1].tolist()
        sipm_bias = int(values[8])

    return sipm_bias


def get_phase(data_dir, old=False):
    # Redundant 
    if old:
        phase = "liquid"
        if "solid" in data_dir or "crystal" in data_dir:
            phase = "solid"
    else:
        phase = "NA"
    
    return phase