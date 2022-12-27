import re

def get_event_window(data_dir):
    event_window = int(re.findall(r'_(\d+)us',data_dir)[0])
    return event_window

def get_vscale(data_dir):
    if data_dir.find("0.5DR") != -1:
        vscale = (500.0/16384.0)
    elif data_dir.find("2DR") != -1:
        vscale = (2000.0/16384.0) # = 0.122 mV/ADCC, vertical scale
    else:
        vscale = (2000.0/16384.0) # default to 2V
    
    return vscale