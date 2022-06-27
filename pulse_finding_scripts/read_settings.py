

def get_event_window(data_dir):

    event_window = -999

    for i in range(1,201):
        if data_dir.find("_"+str(i)+"us") != -1:
            event_window = i

    return event_window


def get_vscale(data_dir):

    if data_dir.find("0.5DR") != -1:
        vscale = (500.0/16384.0)
    elif data_dir.find("2DR") != -1:
        vscale = (2000.0/16384.0) # = 0.122 mV/ADCC, vertical scale
    else:
        vscale = (2000.0/16384.0) # default to 2V

    return vscale