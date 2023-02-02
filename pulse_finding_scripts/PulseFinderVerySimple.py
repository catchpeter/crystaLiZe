import numpy as np

tscale = (8.0/4096.0) 

def PulseFinderVerySimple(waveform, verbose=False):

    wsize = waveform.size

    start_time = 0
    end_time = 0

    max_loc = np.argmax(waveform)
    max_val = max(waveform)
    diff_wf = np.diff(waveform)


    # If pulse is too close to left edge of wf or too small, toss it
    if max_loc < 100: 
        if verbose: print("Pulse too close to waveform edges")
        return start_time, end_time
    elif max_val < 3*0.013: # this is about 2.5 phd
        if verbose: print("Pulse height too small")
        return start_time, end_time


    # Find left bound
    lht = 0.005
    try:
        below_lht_i = (diff_wf[:max_loc] < lht).nonzero()[0]
        lb = max_loc - min(max_loc - below_lht_i)
    except:
        if verbose: print("Cannot find left bound")
        return start_time, end_time

    i = 1
    while waveform[lb] > 0.04*max_val:
        try:
            below_lht_i = (diff_wf[:(max_loc-i*40)] < lht).nonzero()[0]
            lb = max_loc - min(max_loc - below_lht_i)
            i+=1
        except:
            i+=1
        if i > 100: break


    # Find the right bound
    rht = 0 #0.005 #0.01 
    below_rht_i = (diff_wf[max_loc:] > -rht).nonzero()[0]
    if below_rht_i.size == 0:
        rb = wsize-1
        if verbose: print("Cannot find right bound")
    else:
        rb = 2*max_loc + min(below_rht_i - max_loc)

    i = 1
    while waveform[rb] > 0.05*max_val and below_rht_i.size != 0:
        try:
            below_rht_i = (diff_wf[(max_loc+i*100):] > -rht).nonzero()[0]
            rb = 2*max_loc + min(below_rht_i - max_loc) + i*100
            i+=1
        except:
            i+=1
            break
        if i > 100: break

    
    # Basic width cutter
    if rb - lb > 3000:
        if verbose: print("Pulse width is too large")
        return start_time, end_time
    elif rb - lb < 15:
        if verbose: print("Pulse width is too small")
        return start_time, end_time


    start_time = lb
    end_time = rb

    return start_time, end_time