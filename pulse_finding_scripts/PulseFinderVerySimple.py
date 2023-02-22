import numpy as np

tscale = (8.0/4096.0) 

def PulseFinderVerySimple(waveform, lht=0.005, rht=0.0001, verbose=False):

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


    # =============================================================================================
    # Find left bound

    # Use threshold on derivative first 
    try:
        below_lht_i = (diff_wf[:max_loc] < lht).nonzero()[0]
        lb = max_loc - min(max_loc - below_lht_i)
    except:
        if verbose: print("Cannot find left bound")
        return start_time, end_time

    # Make sure derivative didn't find the wrong local minima
    # Use threshold based on fraction of max height of pulse
    if max_val < 20:
        if max_val < 0.25: 
            lht_wf = 0.1*max_val # small S1's
        else: 
            lht_wf = 0.01*max_val # Medium S1's and small S2's
    else:
        lht_wf = 0.04*max_val # Large S1's and S2's
    i = 1
    while waveform[lb] > lht_wf:
        try:
            below_lht_i = (diff_wf[:(max_loc-i*10)] < lht).nonzero()[0]
            lb = max_loc - min(max_loc - below_lht_i)
            i+=1
        except:
            i+=1
        if i > 100: break
    
    # Make sure left bound isn't negative
    if waveform[lb] < 0:
        try: 
            lb = lb + (waveform[lb:max_loc] > 0).nonzero()[0][0]
        except:
            if verbose: print("Error in getting left bound above zero")

    
    # =============================================================================================
    # Find the right bound

    # Use threshold on derivative first 
    below_rht_i = (diff_wf[max_loc:] > -rht).nonzero()[0]
    if below_rht_i.size == 0:
        rb = wsize-1
        if verbose: print("Cannot find right bound")
    else:
        rb = 2*max_loc + min(below_rht_i - max_loc)

    # Make sure derivative didn't find the wrong local minima
    # Use threshold based on fraction of max height of pulse
    if max_val < 20 and max_val > 0.1: # Small S2's
        rht_wf = 0.01*max_val
    else: # Small S1's, large S2's
        rht_wf = 0.05*max_val 
    i = 1
    while waveform[rb] > rht_wf and below_rht_i.size != 0:
        try:
            below_rht_i = (diff_wf[(max_loc+i*100):] > -rht).nonzero()[0]
            rb = 2*max_loc + min(below_rht_i - max_loc) + i*100
            i+=1
        except:
            i+=1
        if i > 100: break

    # Make sure left bound isn't negative
    if waveform[rb] < 0:
        try: 
            rb = max_loc + (waveform[max_loc:rb] > 0).nonzero()[0][-1]
        except:
            if verbose: print("Error in getting left bound above zero")
    


    # =============================================================================================
    # Final cuts

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