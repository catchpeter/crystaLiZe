import numpy as np

# functions to calculate interesting quantities for pulses
# inputs: pulse start and end samples and baseline subtracted waveform numpy array
# outputs: area, max height, max height sample, min height, mean and RMS between the start and end times

def GetPulseArea( p_start, p_end, waveforms_bls ):
    area = -999.
    try:
        area = np.sum( waveforms_bls[p_start:p_end] )
    except ValueError:
        area = -999.
    return area

############################################################

def GetPulseAreaChannel(p_start, p_end, waveform_bls):
    areas = np.zeros_like(waveform_bls)
    try:
        areas = np.sum(waveform_bls[:,p_start:p_end],axis=1)
        return areas
    except:
        return areas

############################################################

def GetPulseMaxSample( p_start, p_end, waveforms_bls ):
    max_height_sample = int(-999)
    if p_start != p_end:
        try:
            max_height_sample = p_start + np.argmax( waveforms_bls[p_start:p_end] )
        except ValueError:
            max_height_sample = int(-999)
    else:
        max_height_sample = int(-999)
    return max_height_sample

############################################################

def GetPulseMaxHeight( p_start, p_end, waveforms_bls ):
    max_height = -999.
    max_height_sample = int(-999)
    if p_start != p_end:
        try:
            max_height_sample = p_start + np.argmax( waveforms_bls[p_start:p_end] )
            max_height = waveforms_bls[max_height_sample]
        except ValueError:
            max_height_sample = int(-999)
            max_height = -999.
    else:
        max_height_sample = int(-999)
        max_height = -999.
    return max_height

############################################################

def GetPulseMaxHeightChannel(p_start, p_end, waveforms_bls):
    height_ch = []
    for i in range(32):
        try:
            height_ch.append(GetPulseMaxHeight(p_start,p_end,waveforms_bls[i,:]))
        except ValueError:
            height_ch.append(0)
    return height_ch

############################################################

def GetPulseMinHeight( p_start, p_end, waveforms_bls ):
    min_height = 999.
    try:
        min_height = np.min( waveforms_bls[p_start:p_end] )
    except ValueError:
        min_height = 999.
    return min_height

############################################################

def GetPulseMeanAndRMS( p_start, p_end, waveforms_bls ):
    
    mean = int(-999)
    rms = int(-999)
    
    try:
        # Calculate mean and rms times:
        
        # First step is to make an array that is the pulse just between the start and end times
        pulse_waveform = waveforms_bls[p_start:p_end]
        
        # Now make an array that contains the sample numbers of the found pulse
        # i.e. its first value is the pulse start and the last is the pulse end
        #pulse_sample_nums = np.arange(p_start_sum_le[0],p_end_sum_le[0]+1)
        pulse_sample_nums = np.arange( p_end - p_start )
        
        # Mean time:
        mean = np.dot( pulse_waveform,pulse_sample_nums ) / np.sum( pulse_waveform )
        
        #RMS time:
        rms = np.sqrt( np.dot( pulse_waveform, np.square(pulse_sample_nums-mean) ) / np.sum( pulse_waveform ) )
        
        #print "start = %d" % p_start
        #print "end = %d" % p_end
        #print "mean = %f" % mean
        #print "rms = %f" % rms
        
    except ValueError:
        mean = int(-999)
        rms = int(-999)
    
    return mean, rms    

############################################################

# function to calculate samples at 10% and 50% max height point looking from left and the right
def GetHeightFractionSamples( p_start, p_end, waveforms_bls ):
    
    hfs_10l = int(999999) # looking from left
    hfs_50l = int(999999) # looking from left
    hfs_90l = int(999999) # looking from left
    hfs_10r = int(-999999) # looking from right
    hfs_50r = int(-999999) # looking from right
    hfs_90r = int(-999999) # looking from right
    
    p_max_height = GetPulseMaxHeight( p_start, p_end, waveforms_bls )

    try:
        # only consider the part of the waveform in this pulse, start of pulse is sample 0
        height_fractions = waveforms_bls[int(p_start):int(p_end)]/p_max_height
        # get first and last index where height is greater than x% of max
        hfs_10 = np.argwhere(height_fractions >= 0.10)
        hfs_10l = hfs_10[0]
        hfs_10r = hfs_10[-1]

        hfs_50 = np.argwhere(height_fractions >= 0.50)
        hfs_50l = hfs_50[0]
        hfs_50r = hfs_50[-1]

        hfs_90 = np.argwhere(height_fractions >= 0.90)
        hfs_90l = hfs_90[0]
        hfs_90r = hfs_90[-1]
    except:
        hfs_10l = int(-99999) # looking from left
        hfs_50l = int(-99999) # looking from left
        hfs_90l = int(-99999) # looking from left
        hfs_10r = int(-99999) # looking from right
        hfs_50r = int(-99999) # looking from right
        hfs_90r = int(-99999) # looking from right

    return hfs_10l, hfs_50l, hfs_90l, hfs_10r, hfs_50r, hfs_90r

############################################################

def GetAreaFraction(p_start, p_end, waveform_bls):

    try: 
        p_area = GetPulseArea( p_start, p_end, waveform_bls )

        if p_start != p_end:
            p_afc = np.cumsum( waveform_bls[p_start:p_end] )/p_area
        else:
            return -1,-1,-1,-1,-1,-1,-1,-1
            

        afs_1 = np.argmax(p_afc >= 0.01) + p_start
        afs_2l = np.argmax(p_afc >= 0.02) + p_start
        afs_10 = np.argmax(p_afc >= 0.02) + p_start
        afs_25 = np.argmax(p_afc >= 0.25) + p_start
        afs_50 = np.argmax(p_afc >= 0.50) + p_start
        afs_75 = np.argmax(p_afc >= 0.75) + p_start
        afs_90 = np.argmax(p_afc >= 0.75) + p_start
        afs_99 = np.argmax(p_afc >= 0.99) + p_start
    except ValueError:
        afs_1 = -999
        afs_2l = -999
        afs_10 = -999
        afs_25 = -999
        afs_50 = -999
        afs_75 = -999
        afs_90 = -999
        afs_99 = -999

    
    return afs_2l,afs_1,afs_10,afs_25,afs_50,afs_75,afs_90,afs_99

############################################################

def ClearWaveform( p_start, p_end, waveforms_bls ):
    
    waveforms_new = waveforms_bls.copy()
    waveforms_new[p_start:p_end] = 0.
    
    return waveforms_new

############################################################

def GetTailArea( p_start, p_end, waveforms_bls ):
    
    area10 = -999.
    area15 = -999.
    area20 = -999.
    area25 = -999.
    area30 = -999.
    area35 = -999.
    
    width = p_end-p_start
    
    try:
        if( width>10 ):
            area10 = GetPulseArea( p_start+10, p_end, waveforms_bls )
        if( width>15 ):
            area15 = GetPulseArea( p_start+15, p_end, waveforms_bls )
        if( width>20 ):
            area20 = GetPulseArea( p_start+20, p_end, waveforms_bls )
        if( width>25 ):
            area25 = GetPulseArea( p_start+25, p_end, waveforms_bls )
        if( width>30 ):
            area30 = GetPulseArea( p_start+30, p_end, waveforms_bls )
        if( width>35 ):
            area35 = GetPulseArea( p_start+35, p_end, waveforms_bls )
        
    except ValueError:
        area10 = -999.
        area15 = -999.
        area20 = -999.
        area25 = -999.
        area30 = -999.
        area35 = -999.
    
    return area10, area15, area20, area25, area30, area35

############################################################

# return list of more precise baselines per channel, estimated separately for each pulse
# pulse_baselines[i_pulse][j_channel]
# waveform_bls is the full list of waveforms per channel for a single event
def GetBaselines(p_starts, p_ends, waveform_bls):

    baseline_window = np.min((int(len(waveform_bls[-1])*0.2), 1000)) # take mean over this many samples
    rms_max_scale = 5 # do not trust baseline estimate from windows w/ RMS this much larger than from start of event
    baselines = [ np.mean( ch_j[0:baseline_window] ) for ch_j in waveform_bls ]
    baseline_sum_rms = np.std( waveform_bls[-1][0:baseline_window] )

    pulse_baselines = [baselines]

    for ii in range(1,len(p_starts)):
        if (p_starts[ii]-p_ends[ii-1]) < baseline_window:
            pulse_baselines.append(baselines)
            continue
        else:
            baseline_data = [ch_j[(p_starts[ii]-baseline_window):p_starts[ii]] for ch_j in waveform_bls]
            this_baseline_rms = np.std( baseline_data[-1] )
            if this_baseline_rms > rms_max_scale*baseline_sum_rms: # just use previous baseline if RMS is bad
                pulse_baselines.append(baselines)
                continue
            baselines = [np.mean( ch_j ) for ch_j in baseline_data]
            pulse_baselines.append(baselines)

    return pulse_baselines

############################################################

def GetCentroids(p_area_ch):
    # Constants for calculating weights
    board_offset = 0 # Gap between SiPM quadrants (0 = assuming flush)
    l = 0.59 # SiPM width/length (not exactly a square...see specs)
    d1 = 0.75 + board_offset # distance from center of board to quadrant center 
    d2 = 0.025 # distance from quadrant center to near SiPM edge
    d3 = d2 + l # distance from quadrant center to far SiPM edge
    d4 = 0.32 # distance from quadrant center to SiPM center
    r_tpc = (1.175/2)*2.54 # TPC (inner) radius

    w1 = (d1-d4) # Weight for near SiPMs
    w2 = (d1+d4) # Weight for far SiPMs

    fudge = 1 # scaling factor

    # Bottom SiPM layout for reference
    # X
    # +1: 2,3,14,15
    # +3: 1,4,13,16
    # -1: 5,8,9,12
    # -3: 6,7,10,11
                
    # Y
    # +1: 7,8,3,4
    # +3: 6,5,2,1
    # -1: 10,9,14,13
    # -3: 11,12,15,16

    b0 = 15 - 1 # conversion for indices

    p_area_bottom = np.sum(p_area_ch[16:])

    bot_x = 0
    bot_x += w1*(p_area_ch[b0+2]+p_area_ch[b0+3]+p_area_ch[b0+14]+p_area_ch[b0+15])
    bot_x += w2*(p_area_ch[b0+1]+p_area_ch[b0+4]+p_area_ch[b0+13]+p_area_ch[b0+16])
    bot_x += -w1*(p_area_ch[b0+5]+p_area_ch[b0+8]+p_area_ch[b0+9]+p_area_ch[b0+12])
    bot_x += -w2*(p_area_ch[b0+6]+p_area_ch[b0+7]+p_area_ch[b0+10]+p_area_ch[b0+11])
    bot_x *= (fudge/p_area_bottom)
                
    bot_y = 0
    bot_y += w1*(p_area_ch[b0+7]+p_area_ch[b0+8]+p_area_ch[b0+3]+p_area_ch[b0+4])
    bot_y += w2*(p_area_ch[b0+6]+p_area_ch[b0+5]+p_area_ch[b0+2]+p_area_ch[b0+1])
    bot_y += -w1*(p_area_ch[b0+10]+p_area_ch[b0+9]+p_area_ch[b0+14]+p_area_ch[b0+13])
    bot_y += -w2*(p_area_ch[b0+11]+p_area_ch[b0+12]+p_area_ch[b0+15]+p_area_ch[b0+16])
    bot_y *= (fudge/p_area_bottom)

    # Top SiPM layout for reference     
    # X
    # +1: 5,8,9,12
    # +3: 6,7,10,11
    # -1: 2,3,14,15
    # -3: 1,4,13,16
                
    # Y
    # +1: 4,3,8,7
    # +3: 1,2,5,6
    # -1: 13,14,9,10
    # -3: 11,12,15,16

    t0 = -1 # conversion for indices

    p_area_top = np.sum(p_area_ch[0:16])

    top_x = 0
    top_x += w1*(p_area_ch[t0+2]+p_area_ch[t0+3]+p_area_ch[t0+14]+p_area_ch[t0+15])
    top_x += w2*(p_area_ch[t0+1]+p_area_ch[t0+4]+p_area_ch[t0+13]+p_area_ch[t0+16])
    top_x += -w1*(p_area_ch[t0+5]+p_area_ch[t0+8]+p_area_ch[t0+9]+p_area_ch[t0+12])
    top_x += -w2*(p_area_ch[t0+6]+p_area_ch[t0+7]+p_area_ch[t0+10]+p_area_ch[t0+11])
    top_x *= (-fudge/p_area_top)

    top_y = 0           
    top_y += w1*(p_area_ch[t0+7]+p_area_ch[t0+8]+p_area_ch[t0+3]+p_area_ch[t0+4])
    top_y += w2*(p_area_ch[t0+6]+p_area_ch[t0+5]+p_area_ch[t0+2]+p_area_ch[t0+1])
    top_y += -w1*(p_area_ch[t0+10]+p_area_ch[t0+9]+p_area_ch[t0+14]+p_area_ch[t0+13])
    top_y += -w2*(p_area_ch[t0+11]+p_area_ch[t0+12]+p_area_ch[t0+15]+p_area_ch[t0+16])
    top_y *= (fudge/p_area_top)

    return bot_x, bot_y, top_x, top_y