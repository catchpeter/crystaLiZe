import numpy as np

def ClassifyPulses(tba, t_rise, n_pulses, p_area):
    max_pulses = np.size(tba)
    classification = np.zeros(max_pulses)

    # May want everything below rise time of 0.2 to be S1-like, can subdivide S1 -> gas vs liquid?
    # Tricky bc discrimination space for gas-like pulses looks pretty different when voltages are on vs off
    # Clearly want to keep green as normal S1, always below ~0.125
    #max_t_rise = (1.5/100/np.log10(2))*(np.sign(p_area)*(np.absolute(p_area)**(np.log10(2)/2.5))) #8 sipms setup
    

    # log_rise = np.log10(t_rise)
    # log_area = np.log10(p_area)
    # log_rise = np.nan_to_num(log_rise, nan = -0.7, posinf = -0.7, neginf = -0.7)
    # log_area = np.nan_to_num(log_area, nan = 3.5, posinf = 3.5, neginf = 3.5)

    # s2_lower_limit = np.minimum(1.6/1.7*(log_area-2.2), -0.92)
    # s1_upper_limit = s2_lower_limit
    #s2_lower_limit = -0.225*(log_area-4)**2-0.3
    #s1_upper_limit = s2_lower_limit
    #s1_upper_limit[(log_area>2.5)*(log_area<4.1)] = -0.8
    riseTimes1s2 = -0.57
    case1 = (tba < 0)*(t_rise < riseTimes1s2) # normal S1s
    case2 = (tba >= 0)*(t_rise < riseTimes1s2) # top-focused S1s; e.g. in gas or LXe above top array
    case3 = (tba > -0.25)*(t_rise >= riseTimes1s2) # normal-ish S2s
    case4 = (tba <= -0.25)*(t_rise >= riseTimes1s2) # Unclear; possible S1/S2 merged pulses

    classification[case1] = 1
    classification[case2] = 2
    classification[case3] = 3
    classification[case4] = 4
    #classification[case5] = 5
    #classification[case6] = 6
    classification[n_pulses:] = 0 # Set pulse class to 0 for empty pulses

    
    return classification

