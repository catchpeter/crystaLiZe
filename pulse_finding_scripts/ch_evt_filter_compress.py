"""
Signal-processing functions for crystalize data processing.
"""
import numpy as np
import c_process as cpr

def filter_channel_event(d, n1=5, n2=21):
    """
    Takes the raw waveform from one channel, from one event, and applies a filter
    which can be used for further processing.  The filter applies a box filter
    twice in series.
    Inputs:
         d: The input waveform. 1D numpy array of dtype float or np.float64. Should be 
            already baseline subtracted, but this is not a requirement.
        n1: The size of the first box filter, in number of samples.  Must be odd
        n2: The size of the second box filter, in number of samples.  Must be odd
    Outputs:
        dfilt: The filtered version of d
    """
    if not (isinstance(d, np.ndarray) and (d.dtype == np.float64)):
        raise TypeError("Input 'd' must be a numpy array of dtype np.float64")
    if d.ndim > 1:
        raise ValueError("Input 'd' must be a 1D array")
    if not isinstance(n1, int):
        raise TypeError("Input 'n1' must be of type int")
    if not isinstance(n2, int):
        raise TypeError("Input 'n2' must be of type int")
    if (n1%2) != 1:
        raise ValueError("n1 must be odd")
    if (n2%2) != 1:
        raise ValueError("n2 must be odd")
    dfilt = cpr.avebox(d, n1)
    dfilt = cpr.avebox(dfilt, n2)
    return dfilt

def baseline_suppress(d, pls_thresh=5., buffL=100, buffR=100, condense_thresh=100):
    """
    Takes a waveform and suppresses the baseline.
    Inputs:
                   d: The input waveform to baseline-suppress
        pulse_thresh: The height threshold that will be applied to look for regions
                      of data to keep.
               buffL: The additional number of samples to the left (i.e. before) to keep
               buffR: The additional number of samples to the right (i.e. after) to keep
     condense_thresh: If two intervals are closer than this number of samples, it will
                      condense them into a single interval to keep.
    Outputs:
              bs_cut: A numpy array of dtype bool, the same size as input 'd'.  The
                      convention is: True means to keep the sample, False means to quash
    """
    if not (isinstance(d, np.ndarray) and (d.dtype == np.float64)):
        raise TypeError("Input 'd' must be a numpy array of dtype np.float64")
    if d.ndim > 1:
        raise ValueError("Input 'd' must be a 1D array")
    if isinstance(pls_thresh, int):
        pls_thresh = float(pls_thresh)
    if not isinstance(pls_thresh, float):
        raise TypeError("Input 'thresh' must be a float scalar")
    if not (isinstance(buffL, int) and isinstance(buffR, int)):
        raise TypeError("'buffL' and 'buffR' must be ints")
    if buffL<0:
        buffL = 0
    if buffR<0:
        buffR = 0
    if not isinstance(condense_thresh, int):
        raise TypeError("'condense_thresh' must be an int")
    
    pods_10n1 = np.concatenate((np.r_[0].astype(np.int8),
        np.diff(((d>pls_thresh)|(d<-pls_thresh)).astype(np.int8))))
    pod_starts, = np.nonzero(pods_10n1 == 1)
    pod_stops,  = np.nonzero(pods_10n1 == -1)
    if len(pod_starts) > len(pod_stops):
        pod_starts = pod_starts[:len(pod_stops)]
    if len(pod_stops) > len(pod_starts):
        pod_stops = pod_stops[:len(pod_starts)]
    pod_starts -= buffL
    pod_stops  += buffR
    pulse_intvls = cpr.condense_intervals(
        np.vstack([pod_starts, pod_stops]), condense_thresh)
    pulse_intvls[pulse_intvls<0] = 0
    pulse_intvls[pulse_intvls>=len(d)] = len(d)-1
    p_cut = cpr.pulse_bool(d,pulse_intvls[0,:], pulse_intvls[1,:])
    return p_cut

def main():
    import matplotlib.pyplot as plt
    d = np.random.randn(int(1e4))
    dfilt = filter_channel_event(d)
    p_cut = baseline_suppress(dfilt,pls_thresh=0.55,buffL=200,buffR=200)
    plt.figure(1,figsize=(14,4.7))
    plt.clf()
    plt.plot(d,'-',lw=.5,color=np.r_[1,1,1]*.4)
    t = np.r_[:len(d)]
    d[~p_cut] = np.nan
    plt.plot(t,d,'-',lw=.5,color=np.r_[.3,.3,.6])
    plt.xlim([0,1e4])
    plt.plot(dfilt,'r-',lw=.75)
    plt.plot(np.r_[plt.xlim(),np.nan,plt.xlim()],
        np.r_[-1,-1,np.nan,1,1]*0.55,'c--',lw=.75)
    plt.xlabel('Time [samples]')
    plt.ylabel('Amplitude [a.u.]')
    plt.minorticks_on()
    plt.title('See stdout for explanation')
    print('Gray is the raw trace', flush=True)
    print('Red is the filtered trace', flush=True)
    print('Cyan is the POD threshold', flush=True)
    print('Blue is the non-baseline-suppressed raw waveform', flush=True)
    plt.show()

if __name__ == "__main__":
    main()






