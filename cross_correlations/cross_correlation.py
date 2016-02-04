import numpy as np
from obspy.core import read



def cross_correlate(trace1,trace2):
    """Returns the correlation coefficient from two obspy traces"""
    return np.corrcoef(trace1.data,trace2.data)
    

def window_trace(st,before_p,after_p):
    """Cuts a window around an obspy trace using a before and after p pick time (sec)"""
    for tr in st:
        p_pick=tr.stats.starttime+tr.stats.sac['t0']
        tr.trim(p_pick-before_p,p_pick+after_p)
    return st

def get_p(st):
    """Get the P-wave trace"""
    p_traces = st.select(component="L")
    
    return p_traces
    
