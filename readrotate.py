import numpy as np

from obspy import read

def calc_geoinc(trace,metric=True):
    """docstring for calc_geoinc"""
    stla=trace.stats.sac.stla
    stlo=trace.stats.sac.stlo
    stdp=trace.stats.sac.stdp
    evla=trace.stats.sac.evla
    evlo=trace.stats.sac.evlo
    evdp=trace.stats.sac.evdp
    
    
    if metric:
        baz=np.rad2deg(np.arctan2((evlo-stlo),(evla-stla)))
        EpiDist = np.sqrt((evlo-stlo)**2. +  (evla-stla)**2.)
        inc = np.rad2deg(np.arctan(EpiDist/ (evdp-stdp)))
        
        HypoDist = np.sqrt((evdp-stdp)**2. + EpiDist**2)
         
    if baz<0.:
        baz=baz+360.
        
    
    azi=np.mod(baz+180.0,360.0)
    inc=np.mod(inc+180.0,180.)

    
    return azi,inc,baz,HypoDist,EpiDist,stdp
    
    
def fixchannels(st):
    """docstring for fixchannels"""
    for trace in st:
        if len(trace.stats.channel)==1:
            trace.stats.channel='BS'+trace.stats.channel
    
    pass
    
def readandrotate(stem,stns='*',ext='[E,N,Z]'):
    """
    Reads in sac files and rotates to ray-frame
    stem:   string providing path to sac files without the extension
            e.g.: 'path/to/sac/sacfile'
    stns:   either '*' for all stations in array, or list of stations
            e.g.: '[001,002,003]'
    ext:    extensions corresponding to 3C data
            e.g.: '[BHE,BHN,BHZ]'
    """
    
    # Read in files
    st = read(stem+'.'+stns+'.'+ext)
    
    # fix the chgannel components if they are not 3 characters long (stupidly required for rotating)
    fixchannels(st)
    
    # copy stream
    strot = st.copy()
    for station in [tr.stats.station for tr in strot.select(component="Z")]:
        sttmp = strot.select(station=station)
        azi,inc,baz,hypodist,epidist,stdp=calc_geoinc(sttmp[0])
        sttmp.rotate('ZNE->LQT',back_azimuth=baz,inclination=inc)

    return strot