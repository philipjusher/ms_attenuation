import numpy as np

from obspy import read

from obspy.signal.polarization import flinn
from obspy.signal.trigger import zDetect, triggerOnset, plotTrigger

import datetime as dt


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
    
    
def picker(st):
    """docstring for picker"""
    
    s_picks=[np.nan]*2
    
    df = st[0].stats.sampling_rate
    npts = st[0].stats.npts
    
    for i,cmp in enumerate(['L','Q','T']):
        tmp=st.select(component=cmp)
        trace=tmp[0]
        # calculate cft of zDetect for 20 sample window
        cft=zDetect(trace,20)
        #mincft=min(cft)
        maxcft=max(cft)
        #meancft=mean(cft)
        #print mincft,maxcft,(maxcft+mincft)/2,meancft
    
    
        if cmp=='L': # P-trace: care less about this one, go ahead and pick
            picks=triggerOnset(cft,0,0)
            
            # plotTrigger(trace, cft, 0, 0)
        
            p_pick=np.nan
            for pick,end in picks:
                if npts-pick<100:
                    break
                else: 
                    stcut = st.copy()
                    # cut trace from pick to 100 samples past pick
                    startcut=stcut[0].stats.starttime+dt.timedelta(seconds= pick/df)
                    endcut=startcut + dt.timedelta(seconds=100/df)
                    stcut = stcut.slice(startcut,endcut)
                    stcut.sort()
                
                    azi,inc,lin,plan=flinn(stcut)
                    
                    #print "inc, lin =", inc, lin
                    if inc<40 and lin>0.6:
                        
                        # good enough for a p-pick
                        p_pick=pick
                        break
            
        
        
        else: # S-traces: only pick if the signal is relatively large
            s_picks[i-1]=np.nan
            if maxcft>6:
                picks=triggerOnset(cft,4,1)
            
                for pick,end in picks:
                    if npts-pick<100:
                        break
                    else:
                        stcut =st.copy()
                        # cut trace from pick to 100 samples past pick
                        startcut=stcut[0].stats.starttime+dt.timedelta(seconds= pick/df)
                        endcut=startcut + dt.timedelta(seconds=100/df)
                        stcut = stcut.slice(startcut,endcut)
                        stcut.sort()
                    
                        azi,inc,lin,plan=flinn(stcut)
                        if inc>50 and lin>0.75:
                            # good enough for an s-pick
                            s_picks[i-1]=pick
                            break

    s_pick=np.nanmin(s_picks)

    return (p_pick-50)/df,[(i-50)/df for i in s_picks]

