import numpy as np
import datetime as dt

import os

import matplotlib.pyplot as plt

from sklearn.cluster import DBSCAN
from sklearn import metrics
from sklearn.datasets.samples_generator import make_blobs
from sklearn.preprocessing import StandardScaler


import pandas as pd

class MicroSeisEvent:
    """docstring for MicroSeisEvent"""
    def __init__(self, evtime, evn, eve, evdp, mag=None, srcrad=None):
        self.evtime = evtime
        self.evn = evn
        self.eve = eve
        self.evdp = evdp
        self.mag = mag
        self.srcrad = srcrad
        
class Station:
    """docstring for Station"""
    def __init__(self, network, name, n, e, dp, cmp1vect, cmp2vect,cmp3vect):
        self.network = network
        self.name = name
        self.n = n
        self.e = e
        self.dp = dp
        self.cmp1vect = cmp1vect
        self.cmp2vect = cmp2vect
        self.cmp3vect = cmp3vect
        
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def randrange(n, vmin, vmax):
    return (vmax - vmin)*np.random.rand(n) + vmin
    
def findcluster(events,maxdist):
    cluster=[ev for ev in events if np.sqrt((ev.eve-largeev.eve)**2+(ev.evn-largeev.evn)**2+(ev.evdp-largeev.evdp)**2)<maxdist]
    
    return cluster


if __name__ == '__main__':		
    HeaderDir = os.path.expanduser('../Data')

    # Get event info		
    fname=os.path.join(HeaderDir,'D-37-H_Event_Catalogue.txt')
    colnames=('date,time,rel_eve,rel_evn,rel_evdp,eve,evn,evdp,mag,srcrad')
    ev = np.genfromtxt(fname,names=colnames,dtype=None)

    events=[]
    for event in ev:
        time=dt.datetime.strptime(event['date']+' '+event['time'],'%d/%m/%y %H:%M:%S')
        tmpevent=MicroSeisEvent(time,event['evn'],event['eve'],event['evdp'],mag=event['mag'],srcrad=event['srcrad'])
        events.append(tmpevent)

    # Get station info    
    fname=os.path.join(HeaderDir,'Sensor-Information.txt')
    colnames=('stname,stn,ste,stdp,sens,tmp1,tmp2,flo,fhi,gain,sensitivity,orient_N,orient_E,orient_D')
    ss = np.genfromtxt(fname,names=colnames,dtype=None)

    stations=[]
    for i in range(0,len(ss),3): # read 1st station then every third station
        if ss[i]['stname'] != ss[i+1]['stname'] or ss[i]['stname'] != ss[i+1]['stname']:
            print "stations dont match"
        else:
            network = ss[i]['stname'][0:-4] # well name but not last 4 characters of string (station number)
            name = ss[i]['stname'][-3:] # last 3 numbers of string (station number)
            e  = ss[i]['ste']
            n  = ss[i]['stn']
            dp = ss[i]['stdp']
            cmp1vect=np.array([-ss[i]['orient_D'],ss[i]['orient_N'],ss[i]['orient_E']])
            cmp2vect=np.array([-ss[i+1]['orient_D'],ss[i+1]['orient_N'],ss[i+1]['orient_E']])
            cmp3vect=np.array([-ss[i+2]['orient_D'],ss[i+2]['orient_N'],ss[i+2]['orient_E']])
            stations.append(Station(network, name, n, e, dp, cmp1vect, cmp2vect,cmp3vect))

    # Sort events by magnitude
    events.sort(key=lambda r: r.mag)
    
    # Choose an event at centre of cluster - reference event
    # Find largest event
    largeev=events[-1]
    # Find second largest event
    #secondev=events[-2]
    print largeev.mag

    # Find hypocentral distance between reference event and receivers
    HypoDist=[]
    for stat in stations:
        # Data used is K and S well only
        if stat.network != "M_Well":
            EpiDist = np.sqrt((largeev.evn-stat.n)**2. +  (largeev.eve-stat.e)**2.)
            HypoDist.append(np.sqrt((largeev.evdp-stat.dp)**2. + EpiDist**2))
    #print min(HypoDist), max(HypoDist)
    maxdist=max(HypoDist)/10.
    
    # Find the cluster of events within maxdist of the reference event
    cluster=findcluster(events,maxdist)
    
    cluster.sort(key=lambda r: r.mag)
    #largesubevts=cluster[-10:-1]
    
    
    # Plot the cluster
    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # ax.scatter([i.eve for i in cluster], [i.evn for i in cluster], [-i.evdp for i in cluster],c='b')
    # ax.set_xlabel('E')
    # ax.set_ylabel('N')
    # ax.set_zlabel('Depth')
    # plt.axis('equal')
    # plt.show()
    

    print len(cluster)
    