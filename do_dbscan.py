from sklearn.cluster import DBSCAN
from sklearn import metrics
import numpy as np
import datetime as dt
import os

import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
        
def plot_clusters(X,db):
    """
    Function to plot clusters found from DBSCAN
    """
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    # Black removed and is used for noise instead.
    unique_labels = set(db.labels_)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
    fig=plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    for k, col in zip(unique_labels, colors):
        #k=-1 are not in a cluster
        #if k == -1:
        #    # Black used for noise.
        #    col = 'k'        
        if k!=-1:
            class_member_mask = (labels == k)
            xy = X[class_member_mask & core_samples_mask]            
            ax.scatter(xy[:, 0], xy[:, 1], xy[:, 2], c=col, marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.title('Estimated number of clusters: %d' % n_clusters_)
    plt.show()
    
        
        
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
        #print tmpevent.evtime
    
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
    
    # Create synthetic data for testing
    #X = np.random.rand(500,3)

    # Array of event locations
    #print np.array([[i.evtime] for i in events])
    
    X = np.array([[i.evtime,i.eve,i.evn,i.evdp] for i in events])
    #print X[:,1:4]

    df=pd.DataFrame(X)
    
    
    # Start with low min_samples (5) to ensure clusters & start with eps=0.1 & increase
    # eps is distance parameter - sets max distance between clusters
    db = DBSCAN(eps=10.0, min_samples=25).fit(X[:,1:4])
    
    print db.labels_[0:20]
    
    #Labels tell you which cluster event is in
    labels = db.labels_
    
    from collections import Counter
    # Cluster labelled -1 are all events not in a cluster
    # Counts number of events in each cluster
    #print Counter(labels),'\n'
    #print df[labels==-1]

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    print('Estimated number of clusters: %d' % n_clusters_)
    
    # Find list of indices for events in each cluster 
    for n in range(n_clusters_):
        list2 = [ind for ind, x in enumerate(labels) if x==n]
    
    print list2
    
    for l in list2:
        ename=events[l].evtime.strftime('%Y%m%d_%H%M%S')
        print ename, events[l].eve, events[l].evn, events[l]. evdp,events[l].mag

        
    
    
    # Plot clusters if wanted
    #plot_clusters(X[:,1:4],db)
    '''
    for label in labels:
        if label==1:
            print np.zeros_like(label, dtype=bool)
        #samples = np.zeros_like(label, dtype=bool)
        #samples[db.core_sample_indices_]=True
    '''
    