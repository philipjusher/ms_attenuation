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
        
def plot_clusters(X,labels):
    """
    Function to plot clusters found from DBSCAN
    """
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
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
    #X = np.random.rand(500,3)

    X = np.array([[i.eve,i.evn,i.evdp] for i in events])
    print type(X)
    print np.shape(X)
    df=pd.DataFrame(X)
    # Start with low min_samples (5) to ensure clusters & start with eps=0.1 & increase
    # eps is distance parameter - sets max distance between clusters
    db = DBSCAN(eps=10.0, min_samples=5).fit(X)
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    
    from collections import Counter
    # Cluster labelled -1 are all events not in a cluster
    print Counter(labels),'\n'
    print df[labels==-1]

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

    print('Estimated number of clusters: %d' % n_clusters_)
    #print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
    #print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
    #print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
    #print("Adjusted Rand Index: %0.3f"
    #      % metrics.adjusted_rand_score(labels_true, labels))
    #print("Adjusted Mutual Information: %0.3f"
    #      % metrics.adjusted_mutual_info_score(labels_true, labels))
    #print("Silhouette Coefficient: %0.3f"
    #      % metrics.silhouette_score(X, labels))

    ##############################################################################
    
    # Plot clusters if wanted
    plot_clusters(X,labels)