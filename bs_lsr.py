# -*- coding: utf-8 -*-
"""
This is a module containing the log spectral ratio method approach to measure attenuation

The module contains the function bs_lsr which measures delta t * using the method from 
Kelly, C. M., Rietbrock, A., Faulkner, D. R., & Nadeau, R. M. (2012). Temporal Changes 
in Attenuation associated with the 2004 M6. 0 Parkfield Earthquake. Journal of
Geophysical Research-Solid Earth, 1–48.

And coded into python by Philip Usher (University of Bristol)
 
"""

import numpy as np
import statsmodels.api as sm
from mtspec import mtspec
from readrotate import readandrotate
from obspy.core import read,Stream
import matplotlib.pyplot as plt
def linearfit(x, y, y_weights):
    """Robust Linear Regression to fit of x and y a straight line and return the gradient with an error
    
    x and y are numpy arrays of the data to be fitted. y_weights is a numpy 
    array of the weights to use during the fitting
    """
    #Add a constant to fit a straight line
    x=sm.add_constant(x, prepend=False)
    #Bring in weights
    y=y*y_weights
    x=x*np.column_stack((y_weights,y_weights))
    
    #Generate robust model using bisquare (aka TukeyBiweights) robust weights       
    rlm_model = sm.RLM(y, x, M=sm.robust.norms.TukeyBiweight())
    
    #Fit model
    rlm_results = rlm_model.fit()
    #Get the 95% confidence interval to use as error
    conf=rlm_results.conf_int()
    grad=rlm_results.params[0]
    grad_err=grad-conf[0][0]

    return grad,grad_err,rlm_results


            
   
def bs_lsr(trace_spectra,ref_spectra,trace_noise_spectra,ref_noise_spectra):
    """Returns delta t star between the two traces (trace and ref)
    
    This is an implementation of the log spectral ratio method to measure delta
    t star.
    Arguments:
        
    trace -- a numpy array of the fft of a window of the trace seismogram
    ref -- a numpy array of the fft of a window of the reference seismogram
    noise_trace -- a numpy array of the fft of a noise window of the trace
    seismogram
    noise_ref -- a numpy array of the fft of a noise window of the reference
    seismogram
    freq -- a numpy array of the corresponding frequencies for the ffts above
    fmin,fmax -- minimum and maximum frequencies for the regression bandwidth
    snr_limit -- a limit for the snr values below this will be ignored.
    
    """
    
    #Calculate the snr, lsr and y weights for the trace and the reference
    snr=trace_spectra/trace_noise_spectra
    snr_ref=ref_spectra/ref_noise_spectra
    lsr=np.log(trace_spectra/ref_spectra)
    weights=(snr+snr_ref)/2 
    
    return lsr,weights
    
def calc_dts(lsr,weights,freq,fmin,fmax,snr_limit=1):
    #Indices at which to do the line fitting if there is not enough points return nan
    ind=(freq>fmin) & (freq<fmax) & (weights>snr_limit) & np.isfinite(lsr) 
    if sum(ind)>6:
        
        grad,grad_err,linefit=linearfit(freq[ind],lsr[ind],weights[ind])              
        meas_tstar=grad/-np.pi
        meas_tstar_error=grad_err/-np.pi
    else:
        meas_tstar=np.nan
        meas_tstar_error=np.nan
    
    return meas_tstar,meas_tstar_error,linefit

def cross_correlate(trace1,trace2):
    """Returns the correlation coefficient from two obspy traces"""
    corr_coef = np.corrcoef(trace1.data,trace2.data)
    xcorr = corr_coef[0,1]
    return xcorr
    
def window_trace(st,before_p,after_p):
    """Cuts a window around an obspy trace using a before and after p pick time (sec)"""
    p_st = Stream()
    for tr in st:
        #Removing traces without p picks
        if tr.stats.sac['t0'] != -12345.0:
            
            p_pick=tr.stats.starttime+tr.stats.sac['t0']
            p_tr = tr.slice(p_pick-before_p,p_pick+after_p)
            p_st.append(p_tr)
    
                
    return p_st

def noise_window_trace(st,window_len_time):
    """Get noise window for st"""
    noise_st = Stream()
    for tr in st:
        noise_tr = tr.slice(tr.stats.starttime,tr.stats.starttime+window_len_time)
        noise_st.append(noise_tr)
    return noise_st
    
def get_p(st):
    """Get the P-wave trace"""
    p_traces = st.select(component="L")
    
    return p_traces
    

def mean_corr(evnm1,evnm2):
    """Calculate the mean correlation between two events"""
    
    st = readandrotate(evnm1,stns='*',ext='[E,N,Z]')
    st2 = readandrotate(evnm2,stns='*',ext='[E,N,Z]')
 
def calc_mtspec(tr):
    """Do mtspec on a trace"""
    nfft = len(tr.data)*2
    spec_den,freq = mtspec(tr.data,tr.stats.delta,3.5,nfft=nfft) 
    amp = np.sqrt(spec_den)
    return amp,freq 
 
def get_and_check_window_len(st,st2):
    len_st=[]
    for tr in st:
        len_st.append(len(tr.data))
    for tr in st2:
        len_st.append(len(tr.data))
        
    assert all(x==len_st[0] for x in len_st)
    
    return len_st[0]
    
def get_fft(evnm,p_before, p_after,stns='001',ext='[E,N,Z]'):
    
    #Reading in the Events
    st = readandrotate(evnm,stns=stns,ext='[E,N,Z]')
    
    
    #Get the P components
    st = get_p(st)

    
    assert len(st) ==1

    
    # Slice of the P-wave windows
    st_p = window_trace(st,p_before,p_after)

    
    assert len(st_p) != 0


    window_len_time = p_before + p_after
    # Calculates the noise window
    # Note uses the beginning of the trace if the P wave is near the beginning it will overlap
    st_n = noise_window_trace(st,window_len_time)
    
    #Calulate power spectral density of P,Noise for event and reference
    fft_p,freq = calc_mtspec(st_p[0])
    fft_n,freq = calc_mtspec(st_n[0])
        
    return(fft_p,fft_n,freq,st,st_p,st_n)
    
def dts_p(evnm1,evnm_ref,fmin,fmax,snr_limit,station='001',p_before=0.1,p_after=0.2,diagnostic_plot=False):
    
    fft_p,fft_n,freq,st,st_p,st_n= get_fft(evnm1,p_before, p_after,stns=station,ext='[E,N,Z]')
    fft_ref_p,fft_ref_n,freq_ref,st_ref,st_ref_p,st_ref_n = get_fft(evnm_ref,p_before, p_after,stns=station,ext='[E,N,Z]')
    
    assert all(freq == freq_ref)
    length = len(fft_p)
    assert all(len(lst) == length for lst in [fft_ref_p,fft_n,fft_ref_n])
    
    
    
    lsr, weights = bs_lsr(fft_p,fft_ref_p,fft_n,fft_ref_n)
    
    meas_tstar,meas_tstar_error,linefit = calc_dts(lsr, weights,freq,fmin,fmax,snr_limit=snr_limit)
    
    if diagnostic_plot:
        plt.figure(figsize=(8,12))
        
        # --------------- #
        plt.subplot(5,2,1)
        plt.title('Event Whole Trace')
        plt.plot(st[0].times(),st[0].data)
        ppick = st[0].stats.sac.t0
        plt.axvspan(ppick-p_before,ppick+p_after,color='r',alpha=0.5)
        plt.xlabel('Time (sec)')
        plt.ylabel('Amplitude')
        # --------------- #
        plt.subplot(5,2,2)
        plt.title('Reference Whole Trace')
        plt.plot(st_ref[0].times(),st_ref[0].data,'g')
        ppick = st_ref[0].stats.sac.t0
        plt.axvspan(ppick-p_before,ppick+p_after,color='r',alpha=0.5)
        plt.xlabel('Time (sec)')
        plt.ylabel('Amplitude')
        
        # --------------- #
        plt.subplot(6,2,3)
        plt.title('Trace Window')
        plt.plot(st_p[0].times(),st_p[0].data,label='P-wave')
        plt.plot(st_n[0].times(),st_n[0].data,'r',label='Noise')
        plt.xlabel('Time (sec)')
        plt.ylabel('Amplitude')
        
        # --------------- #
        plt.subplot(6,2,4)
        plt.title('Reference Window')
        plt.plot(st_ref_p[0].times(),st_ref_p[0].data,'g',label='P-wave')
        plt.plot(st_ref_n[0].times(),st_ref_n[0].data,'r',label='Noise')
        plt.xlabel('Time (sec)') 
        plt.ylabel('Amplitude')
        
        # --------------- #
        plt.subplot(6,2,5)
        plt.title('Trace FFT')
        plt.semilogy(freq,fft_p)
        plt.semilogy(freq,fft_n,'r')
        plt.xlabel('Freq (Hz)')
        plt.ylabel('Amplitude')
        
        # --------------- #
        plt.subplot(6,2,6)
        plt.title('Reference FFT')
        plt.semilogy(freq,fft_ref_p,'g')
        plt.semilogy(freq,fft_ref_n,'r') 
        plt.xlabel('Freq (Hz)')
        plt.ylabel('Amplitude')
        
        # --------------- #
        plt.subplot(6,2,7)
        plt.title('Log-spectral-ratio')
        plt.plot(freq,np.log(fft_p/fft_ref_p))
        plt.axvspan(fmin,fmax,color='r',alpha=0.5)             
        plt.xlabel('Freq (Hz)')
        plt.ylabel('LSR')
        
        ind = (freq>fmin) & (freq<fmax)
        linefit_freq = freq[ind]
        line = linefit_freq*linefit.params[0]+linefit.params[1]
        
        plt.plot(linefit_freq,line,'-w')
        # --------------- #
        plt.subplot(6,2,8)
        plt.title('Signal/Noise')
        snr=fft_p/fft_n
        snr_ref=fft_ref_p/fft_ref_n
        plt.semilogy(freq,snr)
        plt.semilogy(freq,snr_ref,'g')
        plt.axhline(snr_limit,color='r',ls='--')
        plt.xlabel('Freq (Hz)')        
        plt.ylabel('SNR')
        
        # --------------- #
        
        # --------------- #
        plt.subplot(3,1,3)
        plt.title('LSR in frequency band')
        plt.plot(freq,np.log(fft_p/fft_ref_p))
        plt.axvspan(fmin,fmax,color='r',alpha=0.5)             
        plt.xlabel('Freq (Hz)')
        plt.ylabel('LSR')
        plt.plot(linefit_freq,line,'-w')
        extra = 0.1*(fmax-fmin)
        plt.xlim(fmin-extra,fmax+extra)
        #TODO autoscale y axis
        # --------------- #
        plt.tight_layout()
        
    return meas_tstar,meas_tstar_error,linefit
    
        
if __name__=="__main__":
    
    #Define a brune source spectra for two different attenuations
    wo=1000.0
    freq=np.linspace(1,500,33)
    fc=500.0
    x=650.0
    c=3000.0
    t=x/c

    Qref=100.0
    ref_spectra=(wo*np.exp((-np.pi*freq*t)/Qref))/((1+((freq/fc)**2.0)))
    
    Q=80.0
    trace_spectra=(wo*np.exp((-np.pi*freq*t)/Q))/((1+((freq/fc)**2.0)))
    
    dts_theoretical=(t/Q)-(t/Qref)
    
    fmin=200
    fmax=400
    snr_limit=1
    noise_trace=np.ones(np.size(trace_spectra))
    noise_ref=np.ones(np.size(ref_spectra))
    delta_tstar,delta_tstar_err=bs_lsr(trace_spectra,ref_spectra,freq,fmin,fmax,noise_trace,noise_ref,snr_limit)
    
    print 'This should be very small', delta_tstar-dts_theoretical