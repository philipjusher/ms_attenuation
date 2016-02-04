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


def calc_mtspec(obspy_trace):
    """A function to calculate a spectral multi-taper Fourier transform"""
    
    return freq,amp
    
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
    return grad,grad_err


            
   
def bs_lsr(trace_spectra,ref_spectra,freq,fmin,fmax,noise_trace,noise_ref,snr_limit=1):
    """Returns delta t star between the two traces (trace and ref)
    
    This is an implemntation of the log spectral ratio method to measure delta
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
    snr=trace_spectra/noise_trace
    snr_ref=ref_spectra/noise_ref
    lsr=np.log(trace_spectra/ref_spectra)
    y_weights=(snr+snr_ref)/2 
    
    #Indices at which to do the line fitting if there is not enough points return nan
    ind=(freq>fmin) & (freq<fmax) & (snr>snr_limit) & (snr_ref>snr_limit) & np.isfinite(lsr) 
    if sum(ind)>6:
        
        grad,grad_err=linearfit(freq[ind],lsr[ind],y_weights[ind])              
        meas_tstar=grad/-np.pi
        meas_tstar_error=grad_err/-np.pi
    else:
        meas_tstar=np.nan
        meas_tstar_error=np.nan
    
    return meas_tstar,meas_tstar_error
    
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