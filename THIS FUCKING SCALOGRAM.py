# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 18:31:23 2018

@author: ludov
"""

from electroPy import HdF5IO as H5

import numpy


# Class analogsignal = convert HdF5 called with HdF5IO in Analogsignal object  COMES FROM ANALOGSIGNAL CLASS OPENEPHY
class Do_AnalogSignal(object):
         
    def __init__(self, signal=None, channel=None, name=None,time_vector=None, 
        sampling_rate=1., t_start=0., t_stop=None, dt=None, **kargs):

        self.signal = signal
        self.channel = channel
        self.name = name
        self.sampling_rate = float(sampling_rate)
        self.t_start = float(t_start)
        self.t_stop = t_stop
        
        # Default for signal is empty array
        if self.signal is None:
            self.signal = numpy.array([])
        
        # Override sampling rate if dt is specified
        if dt is not None:            
            self.sampling_rate = 1. / dt        
        if self.sampling_rate == 0.:
            raise(ValueError("sampling rate cannot be zero"))
        
        # Calculate self.t_stop
        if self.t_stop is None:
            self.t_stop = self.t_start + len(self.signal)/self.sampling_rate

        # Initialize variable for time array, to be calculated later
        self._t = time_vector

    def compute_time_vector(self) :
        return numpy.arange(len(self.signal), dtype = 'f8')/self.sampling_rate + self.t_start

    def t(self):
        #if self._t==None:
        self._t=self.compute_time_vector()
        return self._t

    def max(self):
        return self.signal.max()

    def min(self):
        return self.signal.min()
    
    def mean(self):
        return numpy.mean(self.signal)
    
    def time_slice(self, t_start, t_stop):

        # Get time axis and also trigger recalculation if necessary
        t = self.t()
        
        # These kinds of checks aren't a good idea because of possible
        # floating point round-off error
        #assert t_start >= t[0], "not enough data on the beginning"
        #assert t_stop <= t[-1], "not enough data on the end"
        assert t_stop > t_start, "t_stop must be > t_start"
        
        # Find the integer indices of self.signal closest to requested limits
        # Do the checks here for 
        i_start = int(numpy.rint((t_start - self.t_start) * self.sampling_rate))
        if i_start < 0:
            print ("warning: you requested data before signal starts")
            i_start = 0
        
        # Add one so that it is inclusive of t_stop
        # In the most generous case, all of the data will be included and
        # i_stop == len(self.signal), which is okay because of fancy indexing
        i_stop = int(numpy.rint((t_stop - self.t_start) * self.sampling_rate)) + 1
        if i_stop > len(self.signal):
            print ("warning: you requested data after signal ended")
            i_stop = len(self.signal)
        
        # Slice the signal
        signal = self.signal[i_start:i_stop]
        
        # Create a new AnalogSignal with the specified data and the correct
        # underlying time-axis
        result = Do_AnalogSignal(signal=signal, sampling_rate=self.sampling_rate, 
            t_start=t[i_start])
        return result

#Calling the Data
path = 'D:/Trials/2018-07-31T16-35-19McsRecording_P3.h5'

data = H5(path) 


ch_0 = data.raw_record()[0]

t_start = data.raw_time()[0]

time_vector = data.raw_time

samp_rate = data.raw_sampling_rate()

anasig = Do_AnalogSignal(signal=ch_0,sampling_rate=samp_rate,t_start=t_start,time_vector=time_vector)

from matplotlib import pyplot

''' NOW THE SHITTY STUF CONVOLUTION--------------------------------------------'''


#COMING FROM CONVOLUTION SCRIPT
from scipy import *
from scipy.signal import resample
from scipy.fftpack import fft, ifft, fftshift

def generate_wavelet_fourier(len_wavelet,
            f_start,
            f_stop,
            deltafreq,
            sampling_rate,
            f0,
            normalisation,
            ):
    """
    Compute the wavelet coefficients at all scales and makes its Fourier transform.
    When different signal scalograms are computed with the exact same coefficients, 
        this function can be executed only once and its result passed directly to compute_morlet_scalogram
        
    Output:
        wf : Fourier transform of the wavelet coefficients (after weighting), Fourier frequencies are the first 
    """
    # compute final map scales
    scales = f0/arange(f_start,f_stop,deltafreq)*sampling_rate
    # compute wavelet coeffs at all scales
    xi=arange(-len_wavelet/2.,len_wavelet/2.)
    xsd = xi[:,newaxis] / scales
    wavelet_coefs=exp(complex(1j)*2.*pi*f0*xsd)*exp(-power(xsd,2)/2.)

    weighting_function = lambda x: x**(-(1.0+normalisation))
    wavelet_coefs = wavelet_coefs*weighting_function(scales[newaxis,:])

    # Transform the wavelet into the Fourier domain
    #~ wf=fft(wavelet_coefs.conj(),axis=0) <- FALSE
    wf=fft(wavelet_coefs,axis=0)
    wf=wf.conj() # at this point there was a mistake in the original script
    
    return wf

def compute_morlet_scalogram(ana,
            f_start=5.,
            f_stop=100.,
            deltafreq = 1.,
            sampling_rate = 200.,
            t_start = -inf, 
            t_stop = inf,
            f0=2.5, 
            normalisation = 0.,
            wf=None
            ):

    # Reduce signal to time limits
    sig=ana.signal[(ana.t()>=t_start)&(ana.t()<=t_stop)]
    
    if sig.size>0:
        if wf is None:
            if ana.sampling_rate != sampling_rate:
                sig=resample(sig,sig.size*sampling_rate/ana.sampling_rate)
            wf = generate_wavelet_fourier(sig.size,max(f_start,deltafreq),min(f_stop,ana.sampling_rate/2.),deltafreq,sampling_rate,f0,normalisation)
        else:
            if sig.size != wf.shape[0]:
                sig=resample(sig,wf.shape[0])

        # Transform the signal into the Fourier domain
        sigf=fft(sig)

        # Convolve (mult. in Fourier space)
        #~ wt_tmp=ifft(wf*sigf[newaxis,:],axis=1)
        wt_tmp=ifft(sigf[:,newaxis]*wf,axis=0)
 
        # shift output from ifft
        wt = fftshift(wt_tmp,axes=[0])
        
    else:
        scales = f0/arange(f_start,f_stop,deltafreq)*sampling_rate
        wt = empty((0,scales.size),dtype='complex')

    return wt
    
'''-------------------------------------------------------------------------'''
#COMING FROM TIMEFREQ SCRIPT
from numpy import inf

# global for caching wf
global cache_for_wf
cache_for_wf = None
global signature_for_wf
signature_for_wf = ''


class TimeFreq():
    doc2="""
    *TimeFreq*
    
    """
    docparam = """
    
    Params:
     :f_start: lkjlkj
    
    """
    
    __doc__ = doc2+docparam
    
    def __init__(self,
                        anaSig,
                        method = 'convolution_freq',
                        f_start=5.,
                        f_stop=20.,
                        deltafreq = 0.5,
                        sampling_rate = 10000.,
                        t_start = -inf, 
                        t_stop = inf,
                        f0=2.5, 
                        normalisation = 0.,
                        **kargs
                        ):
                        
        self.anaSig = anaSig
        self.method = method
        self.f_start=f_start
        self.f_stop=f_stop
        self.deltafreq = deltafreq
        self.sampling_rate = sampling_rate
        self.t_start = t_start
        self.t_stop = t_stop
        self.f0=f0
        self.normalisation = normalisation
        
        self.t_start = max(self.t_start , self.anaSig.t_start)
        self.t_stop = min(self.t_stop , self.anaSig.t()[-1]+1./self.anaSig.sampling_rate )
        
        self._map = None
        self._t = None
        self._f = None
        
        if self.method == 'convolution_freq':
            self._wf = None
            self.subAnaSig = None


    def compute_time_vector(self) :
        return numpy.arange(len(self.subAnaSig.signal), dtype = 'f8')/self.sampling_rate + self.t_start

    def compute_freq_vector(self) :
        return numpy.arange(self.f_start,self.f_stop,self.deltafreq, dtype = 'f8')

    def t(self):
        if self._t==None:
            self._t=self.compute_time_vector()
        return self._t
    
    def f(self):
        if self._f==None:
            self._f=self.compute_freq_vector()
        return self._f
    
        
    @property
    def map(self):
        if self._map is None:
            self.recomputeMap()
        return self._map


    def recomputeMap(self):
        """
        Compute or recompute a map
        """
        if self.subAnaSig is None:
        #~ if True:
            sig=self.anaSig.signal[(self.anaSig.t()>=self.t_start)&(self.anaSig.t()<self.t_stop)]
            if self.sampling_rate != self.anaSig.sampling_rate :
                sig=resample(sig,sig.size*self.sampling_rate/self.anaSig.sampling_rate)
            self.subAnaSig = Do_AnalogSignal( signal = sig,
                                                            sampling_rate = self.sampling_rate,
                                                            t_start = self.t_start,
                                                        )
        
        if self.method == 'convolution_freq':
            if self._wf is None :
                global signature_for_wf
                global cache_for_wf
                signature = '%d %f %f %f %f %f %f' % (self.subAnaSig.signal.size,
                                                                                self.f_start,
                                                                                self.f_stop,
                                                                                self.deltafreq,
                                                                                self.sampling_rate,
                                                                                self.f0,
                                                                                self.normalisation,
                                                                                )
                if signature != signature_for_wf:
                    cache_for_wf= generate_wavelet_fourier(len_wavelet=self.subAnaSig.signal.size,
                                                                                f_start=self.f_start,
                                                                                f_stop=self.f_stop,
                                                                                deltafreq=self.deltafreq,
                                                                                sampling_rate=self.sampling_rate,
                                                                                f0=self.f0,
                                                                                normalisation = self.normalisation,
                                                                                )
                    signature_for_wf = signature
                
                self._wf = cache_for_wf
            
            
            self._map = compute_morlet_scalogram(self.subAnaSig,wf = self._wf )
            
        
    def plotMap(self, ax,
                                    colorbar = True,
                                    cax =None,
                                    orientation='horizontal',
                                    **kargs):
        """
        
        ax : a matplotlib axes
        
        """
        im = ax.imshow(abs(self.map).transpose(),
                                    interpolation='nearest', 
                                    extent=(self.t_start, self.t_stop, self.f_start-self.deltafreq/2., self.f_stop-self.deltafreq/2.),
                                    #origin ='lower' ,
                                    )

        if colorbar:
            if cax is None:
                ax.figure.colorbar(im)
            else:
                ax.figure.colorbar(im,ax = ax, cax = cax ,orientation=orientation)
            
                
        return im
        



#DOING THE SCALO + FIGURE

fig=pyplot.figure()
ax=fig.add_subplot(1,1,1)
THE_TRIAL = TimeFreq(anasig, sampling_rate=samp_rate)
ax.set_title('%s'%path)


THE_TRIAL.plotMap(ax)






