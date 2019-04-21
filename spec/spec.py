import numpy as np
import astropy.units as u
import astropy.constants as const
from scipy import interpolate

__all__ = ['Spectrum']

class Spectrum():
    '''
    Object for implement radio spectrum 

    Parameters
    ----------
    x : '~astropy.Quantity' with frequency or velocity unit
    y : '~numpy.ndarray'
    
    kwargs
    ------
    restfreq : '~astropy.Quantity' with frequency unit
        
    '''

    def __init__(self, x, y, **kwargs):
        if not isinstance(x, u.Quantity):
            raise TypeError('freq must be an astropy.units.Quantity instance')

        if not isinstance(y, u.Quantity):
            if not isinstance(y, np.ndarray):
                raise TypeError('value must be an astropy.units.Quantity or numpy.ndarray instance')
            else:
                y = y * u.one

        if x.shape != y.shape:
            raise ValueError('shape of freq and value are not identical {0} {1}'.format(freq.shape, value.shape))
        
        self.x = x 
        self.y = y 
        
        if x.unit in [u.Hz, u.kHz, u.MHz, u.GHz]:
            self.freq = x
            if 'restfreq' in kwargs.keys():
                self.velo = (kwargs['restfreq'] - self.freq) / kwargs['restfreq'] * const.c
        elif x.unit in [u.m/u.s, u.km/u.s]:
            self.velo = x
            if 'restfreq' in kwargs.keys():
                self.freq = (const.c - self.velo) / const.c * kwargs['restfreq']
        else:
            raise ValueError('x must be in frequency or velocity unit')



    def set_restfreq(self, restfreq):
        '''
        set rest frequency for self instance
        '''
        self.restfreq = restfreq
        if 'freq' in self.__dict__:
            self.velo = (1 - self.freq / restfreq) * const.c
        if 'velo' in self.__dict__:
            self.freq = (1 - self.velo / const.c) * restfreq


    def lsr_corr(self, vlsr):
        '''
        Vlsr correction

        Paramters
        ---------
        vlsr : '~astropy.Quantity'
            the velocity of LSR in respect of observatory in direction of
            the source
        
        Returns
        -------
        spect_shifted : Spectrum
            A Spectrum with shifted x axes
        '''
        if 'restfreq' not in self.__dict__:
            raise ValueError('No rest frequency')

        spec_shifted = self

        spec_shifted.velo = self.velo - vlsr
        spec_shifted.freq = (const.c - spec_shifted.velo) / const.c * self.freq

        return spec_shifted

    
    def regrid(self, x1, kind='linear'):
    '''
    Parameters
    ----------
    x1 : '~astropy.Quantity'
        The frequency or velocity range that the spectrum is to be regridded
        to
    kind : string
        kind of scipy.interpolate
        Specifies the kind of interpolation as a string (‘linear’, ‘nearest’, 
        ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’, ‘next’, where 
        ‘zero’, ‘slinear’, ‘quadratic’ and ‘cubic’ refer to a spline 
        interpolation of zeroth, first, second or third order; ‘previous’ and 
        ‘next’ simply return the previous or next value of the point) or as 
        an integer specifying the order of the spline interpolator to use. 
        Default is ‘linear’.

    Returns
    -------
    spec1 : Spectrum
        Regridded spectrum 
    '''
    if x1.unit in [u.Hz, u.kHz, u.MHz, u.GHz]:
        if 'freq' not in self.__dict__.keys():
            raise ValueError('No frequency in {0}'.format(self.__name__))
        if x1.min() < self.freq.min() or x1.max() > self.freq.max():
            raise ValueError('The new frequency is outside the original')
        x0 = self.freq
    elif x1.unit in [u.m/u.s, u.km/u.s]:
        if 'velo' not in self.__dict__.keys():
            raise ValueError('No velocity in {0}'.format(self.__name__))
        if x1.min() < self.velo.min() or x1.max() > self.velo.max():
            raise ValueError('The new velocity is outside the original')
        x0 = self.velo
    else:
        raise ValueError('new x must be in frequency or velocity unit')

    x1 = x1.to(x0.unit)
    y0 = self.y
    f = interpolate.interp1d(x0, y0, kind=kind)
    y1 = f(x1) 
    spec1 = Spectrum(x1, y1)
    if 'restfreq' in self.__dict__.keys():
        spec1.set_restfreq(self.restfreq)

    return spec1
    
    #def smooth():
    #def interpolate():


