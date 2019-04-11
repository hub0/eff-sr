import numpy as np
import astropy.units as u
import astropy.constants as const

__all__ = ['Spectrum']

class Spectrum():
    '''
    Object for implement radio spectrum 

    Parameters
    ----------
    x : '~astropy.Quantity' with frequency or unit
    y : '~astropy.Quantity' with power or temperature unit
        
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
        
        if x.units in [u.Hz, u.kHz, u.MHz, u.GHz]:
            self.freq = x
            if 'restfreq' in kwargs.keys():
                self.velo =  (1 - self.freq / kwargs['restfreq']) * const.c
        elif x.units in [u.m/u.s, u.km/u.s]:
            self.velo = x
            if 'restfreq' in kwargs.keys():
                self.freq = (1 - self.velo / const.c) * kwargs['restfreq']
        else:
            raise ValueError('x must be in frequency or velocity unit')

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
        spect_shifted : an instance of Spectrum with shifted x axes
        '''
        spec_shifted = self
        spec_shifted.y = self.y
        if 'freq' in self.__dict__:
            spec_shifted.freq = self.freq * vlsr / const.c + self.freq
        if 'velo' in self.__dict__:
            spec_shifted.velo = self.velo - vlsr



