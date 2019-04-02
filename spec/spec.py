import numpy as np
import astropy.units as u

__all__ = ['Spectrum']

class Spectrum():
    '''
    Object for implement radio spectrum 

    Parameters
    ----------
    freq : '~astorpy.Quantity' with frequency unit
    value : '~astropy.Quantity' with power or temperature unit
        
    '''

    def __init__(self, freq, value):
        if not isinstance(freq, u.Quantity):
            raise TypeError('freq must be an astropy.units.Quantity instance')
        if not isinstance(value, u.Quantity):
            raise TypeError('value must be an astropy.units.Quantity instance')
        if freq.shape == value.shape:
            raise ValueError('shape of freq and value are not identical {0} {1}'.format(freq.shape, value.shape))

