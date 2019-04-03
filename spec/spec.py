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
        
        if 'restfreq' in kwargs.keys():

