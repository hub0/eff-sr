import numpy as np
import astropy.units as u
from spec import *

a = np.array([1,2,3,4,5])
x = a * u.GHz
y = a * 2

spec = Spectrum(x, y, restfreq=2*u.GHz)


