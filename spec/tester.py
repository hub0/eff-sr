import numpy as np
import astropy.units as u
from spec import *

a = np.array([1,2,3,4,5])
b = a + np.array([ -0.1, 0.2, -0.1, -0.05, 0.2])
x0 = a * u.GHz
x1 = b * u.GHz
y0 = a * 2
y1 = np.ones(5)

spec0 = Spectrum(x0, y0, restfreq=2*u.GHz)
spec1 = Spectrum(x1, y1)


