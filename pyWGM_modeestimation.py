import numpy as np
import os
import csv
from scipy.constants import speed_of_light,pi

def estnPhi(pillarRadius, wavelength, permCavity):
    # Values were determined via a regression of previous sim values.
    # Formula: nPhi = (n / lambda) * 2 * (a + r * b + r^2 * c) * pi
    # The modes to not propagate along the pillar edges, but at a certain distance,
    # that also changes over radius. a = constant offset, b = linear offset change over r, c = square offset change over r.
    # n = (real) refractive index of the cavity material, lambda = vacuum wavelength 
    a = -1.194e+02
    b = 8.727e-01
    c = 7.137e-06
    d = -2.950e-10

    nPhi = (np.sqrt(np.real(complex(permCavity.replace('i', 'j')))) / wavelength) * 2 * (a + pillarRadius * b + pillarRadius**2 * c + pillarRadius**3 * d) * pi
    # Return floor instead of rounded value, since the sim algo will propagate to higher nPhi
    return int(round(nPhi))

def estGuessY(pillarRadius, nPhi, permCavity):
    a = -1.194e+02
    b = 8.727e-01
    c = 7.137e-06
    d = -2.950e-10

    guess = (nPhi * speed_of_light)/( (np.sqrt(np.real(complex(permCavity.replace('i', 'j'))))) * (a + pillarRadius * b + pillarRadius**2 * c + pillarRadius**3 * d) * 1e-9 )
    return guess