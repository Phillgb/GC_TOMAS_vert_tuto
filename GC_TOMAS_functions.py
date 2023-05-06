#!/usr/bin/python3
# Author : Phillipe Gauvin-Bourdon
# Created : April 03, 2023
# ------------------------------------------------------------------------------
# Import modules
# ------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------------------------
# Function definition
# ------------------------------------------------------------------------------
def GC_sig2press(levels=None, surf=None, top=None):
    '''Function converting GEOS-Chem sigma pressure levels to air pressure
    
    PARAMETERS:
    - levels = Array containning the sigma pressure levels to be converted
    - psurf = Surface pressure value in mb
    - ptop = Top of the atmosphere pressure value in mb
    
    RETURN:
    - press = Numpy array containning the pressure levels in mb'''

    press = []
    for l in levels:
        press.append(top + l * (surf - top))
    press = np.array(press)

    return press

