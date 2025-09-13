#keeping post simulation calculations in a seperate file

import numpy as np

#parameters: two arrays,
# first list of masses in core
#second list of masses in scavenge
def captureEfficiency(core,scavenge):
    """
    Calculate capture efficiency

    Parameters
    ----------
    core : array of integers(masses)
    scavenge : array of integers(masses)

    Returns
    ----------
    efficiency : float
        capture efficiency
    """
    #Calulate mass in each outlet--------------
    #scavenge
    scavenge_mass = np.sum(scavenge)
    #core
    core_mass = np.sum(core)

    #calculate capture efficiency---------------
    efficiency = (scavenge_mass / (scavenge_mass + core_mass))
    return efficiency