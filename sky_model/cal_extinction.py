import numpy as np

def cal_extinction(zen):
    '''
    calculate atmosphere extinction for given zenith angle
    
    The Effelsberg atmospherical model is provided by Peter Kalberla
    ext = 1.012745 ** (-X)  with
    X=(1.1033422 + 0.0828460*sec + 0.0759914*sec*sec + 0.00000391*sec*sec*sec) 
    sec = secans (zenith distance)

    Parameters
    ----------
    zen : float
        zenith angle in degree
    
    Returns
    -------
    ext : float
        extinction coefficiency
    '''

    sec = np.cos(np.deg2rad(zen))
    index = 1.1033422 
            + 0.0828460*sec 
            + 0.0759914*sec*sec 
            + 0.00000391*sec*sec*sec
    ext = 1.012745 ** (-1 * index)

    return ext

