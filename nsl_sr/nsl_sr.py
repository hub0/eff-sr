#!python

def calculate_nsl_sr(obsdir, obstime, pa, spec, tele_loc, nsl_pat, 
                     sky_model, sky_model_spec):
    '''
    calculate near side lobe stray radiation spectrum

    Parameters
    ----------
    obsdir : '~astropy.coordinates.SkyCoord'
        sky direction of the beam
    obstime : '~astropy.time'
        observation time
    pa : float
        position angle of the beam in degree
    spec : '~numpy.ndarray'
        frequency scheme on which the sr spectrum is to be built
    tele_loc : '~astropy.coordinates.EarthLocation'
        location of the telescope
    nsl_pat : '~numpy.ndarray'
        beam pattern of near side lobe with columns
        az az zen zen gain
    sky_model : '~numpy.ndarray'
        healpix array of sky model in ring order
    sky_model_spec : '~numpy.ndarray'
        frequency scheme of sky model

    Returns
    -------
    nsl_sr : '~numpy.ndarray'
        near side lobe stray radiation spectrum with the same frequency
        scheme of sky model 
    '''

    import numpy as np
    from scipy.integrate import quad
    from astropy.time import Time
    import astropy.coordinates as coord
    from astropy.coordinates.representation import UnitSphericalRepresentation as USphR
    from astropy.coordinates.representation import CartesianRepresentation as CartR
    from astropy.coordinates.matrix_utilities import rotation_matrix as rotm
    import astropy.units as u
    import healpy as hp

    nsl_az = np.deg2rad(nsl_pat.T[:2])
    nsl_zen = np.deg2rad(nsl_pat.T[2:4])
    nsl_gain = np.array(nsl_pat.T[5])
    nsl_solan = np.zeros(len(nsl_gain))
    # calculate pattern sector solid angle
    for i in np.arange(len(nsl_gain)):
        az0 = nsl_az[0, i]
        az1 = nsl_az[1, i]
        zen0 = nsl_zen[0, i]
        zen1 = nsl_zen[1, i]
        nsl_solan[i] = quad(lambda zen: np.sin(zen), zen0, zen1)[0] * 
                       (az1 - az0)
    
    # position of the center of sectors
    nsl_az = nsl_az.mean(0)
    nsl_zen = nsl_zen.mean(0)

    # from beam spherical to beam cartesian
    nsl_sph = USphR(nsl_az * u.deg, 90 * u.deg - nsl_zen * u.deg)
    nsl_cart = nsl_sph.to_cartesian()
    nsl_cart = nsl_cart.transform(rotm(pa, axis='z')) # restore the pa
    
    # from beam cartesian to galactic spherical
    obs_dec = obsdir.icrs.dec
    obs_ra = obsdir.icrs.ra
    nsl_cart = nsl_cart.transform(rotm(90 * u.deg - obs_dec, axis='y'))
    nsl_cart = nsl_cart.transform(rotm(180 * u.deg - obs_ra, axis='z'))
    nsl_sph = coord.cartesian_to_spherical(nsl_cart.x, sky_cart.y, sky_cart.z)
    nsl_coord = coord.SkyCoord(nsl_sph[2], nsl_sph[1], frame='icrs',
                               obstime=obstime, location=tele_loc)
    gal = beam_coord.galactic
    l = gal.l.rad
    b = np.pi/2 - gal.b.rad

    # get sky model value
    model_npix = len(sky_model)
    model_nside = hp.npix2nside(model_npix)
    idx = hp.ang2pix(model_nside, b, l)
    nsl_model = sky_model(idx) # 2-d array [position, frequency]
    
    # extinction
    nsl_ext = cal_extinction(nsl_zen)

    # weight
    nsl_wei = cal_weight(nsl_az, nsl_zen)

    # the vlsr-uncorrected SR spectrum
    nsl_sr_uncorr = nsl_model * nsl_gain * nsl_wei * nsl_ext * nsl_solan

    # vlsr correction
    obsvlsr = calculate_lsr_correction(obsdir, tele_loc, obstime) # in km/s
    nsl_sr = spec_shift(nsl_sr_uncorr, obsvlsr)

    # regrid nsl_sr to wanted frequency scheme
    nsl_sr = freq_regrid(nsl_sr, spec)

    return nsl_sr


def cal_extinction(zen):
    '''
    calculate the atmophere extinction for effelsberg

    Parameters
    ----------
    zen: float
        zenith angle

    Returns
    -------
    extinc: float
        the extinction coefficient
    '''
    
