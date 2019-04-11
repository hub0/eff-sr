def calculate_lsr_correction(obsdir, obsloc, obstime):
    '''
    Calculate Vlsr velocity towards certain direction at certain time

    Parameters
    ----------
    obsdir : '~astropy.coordinates.SkyCoord'
        Direction of observation
    obsloc : '~astropy.coordinates.EarthLocation'
        Location of the observatory
    obstime : '~astropy.time.Time'
        Observation time

    Returns
    -------
    vlsr_corr : float
       the velocity of LSR in respect of observatory??? in direction of
       the obsdir 
    '''

    #from astropy.coordinates import SkyCoord, EarthLocation
    #from astropy import units as u

    #effberg = EarthLocation.from_geodetic(
    #    lat=50.52483 * u.deg, lon=6.88361 * u.deg, height=416.7 * u.m
    #                                  )
    #sc = SkyCoord(ra=ra2000 * u.deg, dec=dec2000 * u.deg)

    #obstime = Time(mjd, format='mjd')
    barycorr = obsdir.radial_velocity_correction(
        'barycentric', obstime=obstime, location=obsloc
                                               )
    #barycorr.to(u.km / u.s)

    cirs = obsdir.cirs

    # this value is used for Effberg data
    VLSR_cirs = [0.29, -17.31726, 10.00141] * u.km / u.s

    v_correct = (
        VLSR_cirs[0] * np.cos(cirs.dec.rad) * np.cos(cirs.ra.rad) +
        VLSR_cirs[1] * np.cos(cirs.dec.rad) * np.sin(cirs.ra.rad) +
        VLSR_cirs[2] * np.sin(cirs.dec.rad)
        )
    vlsr_corr = -barycorr.to(u.km / u.s) - v_correct

    return vlsr_corr.to(u.km / u.s).value
