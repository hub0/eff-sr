#!python
import numpy as np
from astropy.time import Time
from astropy.time import TimeDelta
from astropy.io import fits
import astropy.coordinates as coord
from astropy.coordinates.representation import UnitSphericalRepresentation as USphR
from astropy.coordinates.representation import CartesianRepresentation as CartR
from astropy.coordinates.matrix_utilities import rotation_matrix as rotm
import astropy.units as u
import healpy as hp
import matplotlib.pyplot as plt

SRMAP_NSIDE = 64 
MODEL_NSIDE = 128

########################
# set observation time #
########################
obstime = Time('2018-12-17 13:30:26.399', format='iso', scale='utc')
# time difference between local time and UTC
dt = TimeDelta(8 * 3600, format='sec', scale='tai')

# convert to UTC
obstime = obstime - dt

#############################################
# set location & position angle of the beam #
#############################################

beam_pa = 17 * u.deg # beam position angle in degree (east of north)
fast_loc = coord.EarthLocation(106 * u.deg + 51 * u.arcmin + 24.0 * u.arcsec,
        25 * u.deg + 39 * u.arcmin + 10.6 * u.arcsec,
        1110.0288 * u.m)

####################################
# read in beam pattern & sky model #
####################################
beam_pat = np.loadtxt('effberg_near_sl.dat', comments='#')
beam_pat = beam_pat.reshape(len(beam_pat), 3, 2).mean(-1).T
#beam_sph = np.deg2rad(beam_pat[0:2])
beam_usph = USphR(beam_pat[0] * u.deg, 90 * u.deg - beam_pat[1] * u.deg)
beam_gain = beam_pat[2]
# rows: position angle (phi), radius (theta), coe 

# read in sky model
nhi_hdu = fits.open('../NHI_HPX.fits')
nhi_data = nhi_hdu[1].data['NHI']

# down grade the nhi model
nhi_data = hp.ud_grade(nhi_data, MODEL_NSIDE)

################################################
# convert beam pattern to equatorial cartesian #
################################################

# from beam spherical to beam cartesian
beam_cart = beam_usph.to_cartesian()

# from beam cartesian to equatorial cartesian (rotational transform)
# rotate p with z-axis
beam_cart = beam_cart.transform(rotm(beam_pa, axis='z'))

###############################################
# build the frame of near sl SR map (healpix) #
###############################################
theta, phi= hp.pix2ang(SRMAP_NSIDE, np.arange(hp.nside2npix(SRMAP_NSIDE)))
glat = np.pi/2 - theta 
glon = phi 

obsdir = coord.SkyCoord(glon * u.rad, glat * u.rad, frame='galactic', 
        obstime=obstime, location=fast_loc)

ra = obsdir.icrs.ra
dec = obsdir.icrs.dec
nhi_sum = np.array([])

#####################################################
# calculate nhi_sum for each obervational direction #
#####################################################
#for i in np.arange(hp.nside2npix(SRMAP_NSIDE)):
for obs_ra, obs_dec in zip(ra, dec):
    sky_cart = beam_cart.transform(rotm(90 * u.deg - obs_dec, axis='y'))
    sky_cart = sky_cart.transform(rotm(180 * u.deg - obs_ra, axis='z'))
    sky_sph = coord.cartesian_to_spherical(sky_cart.x, sky_cart.y, sky_cart.z)
    beam_coord = coord.SkyCoord(sky_sph[2], sky_sph[1], frame='icrs', 
            obstime=obstime, location=fast_loc)
    gal = beam_coord.galactic
    l = gal.l.rad
    b = np.pi/2 - gal.b.rad 
    idx = hp.ang2pix(MODEL_NSIDE, b, l) 

    nhi_sum = np.append(nhi_sum, sum(nhi_data[idx] * beam_gain))
    #nhi_sum = np.append(nhi_sum, idx[0])

############################
# save and plot the result #
############################
hp.write_map('sr.fits', nhi_sum, overwrite=True)

fig = plt.figure()
ax = fig.add_subplot(111)
glon_deg = np.rad2deg(glon)
glon_deg = np.where(glon_deg < 180, glon_deg, glon_deg - 360)
glat_deg = np.rad2deg(glat)
ax.scatter(glon_deg, glat_deg, c=nhi_sum)
ax.set_ylim(-40, 40)
ax.set_xlim(-180, 180)
#ax.set_aspect('equal')
fig.show()
fig.savefig('near_sl_sr_map{0}.pdf'.format(MODEL_NSIDE))

