import os

import numpy as np
import pylab as pl
import healpy as H
from astropy.io import fits

def bin_catalog(ra_rad, dec_rad, redshift, nside, z_left, z_right):
    npix = H.nside2npix(nside)
    nnu = len(z_left)
    # from Ra/dec to galactic
    rotate = H.rotator.Rotator(coord=['C','C'])
    theta_gal, phi_gal = rotate(dec_rad, ra_rad)

    gal_ind = H.pixelfunc.ang2pix(nside, theta_gal, phi_gal, 
                                       nest=False)

    # spatial density
    gal_spatial = np.bincount(gal_ind, minlength=npix)

    # z binning
    gal_ring = np.zeros(shape=(npix, nnu))
    gal_counts = np.zeros_like(z_left)
    for ind in range(nnu):
        in_bin = np.logical_and(redshift > z_left[ind], redshift < z_right[ind])
        gal_bin = gal_ind[in_bin]
        gal_ring[:,ind] = np.bincount(gal_bin, minlength=npix)
        gal_counts[ind] = len(gal_bin)

    # make a separable selection function
    nbar = gal_spatial[:, None] * gal_counts[None, :]
    nbar *= np.sum(gal_counts) / np.sum(nbar)

    overdensity = (gal_ring - nbar) / nbar

    return overdensity, nbar, gal_counts, gal_spatial

def load_qso(root_data):
    #qso_sample = root_data + 'BOSS/dr10/boss/qso/DR10Q/DR10Q_v2.fits'
    # http://data.sdss3.org/sas/dr12/boss/qso/DR12Q/
    qso_sample = os.path.join(root_data,'BOSS/dr12/boss/qso/DR12Q/DR12Q.fits')

    gal_data = fits.open(qso_sample)
    gal_dec = -gal_data[1].data['DEC'][:]* np.pi / 180. + np.pi / 2.
    gal_ra = gal_data[1].data['RA'][:] * np.pi / 180.
    gal_z = gal_data[1].data['Z_PCA'][:]

    return gal_ra, gal_dec, gal_z

def load_twompz(root_data):
    # SELECT * FROM TWOMPZ..twompzPhotoz
    # http://surveys.roe.ac.uk/ssa/sql.html
    twompz_sample = os.path.join(root_data,"2MASS_combined_photoz/results15_14_50_41_9.fits.gz")

    gal_data = fits.open(twompz_sample)
    gal_dec = -gal_data[1].data['DEC'][:] + np.pi / 2.
    gal_ra = gal_data[1].data['RA'][:]
    gal_z = gal_data[1].data['ZPHOTO'][:]

    return gal_ra, gal_dec, gal_z

def load_boss(root_data):
    cmassn_sample = os.path.join(root_data,'BOSS/dr10/boss/lss/galaxy_DR10v8_CMASS_North-specObj.fits.gz')
    rcmassn_sample = os.path.join(root_data,'BOSS/dr10/boss/lss/random2_DR10v8_CMASS_North.fits.gz')
    lowzn_sample = os.path.join(root_data,'BOSS/dr10/boss/lss/galaxy_DR10v8_LOWZ_North-specObj.fits.gz')
    rlowzn_sample = os.path.join(root_data,'BOSS/dr10/boss/lss/random2_DR10v8_LOWZ_North.fits.gz')

    cmasss_sample = os.path.join(root_data,'BOSS/dr10/boss/lss/galaxy_DR10v8_CMASS_South-specObj.fits.gz')
    rcmasss_sample = os.path.join(root_data,'BOSS/dr10/boss/lss/random2_DR10v8_CMASS_South.fits.gz')
    lowzs_sample = os.path.join(root_data,'BOSS/dr10/boss/lss/galaxy_DR10v8_LOWZ_South-specObj.fits.gz')
    rlowzs_sample = os.path.join(root_data,'BOSS/dr10/boss/lss/random2_DR10v8_LOWZ_South.fits.gz')

    gal_dec = np.array([])
    gal_ra = np.array([])
    gal_z = np.array([])
    
    #ra_key = 'RA'
    #dec_key = 'DEC'
    #z_key = 'Z'
    
    ra_key = 'PLUG_RA'
    dec_key = 'PLUG_DEC'
    z_key = 'Z'

    gal_data = fits.open(lowzn_sample)
    gal_dec = np.append(gal_dec, -gal_data[1].data[dec_key][:] * np.pi / 180. + np.pi / 2.)
    gal_ra = np.append(gal_ra, gal_data[1].data[ra_key][:] * np.pi / 180.)
    gal_z = np.append(gal_z, gal_data[1].data[z_key][:])
    gal_data.close()

    gal_data = fits.open(lowzs_sample)
    gal_dec = np.append(gal_dec, -gal_data[1].data[dec_key][:] * np.pi / 180. + np.pi / 2.)
    gal_ra = np.append(gal_ra, gal_data[1].data[ra_key][:] * np.pi / 180.)
    gal_z = np.append(gal_z, gal_data[1].data[z_key][:])
    gal_data.close()

    gal_data = fits.open(cmassn_sample)
    gal_dec = np.append(gal_dec, -gal_data[1].data[dec_key][:] * np.pi / 180. + np.pi / 2.)
    gal_ra = np.append(gal_ra, gal_data[1].data[ra_key][:] * np.pi / 180.)
    gal_z = np.append(gal_z, gal_data[1].data[z_key][:])
    gal_data.close()

    gal_data = fits.open(cmasss_sample)
    gal_dec = np.append(gal_dec, -gal_data[1].data[dec_key][:] * np.pi / 180. + np.pi / 2.)
    gal_ra = np.append(gal_ra, gal_data[1].data[ra_key][:] * np.pi / 180.)
    gal_z = np.append(gal_z, gal_data[1].data[z_key][:])
    gal_data.close()

    return gal_ra, gal_dec, gal_z

if __name__ == '__main__':

    root_data = './maps/'

    ra, dec, z = load_boss(root_data)

    data = bin_catalog(ra,dec,z,128,[0],[100])

    print len(ra), np.sum(data[3])

    H.write_map('maps/BOSS_dr10_lss.fits', data[3])
   
    H.mollview(data[3])
    pl.show()

