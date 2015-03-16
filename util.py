# pylint: disable=E1101
# pylint: disable=C0325

'''This module provides utility functions needed for the footprint code.'''

import numpy as np
import healpy as H
import astropy.wcs as wcs
from scipy.sparse import coo_matrix


def wcs_to_healpix(hdulist, nside):
    '''Converts data in an opened FITS file from a generic WCS to Healpix'''

    wcs1 = wcs.WCS(hdulist[0].header)

    data = hdulist[0].data

    print("x/ysize = ", data.shape[0], data.shape[1])

    pixcrd = [(y, x) for x in xrange(data.shape[0]) for y in
              xrange(data.shape[1])]

#   pix2world to get ra/dec values of each pixel
    print("wcs_pix2world")
    world = wcs1.wcs_pix2world(pixcrd, 1)

    print("np.array")
    pixcrd = np.array(pixcrd)

#   needed since in some ACT data the length of the map is >360 degrees
    idx = np.isfinite(world[:, 0])

    pixcrd = pixcrd[idx, :]
    world = world[idx, :]

#   Use Healpy to do ang2pix for Healpix pixel numbers
    theta = np.pi/2.0 - np.radians(world[:, 1])
    phi = np.radians(world[:, 0])
    print("ang2pix")
    pixnum = H.ang2pix(nside, theta, phi)

    print("hpx_data")
    npix = H.nside2npix(nside)
    row = pixnum
    col = np.zeros_like(row)

#   coo_matrix will sum entries that have multiple elements
    hpx_data = coo_matrix((data[(pixcrd[:, 1], pixcrd[:, 0])], (row, col)),
                          shape=np.array([npix, 1])).A.transpose()
    hpx_data.shape = (npix)

    return hpx_data
