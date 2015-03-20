#!/usr/bin/env python

import argparse

from astropy.io import fits
import healpy as H
import numpy as np

import util

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Convert hitmaps to Healpix format to be used in footprint code')
    parser.add_argument('-i', dest='fn_in', action='store', nargs='+', help='Input maps')
    parser.add_argument('-o', dest='fn_out', action='store', help='Output map', default='map_out.fits')
    parser.add_argument('-nside', dest='nside', action='store', type=int, help='Nside for output map', default=512)

    args = parser.parse_args()

    print "Converting", len(args.fn_in), "maps to a single Healpix map"

    hpx_map = np.zeros(H.nside2npix(args.nside))
    for fn in args.fn_in:
        hdulist = fits.open(fn)
        print hdulist
        hpx_map += util.wcs_to_healpix(hdulist,args.nside)

    H.write_map(args.fn_out,hpx_map)

