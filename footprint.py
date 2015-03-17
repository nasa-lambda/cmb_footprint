# pylint: disable=E1101
# pylint: disable=C0325

'''The purpose of this module is to allow us to plot footprints of different
experiments upon an input background map.
'''

import numpy as np
import pylab as pl
import healpy as H
import matplotlib.cm as cm
from astropy.io import fits

import util


class SurveyStack(object):
    '''The SurveyStack class allows us to overlay experiment hitmaps upon
    an input background. The hitmaps have their transparency associated with
    the number of hits (more hits = less transparent).
    '''

    def __init__(self, background, xsize=800, nside=None, fignum=None):
        self.xsize = xsize
        self.fn_background = background
        self.fig = pl.figure(fignum)

#       Q,U should be fields 1 and 2 for a normal Healpix map
        maps = H.read_map(background, field=(1, 2))

        if nside is None:
            nside = H.npix2nside(len(maps[0]))

        self.nside = nside

        polamp = np.sqrt(maps[0]**2 + maps[1]**2)

        cm.Greys.set_under(alpha=0.0)
        H.mollview(polamp, title='Experiment Footprints', xsize=1600,
                   coord=['G', 'C'], fig=self.fig.number, cmap=cm.Greys,
                   max=0.001, min=0.0, notext=True, cbar=None, flip='astro')
        H.graticule(dpar=15.0, dmer=15.0, coord='C')

    def read_hpx_maps(self, fns):
        '''Read in one or more healpix maps and add them together. Must input
        an array of strings even if only inputting a single map.

        Parameters
        ----------
        fns: list of strings
            The filenames for the healpix maps to read in.

        Returns
        -------
        hpx_map: array-like
            A healpix map'''

        hpx_map = np.zeros(H.nside2npix(self.nside))
        for fn_tmp in fns:
            tmp_map = H.read_map(fn_tmp)
            hpx_map += H.ud_grade(tmp_map, self.nside)

        return hpx_map

    def read_wcs_maps(self, fns):
        '''Read in WCS FITS files and convert them to a Healpix map.

        Parameters
        ----------
        fns: list of strings
            The filenames for the WCS maps to read in.

        Returns
        -------
        hpx_map: array-like
            The WCS maps read in and converted to healpix format
        '''

        hpx_map = np.zeros(H.nside2npix(self.nside))
        for fn_tmp in fns:
            hdulist = fits.open(fn_tmp)
            hpx_map += util.wcs_to_healpix(hdulist, self.nside)
            hdulist.close()

        return hpx_map

    def superimpose_hpxmap(self, hpx_map, color='red'):
        '''Superimpose a Healpix map on the background map.

        Parameters
        ----------
        hpx_map: array-like
            The hpx_map to superimpose upon the background.

        color: string or array-like with shape (3,)
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet.

        Notes
        -----
        The input healpix map will replace zeros in the map with NaNs to make
        those points completely transparent.
        '''

        idx_nan = (hpx_map == 0)
        hpx_map[idx_nan] = np.NaN

        cm1 = get_color_map(color)

        H.mollview(hpx_map, title='', xsize=1600, coord='C',
                   fig=self.fig.number, cmap=cm1, notext=True, cbar=None,
                   flip='astro')

#       This code is what is needed if you can't call mollview()
#       extent = (0.02, 0.05, 0.96, 0.9)

#       ax = H.projaxes.HpxMollweideAxes(self.fig, extent, coord='C',
#                                        rot=None, format='%g',
#                                        flipconv='astro')
#       ax.projmap(hpx_map, xsize=1600, coord='C', cmap=cm1)
#       self.fig.add_axes(ax)

    def superimpose_fits(self, fns, color='red', maptype='WCS'):
        '''Superimpose the footprint of an experiment on the background image.
        Can be a single fits file or a list of them that will be added
        together.

        Parameters
        ----------
        fns: list of strings
            The filenames for the maps to read in

        color: string or array-like with shape (3,)
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet.

        maptype: string
            'WCS' or 'HPX' describing the type of map in the FITS file.
        '''

        if maptype == 'WCS':
            hpx_map = self.read_wcs_maps(fns)
        elif maptype == 'HPX':
            hpx_map = self.read_hpx_maps(fns)

        self.superimpose_hpxmap(hpx_map, color=color)

    def superimpose_boundary_cen(self, radec_cen, radec_size, color='red'):
        '''Superimpose the footprint of an experiment on the background image
        by giving input radec boundaries for the map. Boundaries are defined
        as the center and "radius" in ra/dec.

        Parameters
        ----------
        radec_cen: array-like with shape (2,)
            The center ra/dec of the survey (degrees)

        radec_size: array-like with shape (2,)
            The "radius" of the survey in ra/dec (degrees)

        color: string or array-like with shape (3,)
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet.
        '''

        corner1 = (radec_cen[0]+radec_size[0], radec_cen[1]+radec_size[1])
        corner2 = (radec_cen[0]+radec_size[0], radec_cen[1]-radec_size[1])
        corner3 = (radec_cen[0]-radec_size[0], radec_cen[1]-radec_size[1])
        corner4 = (radec_cen[0]-radec_size[0], radec_cen[1]+radec_size[1])

        corners = (corner1, corner2, corner3, corner4)

        self.superimpose_boundary_corners(corners, color=color)

    def superimpose_boundary_corners(self, radec_corners, color='red'):
        '''Superimpose the footprint of an experiment on the background image
        by giving the ra/dec corners of the image. The enclosed survey
        footprint is generated by calling healpy.query_polygon.

        Parameters
        ----------
        radec_corners: array-like with shape (n,2)
            The n corners of the survey footprint in degrees.

        color: string or array-like with shape (3,)
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet.
        '''

        radec_corners = np.array(radec_corners)

        thetas = np.pi/2 - np.radians(radec_corners[:, 1])
        phis = np.radians(radec_corners[:, 0])

        vecs = H.ang2vec(thetas, phis)

        ipix = H.query_polygon(self.nside, vecs)

        hpx_map = np.zeros(H.nside2npix(self.nside))
        hpx_map[ipix] = 1.0

        self.superimpose_hpxmap(hpx_map, color=color)

    def superimpose_experiment(self, experiment_name, color='red'):
        '''Superimpose a specific experiment whose Healpix footprints we have
        pregenerated.

        Parameters
        ----------
        experiment_name: string
            Currently only 'ACT' is valid

        color: string or array-like with shape (3,)
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet.
        '''

        if experiment_name == 'ACT':
            fns = ['maps/ACT_148_equ_hits_hpx.fits',
                   'maps/ACT_148_south_hits_hpx.fits']
        else:
            print('We do not have Healpix maps for this experiment')

        self.superimpose_fits(fns, color=color, maptype='HPX')


def get_color_map(color):
    '''Generate a LinearSegmentedColormap with a single color and varying
    transparency. Bad values and values below the lower limit are set to be
    completely transparent.

    Parameters
    ----------
    color: string or array-like with shape (3,)
        The color to use when overlaying the survey footprint. Either a
        string that can be input to colorConverter.to_rgb() or rgb triplet.

    Returns
    -------
    colormap1: LinearSegmentedColormap
        A color map that is a single color but varies its transparency
        from 0.5 to 1. Bad values and values below the minimum are completely
        transparent.
    '''

    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.colors import colorConverter

#   if the type is a string it is assumed to be an input to allow us
#   to get an rgb value. If it is not a string and length is 3, it is
#   assumed to be an actual rgb value
    if isinstance(color, str):
        rgb = colorConverter.to_rgb(color)
    elif len(color) == 3:
        rgb = color

    cdict = {'red':   [(0, rgb[0], rgb[0]),
                       (1, rgb[0], rgb[0])],
             'green': [(0, rgb[1], rgb[1]),
                       (1, rgb[1], rgb[1])],
             'blue':  [(0, rgb[2], rgb[2]),
                       (1, rgb[2], rgb[2])],
             'alpha': [(0, 0.5, 0.5),
                       (1, 1, 1)]}

    colormap1 = LinearSegmentedColormap('FootprintCM', cdict)
    colormap1.set_bad(alpha=0.0)
    colormap1.set_under(alpha=0.0)

    return colormap1
