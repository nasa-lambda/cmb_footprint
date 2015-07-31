# pylint: disable=E1101
# pylint: disable=C0325

'''
=============================================================
footprint.py : Survey footprint related classes and functions
=============================================================

This module provides the class which we use to generate a survey footprint
'''

from __future__ import print_function

import os
import inspect

import numpy as np
import pylab as pl
import healpy as H
import matplotlib.cm as cm

import cmb_footprint.util as util
from cmb_footprint.config_handler import ConfigHandler


class SurveyStack(object):
    '''The SurveyStack class allows us to overlay survey hitmaps upon
    an input background. The hitmaps have their transparency associated with
    the number of hits (more hits = less transparent).
    '''

    def __init__(self, background, xsize=1600, nside=None, fignum=None,
                 projection=H.mollview, coord_bg='G', coord_plot='C',
                 rot=None, partialmap=False, latra=None, lonra=None,
                 config='footprint.cfg', map_path=None,
                 download_config=False):

        self.xsize = xsize
        self.fig = pl.figure(fignum)
        self.coord_plot = coord_plot
        self.partialmap = partialmap
        self.latra = latra
        self.lonra = lonra
        self.rot = rot
        self.mapview = projection
        self.cbs = []

        if map_path is None:
            full_path = inspect.getfile(inspect.currentframe())
            abs_path = os.path.split(full_path)[0]
            map_path = os.path.join(abs_path, 'maps/')

        self.config = ConfigHandler(config, map_path, nside=nside,
                                    download_config=download_config)

#       Could also just call load_survey which will call get_background
        if isinstance(background, str):
            background, coord_bg = self.config.load_survey(background)

        if nside is None:
            nside = H.npix2nside(len(background))

        self.nside = nside

        coord = [coord_bg, coord_plot]

        cm.Greys.set_under(alpha=0.0)

        title = 'Survey Footprints'

        if self.partialmap:
            sub = (1, 1, 1)
            margins = (0.01, 0.025, 0.01, 0.03)
            H.cartview(background, title=title, xsize=self.xsize, coord=coord,
                       fig=self.fig.number, cmap=cm.Greys, norm='log',
                       notext=True, rot=rot, flip='astro', min=1.0, max=5000,
                       latra=latra, lonra=lonra, sub=sub, margins=margins)
            self.fig.delaxes(self.fig.axes[-1])
        else:
            self.mapview(background, title=title, xsize=self.xsize,
                         coord=coord, fig=self.fig.number, cmap=cm.Greys,
                         norm='log', min=1.0, max=5000, notext=True,
                         cbar=None, rot=rot, flip='astro')

        H.graticule(dpar=30.0, dmer=30.0, coord='C', verbose=False)

    def superimpose_hpxmap(self, hpx_map, label, color='red', coord_in='C'):
        '''Superimpose a Healpix map on the background map.

        Parameters
        ----------
        hpx_map : array-like
            The hpx_map to superimpose upon the background.

        label : string
            The label to put on the colorbar for this footprint.

        color : string or array-like with shape (3,)
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet.

        coord_in : 'C', 'G', or 'E'
            The coordinate system of the input healpix map.

        Notes
        -----
        The input healpix map will have zeros replaced with NaNs to make
        those points completely transparent.
        '''

        idx_nan = (hpx_map == 0)
        hpx_map /= np.max(hpx_map)
        hpx_map[idx_nan] = np.NaN

        cm1 = util.get_color_map(color)

        coord = [coord_in, self.coord_plot]

        if self.partialmap:
#           Colorbar is added to this and then deleted to make sure there is
#           room at the bottom of the map for the labels. Margins are to make
#           sure the title is not partially off the figure for a square map
            sub = (1, 1, 1)
            margins = (0.01, 0.025, 0.01, 0.03)
            map_tmp = H.cartview(hpx_map, title='', xsize=1600, coord=coord,
                                 fig=self.fig.number, cmap=cm1, notext=True,
                                 flip='astro', rot=self.rot, latra=self.latra,
                                 lonra=self.lonra, sub=sub, margins=margins,
                                 return_projected_map=True)
            idx = np.isfinite(map_tmp)
            add_cb = len(map_tmp[idx]) > 0
            self.fig.delaxes(self.fig.axes[-1])
        else:
            self.mapview(hpx_map, title='', xsize=1600, coord=coord,
                         cbar=None, fig=self.fig.number, cmap=cm1,
                         notext=True, flip='astro', rot=self.rot)
            add_cb = True

        if add_cb:
#           First add the new colorbar axis to the figure
            im0 = self.fig.axes[-1].get_images()[0]
            box = self.fig.axes[0].get_position()
            ax_color = pl.axes([len(self.cbs), box.y0-0.1, 0.05, 0.05])
            self.fig.colorbar(im0, cax=ax_color, orientation='horizontal',
                              label=label, values=[2, 3])

            self.cbs.append(ax_color)

#           Readjust the location of every colorbar
            ncb = len(self.cbs)

            left = 1.0 / (2.0*ncb) - 0.025
            for ax_tmp in self.cbs:
                ax_tmp.set_position([left, box.y0-0.1, 0.05, 0.05])
                left += 1.0 / ncb

    def superimpose_fits(self, fns, label, color='red', maptype='WCS',
                         coord_in='C'):
        '''Superimpose the footprint of a survey on the background image.
        Can be a single fits file or a list of them that will be added
        together.

        Parameters
        ----------
        fns : list of strings
            The filenames for the maps to read in

        color : string or array-like with shape (3,)
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet.

        maptype : string
            'WCS' or 'HPX' describing the type of map in the FITS file.

        coord_in : 'C', G', or 'E'
            Coordinate system of the input map
        '''

        if maptype == 'WCS':
            hpx_map = util.read_wcs_maps(fns, 256)
        elif maptype == 'HPX':
            hpx_map = util.read_hpx_maps(fns)

        self.superimpose_hpxmap(hpx_map, label, color=color, coord_in=coord_in)

    def superimpose_bound_cen(self, center, size, label,
                              color='red', coord_in='C'):
        '''Superimpose the footprint of a survey on the background image
        by giving input radec boundaries for the map. Boundaries are defined
        as the center and "radius" in ra/dec.

        Parameters
        ----------
        center : array-like with shape (2,)
            The center of the survey (degrees). ra/dec, gall/galb, etc.

        size : array-like with shape (2,)
            The length of the edge of the rectangle in degrees

        label : string
            The label to put on the colorbar for this survey

        color : string or array-like with shape (3,)
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet. Default : 'red'

        coord_in : 'C', 'G', or 'E'
            The coordinate system of the input parameters. 'C' would mean
            input values are in ra,dec. Default : 'C'
        '''

        hpx_map = util.gen_map_centersize(center, size, self.nside)

        self.superimpose_hpxmap(hpx_map, label, color=color,
                                coord_in=coord_in)

    def superimpose_bound_circ(self, center, radius, label,
                               color='red', coord_in='C'):
        '''Superimpose the footprint of a survey on the background image
        by giving an input center ra/dec and a radius of a disc.

        Parameters
        ----------
        center : array-like with shape (2,)
            The center of the survey (degrees). ra/dec, gall/galb, etc.

        radius : float
            The radius of the disc (degrees)

        label : string
            The label to put on the colorbar for this survey

        color : string or array-like with shape (3,)
            The color for the survey. Either a string recognized by
            matplotlib or a rgb triplet. Default : 'red'

        coord_in : 'C', 'E', or 'G'
            The coordinate system of the input values. Default : 'C'
        '''

        hpx_map = util.gen_map_disc(center, radius, self.nside)

        self.superimpose_hpxmap(hpx_map, label, color=color,
                                coord_in=coord_in)

    def superimpose_bound_vtx(self, vertices, label, color='red',
                              coord_in='C'):
        '''Superimpose the footprint of a survey on the background image
        by giving the ra/dec corners of the image. The enclosed survey
        footprint is generated by calling healpy.query_polygon.

        Parameters
        ----------
        vertices : array-like with shape (n,2)
            The n corners of the survey footprint in degrees.

        label : string
            The label to put on the colorbar for this survey

        color : string or array-like with shape (3,)
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet. Default : 'red'

        coord_in : 'C', 'E', or 'G'
            The coordinate system of the input vertices. Default : 'C'
        '''

        hpx_map = util.gen_map_polygon(vertices, self.nside)

        self.superimpose_hpxmap(hpx_map, label, color=color,
                                coord_in=coord_in)

    def superimpose_survey(self, survey_name, color='red',
                           label=None):
        '''Superimpose a specific survey whose Healpix footprints we have
        pregenerated and are listed in the configuration file

        Parameters
        ----------
        survey_name : string
            Name of survey. Valid values are section names in the
            configuration file.

        color : string or array-like with shape (3,)
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet.

        label : string
            The label for the survey. If none, survey_name is used as
            the label.
        '''

        hpx_map, coord = self.config.load_survey(survey_name)

        if label is None:
            label = survey_name

        self.superimpose_hpxmap(hpx_map, label, color=color,
                                coord_in=coord)

    def superimpose_rect_outline(self, lonra, latra, color='red',
                                 label=None):

        linelon = np.linspace(lonra[0], lonra[1], num=1000)
        linelat = latra[0]*np.ones_like(linelon)
        H.projplot(linelon, linelat, '.', lonlat=True, markersize=1,
                   color=color)

        linelon = np.linspace(lonra[0], lonra[1], num=1000)
        linelat = latra[1]*np.ones_like(linelon)
        H.projplot(linelon, linelat, '.', lonlat=True, markersize=1,
                   color=color)

        linelat = np.linspace(latra[0], latra[1], num=1000)
        linelon = lonra[0]*np.ones_like(linelat)
        H.projplot(linelon, linelat, '.', lonlat=True, markersize=1,
                   color=color)

        linelat = np.linspace(latra[0], latra[1], num=1000)
        linelon = lonra[1]*np.ones_like(linelat)
        H.projplot(linelon, linelat, '.', lonlat=True, markersize=1,
                   color=color)

        add_cb = False
        if add_cb:
#           First add the new colorbar axis to the figure
            im0 = self.fig.axes[-1].get_images()[0]
            box = self.fig.axes[0].get_position()
            ax_color = pl.axes([len(self.cbs), box.y0-0.1, 0.05, 0.05])
            self.fig.colorbar(im0, cax=ax_color, orientation='horizontal',
                              label=label, values=[2, 3])

            self.cbs.append(ax_color)

#           Readjust the location of every colorbar
            ncb = len(self.cbs)

            left = 1.0 / (2.0*ncb) - 0.025
            for ax_tmp in self.cbs:
                ax_tmp.set_position([left, box.y0-0.1, 0.05, 0.05])
                left += 1.0 / ncb
