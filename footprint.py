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

import cmb_footprint.visufunc_ext as vf
import cmb_footprint.util as util
from cmb_footprint.config_handler import ConfigHandler


class SurveyStack(object):
    '''The SurveyStack class allows us to overlay survey hitmaps upon
    an input background. The hitmaps have their transparency associated with
    the number of hits (more hits = less transparent).
    '''

    def __init__(self, background, nside=None, fignum=None,
                 projection='mollweide', coord_bg='G', coord_plot='C',
                 partialmap=False, config='footprint.cfg',
                 map_path=None, download_config=False,
                 title='Survey Footprints', cbar=None, min=1.0,
                 max=5000.0, log=True, unit='', **kwds):

        self.fig = pl.figure(fignum)
        self.coord_plot = coord_plot
        self.partialmap = partialmap

        self.kwds = kwds
        self.cbs = []

        if projection == 'mollweide':
            self.mapview = H.mollview
            self.mapcontour = vf.mollcontour
        elif projection == 'cartesian':
            self.mapview = H.cartview
            self.mapcontour = vf.cartcontour
        elif projection == 'orthographic':
            self.mapview = H.orthview
            self.mapcontour = vf.orthcontour
        elif projection == 'gnomonic':
            self.mapview = H.gnomview
            self.mapcontour = vf.gnomcontour

        if map_path is None:
            full_path = inspect.getfile(inspect.currentframe())
            abs_path = os.path.split(full_path)[0]
            map_path = os.path.join(abs_path, 'maps/')

        self.config = ConfigHandler(config, map_path, nside=nside,
                                    download_config=download_config)

#       Could also just call load_survey which will call get_background
        if isinstance(background, str):
            bgmap, coord_bg, unit2 = self.config.load_survey(background,
                                                             get_unit=True)
            background = bgmap[0]
            if unit2 is not None:
                unit = unit2

        if nside is None:
            nside = H.npix2nside(len(background))

        self.nside = nside

        coord = [coord_bg, coord_plot]

        cm.Greys.set_under(alpha=0.0)

        if log:
            min = np.log(min)
            max = np.log(max)
            unit = r'$\log($' + unit + r'$)$'
            background = np.log(background)

        if self.partialmap:
            sub = (1, 1, 1)
            margins = (0.01, 0.025, 0.01, 0.03)
            H.cartview(background, title=title, coord=coord,
                       fig=self.fig.number, cmap=cm.Greys,
                       notext=True, flip='astro', min=min, max=max,
                       sub=sub, margins=margins, **kwds)
            self.fig.delaxes(self.fig.axes[-1])
        else:
            self.mapview(background, title=title,
                         coord=coord, fig=self.fig.number, cmap=cm.Greys,
                         min=min, max=max, notext=True,
                         cbar=True, flip='astro', unit=unit, **kwds)
            if not cbar:
                self.fig.delaxes(self.fig.axes[-1])

        H.graticule(dpar=30.0, dmer=30.0, coord='C', verbose=False)

    def superimpose_hpxmap(self, hpx_map, label, color='red', coord_in='C',
                           cbar=True):
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
            map_tmp = H.cartview(hpx_map, title='',
                                 coord=coord, fig=self.fig.number, cmap=cm1,
                                 notext=True, flip='astro', sub=sub,
                                 margins=margins, return_projected_map=True,
                                 **self.kwds)
            idx = np.isfinite(map_tmp)
            cbar = len(map_tmp[idx]) > 0
        else:
            self.mapview(hpx_map, title='', coord=coord,
                         cbar=True, fig=self.fig.number, cmap=cm1,
                         notext=True, flip='astro', **self.kwds)
        self.fig.delaxes(self.fig.axes[-1])

        if cbar:
#           First add the new colorbar axis to the figure
            im0 = self.fig.axes[-1].get_images()[0]
            box = self.fig.axes[0].get_position()
            ax_color = pl.axes([len(self.cbs), box.y0-0.1, 0.05, 0.05])
            self.fig.colorbar(im0, cax=ax_color, orientation='horizontal',
                              label=label, values=[1, 1])

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
                           label=None, cbar=True):
        '''Superimpose a specific survey whose Healpix footprints we have
        pregenerated and are listed in the configuration file

        Parameters
        ----------
        survey_name : string
            Name of survey. Valid values are section names in the
            configuration file.

        color : string or array-like with shape (3,), optional
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet.

        label : string
            The label for the survey. If none, survey_name is used as
            the label.
        '''

        hpx_maps, coord = self.config.load_survey(survey_name)
        hpx_map = combine_maps(hpx_maps)

        if label is None:
            label = survey_name

        self.superimpose_hpxmap(hpx_map, label, color=color,
                                coord_in=coord, cbar=cbar)

    def superimpose_survey_outline(self, survey_name, color='red',
                                   label=None, cbar=True):
        '''Superimpose an outline of a survey

        Parameters
        ----------
        survey_name : string
            The name of the survey in the configuration file

        color : string or array-like with shape (3,), optional
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet. Default = 'red'

        label : string, optional
            The name to use when labeling this survey on the footprint.
            If not input, we will use the survey name.

        Notes
        -----
        This function is for survey footprints that are defined in the
        configuration file as opposed loading a healpix map
        '''

        vtxs, coord = self.config.load_survey_outline(survey_name)

        if label is None:
            label = survey_name

        if isinstance(coord, list):
            for vtx1, coord1 in zip(vtxs[:-1], coord[:-1]):
                self.superimpose_polygon_outline(vtx1, label, color=color,
                                                 coord_in=coord1, cbar=False)
            vtxs = vtxs[-1]
            coord = coord[-1]

        self.superimpose_polygon_outline(vtxs, label, color=color,
                                         coord_in=coord, cbar=cbar)

    def superimpose_survey_contour(self, survey_name, color='red',
                                   label=None, frac=0.85, **kwds):
        '''Superimpose an outline of a survey.

        Parameters
        ----------
        survey_name : string
            The name of the survey in the configuration file

        color : string or array-like with shape (3,), optional
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet. Default = 'red'

        label : string, optional
            The name to use when labeling this survey on the footprint.
            If not input, we will use the survey name.

        Notes
        -----
        This function is for entries that load Healpix maps. We draw
        contours instead of plotting an image of the map. For entries that
        are not Healpix maps, but define the survey region, try
        superimpose_survey_outline(...)
        '''

        hpx_maps, coord = self.config.load_survey(survey_name)

        if label is None:
            label = survey_name

        for hpx_map in hpx_maps[:-1]:
            self.superimpose_hpxmap_contour(hpx_map, label, color=color,
                                            coord_in=coord, cbar=False,
                                            **kwds)

        self.superimpose_hpxmap_contour(hpx_maps[-1], label, color=color,
                                        coord_in=coord, cbar=True, frac=frac,
                                        **kwds)

    def superimpose_hpxmap_contour(self, hpx_map, label, color='red',
                                   coord_in='C', cbar=True, frac=0.85,
                                   smooth_map=None):
        '''Superimpose a contour of an input healpix map.

        Parameters
        ----------
        hpx_map : array-like
            The input healpix map

        label : string
            The name to use as a label for the input map

        color : string or array-like with shape (3,), optional
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet. Default = 'red'

        coord_in : 'C', 'E', or 'G', optional
            The coordinate system of the input map. Default = 'C'.

        cbar : boolean, optional
            Whether to add a colorbar labeling the input map. Default = true.

        frac : float, optional
            The contour level will be drawn containing `frac' levels of observation time

        smooth_map : float
            FWHM to smooth the input map (in arcminutes)
        '''

        idx_nan = (hpx_map == 0)

#       Smoothing makes it more likely that contours don't have holes in them,
#       but it takes some time to smooth each map
        if smooth_map:
            hpx_map = H.smoothing(hpx_map, fwhm=np.radians(smooth_map/60.0),
                                  verbose=False)

        hpx_map /= np.max(hpx_map)
        hpx_map[idx_nan] = np.NaN

        cm1 = util.get_color_map(color)

        coord = [coord_in, self.coord_plot]

        level = determine_level(hpx_map, frac)

        if self.partialmap:
#           Colorbar is added to this and then deleted to make sure there is
#           room at the bottom of the map for the labels. Margins are to make
#           sure the title is not partially off the figure for a square map
            sub = (1, 1, 1)
            margins = (0.01, 0.025, 0.01, 0.03)
            map_tmp = H.cartcontour(hpx_map, 5, title='', coord=coord,
                                    fig=self.fig.number, cmap=cm1, notext=True,
                                    flip='astro', latra=self.latra,
                                    lonra=self.lonra, sub=sub, margins=margins,
                                    return_projected_map=True, **self.kwds)
            idx = np.isfinite(map_tmp)
            if cbar:
                cbar = len(map_tmp[idx]) > 0
        else:
            self.mapcontour(hpx_map, [-0.1, level], title='', coord=coord,
                            cbar=True, fig=self.fig.number, cmap=cm1,
                            notext=True, flip='astro', **self.kwds)
        self.fig.delaxes(self.fig.axes[-1])

        if cbar:
    #       Temporary axis with a Healpix map so I can get the correct color
    #       for the colorbar
            cm1 = util.get_color_map(color)
            coord = [coord_in, self.coord_plot]
            hpx_map = np.ones(12*32**2)
            self.mapview(hpx_map, title='', coord=coord,
                         cbar=None, fig=self.fig.number, cmap=cm1,
                         notext=True, flip='astro', **self.kwds)

#           First add the new colorbar axis to the figure
            im0 = self.fig.axes[-1].get_images()[0]
            box = self.fig.axes[0].get_position()
            ax_color = pl.axes([len(self.cbs), box.y0-0.1, 0.05, 0.05])
            self.fig.colorbar(im0, cax=ax_color, orientation='horizontal',
                              label=label, values=[1, 1])

            self.cbs.append(ax_color)

            self.fig.delaxes(self.fig.axes[-2])

#           Readjust the location of every colorbar
            ncb = len(self.cbs)

            left = 1.0 / (2.0*ncb) - 0.025
            for ax_tmp in self.cbs:
                ax_tmp.set_position([left, box.y0-0.1, 0.05, 0.05])
                left += 1.0 / ncb

    def superimpose_polygon_outline(self, vertices, label, color='red',
                                    coord_in='C', cbar=True):
        '''Superimpose an outline of a survey given input vertices

        Parameters
        ----------
        vertices: array-like (nvtxs, 2)
            The vertices of the polygon

        label : string
            The label for the survey

        color : string or array-like with shape (3,)
            The color to use when overlaying the survey footprint. Either a
            string or rgb triplet.

        coord_in : 'C', 'E', or 'G'
            The coordinate system for the input vertices

        cbar : boolean
            Whether to add a colorbar corresponding to this polygon or not
        '''

        lons = vertices[:, 0]
        lons = np.append(lons, lons[0])

        lats = vertices[:, 1]
        lats = np.append(lats, lats[0])

        nvertices = len(lons)

#       Loop over all vertices and generate lines between adjacent vertices
#       in list. This is to ensure the lines are drawn.
        linelon = np.array([])
        linelat = np.array([])
        for i in range(nvertices-1):
            tmplon = np.linspace(lons[i], lons[i+1], num=1000)
            tmplat = np.linspace(lats[i], lats[i+1], num=1000)
            linelon = np.append(linelon, tmplon)
            linelat = np.append(linelat, tmplat)

        H.projplot(linelon, linelat, lonlat=True, markersize=1,
                   color=color)

        if cbar:
#           Temporary axis with a Healpix map so I can get the correct color
#           for the colorbar
            cm1 = util.get_color_map(color)
            coord = [coord_in, self.coord_plot]
            hpx_map = np.ones(12*32**2)
            self.mapview(hpx_map, title='', coord=coord,
                         cbar=True, fig=self.fig.number, cmap=cm1,
                         notext=True, flip='astro', **self.kwds)
            self.fig.delaxes(self.fig.axes[-1])

#           First add the new colorbar axis to the figure
            im0 = self.fig.axes[-1].get_images()[0]
            box = self.fig.axes[0].get_position()
            ax_color = pl.axes([len(self.cbs), box.y0-0.1, 0.05, 0.05])
            self.fig.colorbar(im0, cax=ax_color, orientation='horizontal',
                              label=label, values=[1, 1])

            self.cbs.append(ax_color)

#           Delete the temporary map
            self.fig.delaxes(self.fig.axes[-2])

#           Readjust the location of every colorbar
            ncb = len(self.cbs)

            left = 1.0 / (2.0*ncb) - 0.025
            for ax_tmp in self.cbs:
                ax_tmp.set_position([left, box.y0-0.1, 0.05, 0.05])
                left += 1.0 / ncb

def combine_maps(hpx_maps):
    '''Code to combine an array of maps.

    Parameters
    ----------
    hpx_maps : list
        A list of healpix maps we want to combine into a single map

    Notes
    -----
    This is called when we plot images of the surveys. It is not called
    when we plot outlines
    '''

    map_comb = np.zeros_like(hpx_maps[0])

    for hpx_map in hpx_maps:
        map_comb += hpx_map

    map_comb /= np.max(map_comb)

    return map_comb

def determine_level(hpx_map, obs_frac, time=True):
    '''Determine the contour level than contains the obs_frac of the total
    observation time.

    Parameters
    ----------
    hpx_map : array-like
        The input healpix map. This should be a survey footprint of a
        single patch.

    obs_frac: float
        The fraction of observation time that we want the contour to enclose.

    time : boolean, optional
        If true, obs_frac is fraction of observation time. If false,
        obs_frac is fraction of observation area. Default = True
    '''

    idx = hpx_map > 0

    vals = hpx_map[idx]

    nvals = len(vals)

    vals_sort = np.sort(vals)

    if time is False:
        level = vals_sort[int((1-obs_frac)*(nvals-1))]
        return level

    idx0 = 0
    idx1 = nvals-1

    totvals = np.sum(vals_sort)

    while 1:
        idx2 = int(0.5*(idx0+idx1))
        vals_sum = np.sum(vals_sort[idx2:])
        vals_frac = vals_sum / totvals

        if (idx2-idx0) < 5:
            return vals_sort[idx2]
        elif np.abs(vals_frac - obs_frac) < 0.01:
            return vals_sort[idx2]

        if vals_frac < obs_frac:
            idx1 = idx2
        else:
            idx0 = idx2

