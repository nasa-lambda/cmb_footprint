# pylint: disable=E1101
# pylint: disable=C0325

'''
====================================================================
config_handler.py : Configuration file related classes and functions
====================================================================

This module provides the class that interprets the configuration file and
generates or reads in healpix maps corresponding to survey footprints.
'''

from __future__ import print_function

# For differing imports between Python2 and Python3
try:
    import ConfigParser
    from urllib2 import urlopen
except ImportError:
    import configparser as ConfigParser
    from urllib.request import urlopen

import os
import hashlib

import healpy as H
import numpy as np
from astropy.coordinates import SkyCoord

import cmb_footprint.util as util

class ConfigHandler(object):
    '''Class that handles reading the configuration file and generating the
    footprints described by the configuration file'''

    def __init__(self, config_fn, map_path, nside=256, download_config=False):
        self.config_fn = config_fn
        self.map_path = map_path

        if nside is None:
            nside = 256

        self.nside = nside

        if download_config:
            self.get_config()

        self.config = ConfigParser.ConfigParser()
        self.config.read(config_fn)
        self.check_all()

    def get_config(self):
        '''Download the configuration file from LAMBDA to the filename given
        in self.config_fn in the local directory.

        Notes
        -----
        Remote file is footprint.txt instead of footprint.cfg because server
        gives an error if you try to download a .cfg file whereas it allows
        you to download a .txt file.
        '''

        url_pre = 'http://lambda.gsfc.nasa.gov/data/footprint-maps/'

        local_path = os.path.join('./', self.config_fn)
        url = os.path.join(url_pre, 'footprint.txt')

        print("Downloading configuration file")

        file_chunk = 16 * 1024
        req = urlopen(url)
        with open(local_path, 'wb') as fp1:
            while True:
                chunk = req.read(file_chunk)
                if not chunk:
                    break
                fp1.write(chunk)

    def check_all(self):
        '''Check if we have all the correct hitmaps downloaded locally
        for surveys in which we provide Healpix files'''

        for section in self.config.sections():
            handler = self.config.get(section, 'handler')
            if handler == 'hpx_file':
                self._download_files(section)

    def load_survey(self, survey_name):
        '''Load or generate the healpix map for a survey defined in the
        configuration file

        Parameters
        ----------
        survey_name : string
            A survey defined in the configuration file

        Returns
        -------
        hpx_map : array-like
            Healpix map with the hitmap for the requested survey
        '''

        handler = self.config.get(survey_name, 'handler')

        func = getattr(self, 'get_'+handler)
        hpx_map, coord = func(survey_name)

        return hpx_map, coord

    def _download_files(self, survey_name):
        '''Return the local paths to the files associated with a survey.
        If the files do not exists or their MD5 checksums do not match what is
        stored in the configuration file, download new files.

        Parameters
        ----------
        survey_name : string
            The survey for which we want the filenames of the hitmaps.

        Returns
        -------
        fns_out : string
            The local filenames

        '''

        urls = self.config.get(survey_name, 'url')
        cksums = self.config.get(survey_name, 'checksum')

        urls = urls.split(',')
        cksums = cksums.split(',')

        fns_out = []
#       Compare the checksum of the local file to the value in the
#       configuration file. If they don't match (or local file does not
#       exist), download the file. If the checksum of the downloaded file does
#       not match the checksum in the configuration file something is wrong
        for url, cksum_cfg in zip(urls, cksums):
            fn_tmp = os.path.split(url)[1]
            download_file = False
            local_path = os.path.join(self.map_path, fn_tmp)
            if os.path.exists(local_path):
                cksum_file = hashlib.md5(open(local_path, 'rb').read()).hexdigest()
                if cksum_cfg != cksum_file:
                    download_file = True
            else:
                download_file = True

            if download_file:
                print("Downloading map for", survey_name)

                if not util.download_url(url, cksum_cfg, local_path):
                    raise ValueError('Could not download file')

            fns_out.append(local_path)

        return fns_out

    def get_hpx_file(self, survey_name):
        '''Load the healpix file associated with a given survey.

        Handler in .cfg file: "hpx_file"

        Parameters
        ----------
        survey_name : string
            Configuration file block specifying survey parameters

        Returns
        -------
        hpx_map : array-like
            The healpix map associated with the survey
        '''

        fns = self._download_files(survey_name)
        coord = self.config.get(survey_name, 'coord')

        hpx_map = H.read_map(fns[0], verbose=False)
        nside = H.npix2nside(len(hpx_map))
        for fn_tmp in fns[1:]:
            tmp_map = H.read_map(fn_tmp, verbose=False)
            hpx_map += H.ud_grade(tmp_map, nside)

        return hpx_map, coord

    def get_polygon(self, survey_name):
        '''Generates a healpix map for a given survey given vertices of
        a polygon.

        Handler in .cfg file: "polygon"

        Parameters
        ----------
        survey_name : string
            Configuration file block specifying survey parameters

        Returns
        -------
        hpx_map : array-like
            The healpix map associated with the survey

        Notes
        -----
        The configuration file must contain 'vertex1', ..., 'vertexXX'
        entries containing the ra,dec locations of each vertex. Nside of
        the output map is the same as the nside specified in the
        initialization of this class.
        '''

        lons = []
        lats = []

        i = 1
        while True:
            radec_point = 'vertex'+str(i)
            try:
                radec_val = self.config.get(survey_name, radec_point)
            except ConfigParser.NoOptionError:
                break

            lonlat_val = radec_val.split(',')
            lonlat_val = [tmp.strip() for tmp in lonlat_val]
            lonlat_val = SkyCoord(lonlat_val[0], lonlat_val[1])
            lons.append(lonlat_val.ra.deg)
            lats.append(lonlat_val.dec.deg)
            i += 1

        vtxs = np.transpose([lons, lats])

        coord = self.config.get(survey_name, 'coord')

        hpx_map = util.gen_map_polygon(vtxs, self.nside)

        return hpx_map, coord

    def get_disc(self, survey_name):
        '''Generates a healpix map for a given survey given a center point
        and a disc radius.

        Handler in .cfg file: "disc"

        Parameters
        ----------
        survey_name : string
            Configuration file block specifying survey parameters

        Returns
        -------
        hpx_map : array-like
            The healpix map associated with the survey

        Notes
        -----
        The configuration file must contain a 'center' and a 'radius'
        entry. The nside of the output map is the same as the nside specified
        in the initializaiton of this class.
        '''

        center = self.config.get(survey_name, 'center').split(',')
        center = [tmp.strip() for tmp in center]
        center = SkyCoord(center[0], center[1])
        center = [center.ra.deg, center.dec.deg]

#       This calculation is so that radius can be input in different
#       coordinates (i.e. deg, arcminutes, etc.)
        radius = self.config.get(survey_name, 'radius')
        tmp = SkyCoord('0d', radius)
        radius = np.abs(tmp.dec.deg)

        coord = self.config.get(survey_name, 'coord')

        hpx_map = util.gen_map_disc(center, radius, self.nside)

        return hpx_map, coord

    def get_rectangle(self, survey_name):
        '''Generates a healpix map for a given survey given a center point
        and edge length of the rectanlge.

        Handler in .cfg file: "rectangle"

        Parameters
        ----------
        survey_name : string
            Configuration file block specifying survey parameters

        Returns
        -------
        hpx_map : array-like
            The healpix map associated with the survey

        Notes
        -----
        The configuration file must contain a 'center', a 'size', and a
        'coord' entry. The 'size' entry must contain the lon,lat length of the
        edge of the rectangle. The nside of the output map is the same as the
        nside specified in the initialization of this class.
        '''

        center = self.config.get(survey_name, 'center').split(',')
        center = [tmp.strip() for tmp in center]
        center = SkyCoord(center[0], center[1])
        center = [center.ra.deg, center.dec.deg]

#       edge length. The corners of the box are center +- length/2.0
        length = self.config.get(survey_name, 'size').split(',')
        length = [tmp.strip() for tmp in length]
        length = SkyCoord(length[0], length[1])
        length = [length.ra.deg, length.dec.deg]

        coord = self.config.get(survey_name, 'coord')

        hpx_map = util.gen_map_centersize(center, length, self.nside)

        return hpx_map, coord

    def get_rectangle_bounds(self, survey_name):
        '''Generates a healpix map for a given survey given a range of
        RA/Dec.

        Handler in .cfg file: "rectangle_bounds"

        Parameters
        ----------
        survey_name : string
            Configuration file block specifying survey parameters

        Returns
        -------
        hpx_map : array-like
            The healpix map associated with the survey

        Notes
        -----
        The configuration file must contain an 'ra_range' and 'dec_range'.
        The nside of the output map is the same as the nside specified in the
        initialization of this class.
        '''

        ra_range = self.config.get(survey_name, 'ra_range').split(',')
        ra_range = [tmp.strip() for tmp in ra_range]

        dec_range = self.config.get(survey_name, 'dec_range').split(',')
        dec_range = [tmp.strip() for tmp in dec_range]

        bottom_left = SkyCoord(ra_range[0], dec_range[0])
        top_right = SkyCoord(ra_range[1], dec_range[1])

#       corner1 = (top_right.ra.deg, top_right.dec.deg)
#       corner2 = (top_right.ra.deg, bottom_left.dec.deg)
#       corner3 = (bottom_left.ra.deg, bottom_left.dec.deg)
#       corner4 = (bottom_left.ra.deg, top_right.dec.deg)
#       vertices = (corner1, corner2, corner3, corner4)

        lonra = [bottom_left.ra.deg, top_right.ra.deg]
        latra = [bottom_left.dec.deg, top_right.dec.deg]

#       hpx_map = util.gen_map_polygon(vertices, self.nside)
        hpx_map = util.gen_map_rectangle(lonra, latra, self.nside)

        coord = self.config.get(survey_name, 'coord')

        return hpx_map, coord

    def get_combination(self, survey_name):
        '''Generates a healpix map for a given survey that is a
        combination of multiple other surveys in the configuration
        file.

        Handler in .cfg file: "combination"

        Parameters
        ----------
        survey_name : string
            Configuration file block specifying survey parameters

        Returns
        -------
        hpx_map : array-like
            The healpix map associated with the survey

        Notes
        -----
        The configuration file must contain a 'components' entry which lists
        the different configuration file entries, separated by a ',' that make
        up this survey.
        '''

        components = self.config.get(survey_name, 'components').split(',')
        components = [tmp.strip() for tmp in components]

        hpx_map, coord = self.load_survey(components[0])
        for component in components[1:]:
            tmp_map, coord_tmp = self.load_survey(component)
            if coord_tmp != coord:
#               Could we do something to rotate the maps using Healpy?
                raise ValueError('All maps in combination must be in the same coordinate system')
            hpx_map += tmp_map

        return hpx_map, coord
