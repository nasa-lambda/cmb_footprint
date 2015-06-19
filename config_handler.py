# pylint: disable=E1101
# pylint: disable=C0325

'''
====================================================================
config_handler.py : Configuration file related classes and functions
====================================================================

This module provides the class that interprets the configuration file and
generates or reads in healpix maps corresponding to the survey footprints
of each experiment
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
        for experiments in which we provide Healpix files'''

        for section in self.config.sections():
            handler = self.config.get(section, 'handler')
            if handler == 'hpx_file':
                self.download_files(section)

    def load_experiment(self, experiment_name):
        '''Load or generate the healpix map for an experiment defined in the
        configuration file

        Parameters
        ----------
        experiment_name : string
            An experiment defined in the configuration file

        Returns
        -------
        hpx_map : array-like
            Healpix map with the hitmap for the requested experiment
        '''

        handler = self.config.get(experiment_name, 'handler')

        func = getattr(self, 'get_'+handler)
        hpx_map, coord = func(experiment_name)

        return hpx_map, coord

    def download_files(self, experiment_name):
        '''Return the local paths to the files associated with an experiment.
        If the files do not exists or their MD5 checksums do not match what is
        stored in the configuration file, download new files.

        Parameters
        ----------
        experiment_name : string
            The experiment for which we want the filenames of the hitmaps. It
            must be listed as 'hpx_file' in the configuration file

        Returns
        -------
        fns_out : string
            The local filenames

        '''

        urls = self.config.get(experiment_name, 'url')
        cksums = self.config.get(experiment_name, 'checksum')

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
                print("Downloading map for", experiment_name)

                if not util.download_url(url, cksum_cfg, local_path):
                    raise ValueError('Could not download file')

            fns_out.append(local_path)

        return fns_out

    def get_hpx_file(self, experiment_name):
        '''Load the healpix file associated with a given experiment. Must be
        listed as 'hpx_file' in the configuration file.

        Parameters
        ----------
        experiment_name : string
            The experiment for which we want the hitmaps. It
            must be listed as 'hpx_file' in the configuration file.

        Returns
        -------
        hpx_map : array-like
            The healpix map associated with the experiment
        '''

        fns = self.download_files(experiment_name)
        coord = self.config.get(experiment_name, 'coord')

        hpx_map = H.read_map(fns[0], verbose=False)
        nside = H.npix2nside(len(hpx_map))
        for fn_tmp in fns[1:]:
            tmp_map = H.read_map(fn_tmp, verbose=False)
            hpx_map += H.ud_grade(tmp_map, nside)

        return hpx_map, coord

    def get_polygon(self, experiment_name):
        '''Generates a healpix map for a given experiment given vertices of
        a polygon. Must be listed as 'radec_polygon' in the configuration
        file.

        Parameters
        ----------
        experiment_name : string
            The experiment for which we want the hitmaps. It must be listed
            as 'radec_polygon' in the configuration file.

        Returns
        -------
        hpx_map : array-like
            The healpix map associated with the experiment

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
                radec_val = self.config.get(experiment_name, radec_point)
            except ConfigParser.NoOptionError:
                break

            lonlat_val = radec_val.split(',')
            lonlat_val = [tmp.strip() for tmp in lonlat_val]
            lonlat_val = SkyCoord(lonlat_val[0], lonlat_val[1])
            lons.append(lonlat_val.ra.deg)
            lats.append(lonlat_val.dec.deg)
            i += 1

        vtxs = np.transpose([lons, lats])

        coord = self.config.get(experiment_name, 'coord')

        hpx_map = util.gen_map_polygon(vtxs, self.nside)

        return hpx_map, coord

    def get_disc(self, experiment_name):
        '''Generates a healpix map for a given experiment given a center point
        and a disc radius. Must be listed as 'radec_disc' in the configuration
        file.

        Parameters
        ----------
        experiment_name : string
            The experiment for which we want the hitmaps. It must be listed
            as 'radec_disc' in the configuration file.

        Returns
        -------
        hpx_map : array-like
            The healpix map associated with the experiment

        Notes
        -----
        The configuration file must contain a 'radec_cen' and a 'radius'
        entry. The nside of the output map is the same as the nside specified
        in the initializaiton of this class.
        '''

        center = self.config.get(experiment_name, 'center').split(',')
        center = [tmp.strip() for tmp in center]
        center = SkyCoord(center[0], center[1])
        center = [center.ra.deg, center.dec.deg]

#       This calculation is so that radius can be input in different
#       coordinates (i.e. deg, arcminutes, etc.)
        radius = self.config.get(experiment_name, 'radius')
        tmp = SkyCoord('0d', radius)
        radius = np.abs(tmp.dec.deg)

        coord = self.config.get(experiment_name, 'coord')

        hpx_map = util.gen_map_disc(center, radius, self.nside)

        return hpx_map, coord

    def get_rectangle(self, experiment_name):
        '''Generates a healpix map for a given experiment given a center point
        and edge length of the rectanlge. Must be listed as 'radec_rect' in
        the configuration file.

        Parameters
        ----------
        experiment_name : string
            The experiment for which we want the hitmaps. It must be listed
            as 'radec_rect' in the configuration file.

        Returns
        -------
        hpx_map : array-like
            The healpix map associated with the experiment

        Notes
        -----
        The configuration file must contain a 'center', a 'size', and a
        'coord' entry. The 'size' entry must contain the lon,lat length of the
        edge of the rectangle. The nside of the output map is the same as the
        nside specified in the initialization of this class.
        '''

        center = self.config.get(experiment_name, 'center').split(',')
        center = [tmp.strip() for tmp in center]
        center = SkyCoord(center[0], center[1])
        center = [center.ra.deg, center.dec.deg]

#       edge length. The corners of the box are center +- length/2.0
        length = self.config.get(experiment_name, 'size').split(',')
        length = [tmp.strip() for tmp in length]
        length = SkyCoord(length[0], length[1])
        length = [length.ra.deg, length.dec.deg]

        coord = self.config.get(experiment_name, 'coord')

        hpx_map = util.gen_map_centersize(center, length, self.nside)

        return hpx_map, coord

    def get_combination(self, experiment_name):
        '''Generates a healpix map for a given experiment that is a
        combination of multiple other experiments in the configuration
        file.

        Parameters
        ----------
        experiment_name : string
            The experiment for which we want the hitmaps. It must be listed
            as 'combination' in the configuration file.

        Returns
        -------
        hpx_map : array-like
            The healpix map associated with the experiment

        Notes
        -----
        The configuration file must contain a 'components' entry which lists
        the different configuration file entries, separated by a ',' that make
        up this experiment.
        '''

        components = self.config.get(experiment_name, 'components').split(',')
        components = [tmp.strip() for tmp in components]

        hpx_map, coord = self.load_experiment(components[0])
        for component in components[1:]:
            tmp_map, coord_tmp = self.load_experiment(component)
            if coord_tmp != coord:
#               Could we do something to rotate the maps using Healpy?
                raise ValueError('All maps in combination must be in the same coordinate system')
            hpx_map += tmp_map

        return hpx_map, coord
