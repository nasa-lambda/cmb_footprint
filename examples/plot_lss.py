#!/usr/bin/env python

import numpy as np
import pylab as pl
import healpy as H

# append the cmb_footprint source path
# you can add this to your PYTHONPATH also
import os
import sys
separator = os.sep
fullpath = os.getcwd().split(separator)
source_path = separator.join(fullpath[0:-2])
sys.path.append(source_path)
from cmb_footprint import footprint

if __name__ == '__main__':

#   fn_background = '../maps/Planck/COM_CompMap_DustPol-commander_1024_R2.00.fits'
#   data = H.read_map(fn_background,field=(0,1))
#   background_map = np.sqrt(data[0]**2 + data[1]**2)
#   fp = footprint.SurveyStack(background_map, projection=H.mollview,
#                              coord_plot='C')

    fp = footprint.SurveyStack('PLANCK-DUSTPOL', projection=H.mollview, coord_plot='C', rot=[0,0],
                               config='../footprint.cfg')

    fp.superimpose_survey('BOSS-LSS-RANDOM2',color='green')
    fp.superimpose_survey('BOSS-QSO',color='red')

    fp.superimpose_survey('HSC-FALL1-EQ',color='blue')
    fp.superimpose_survey('HSC-FALL2-EQ',color='orange')

    pl.show()
