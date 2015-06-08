#!/usr/bin/env python

import numpy as np
import pylab as pl
import healpy as H

from cmb_footprint import footprint

if __name__ == '__main__':

#   fn_background = '../maps/Planck/COM_CompMap_DustPol-commander_1024_R2.00.fits'

#   data = H.read_map(fn_background,field=(0,1))
#   background_map = np.sqrt(data[0]**2 + data[1]**2)

#   fp = footprint.SurveyStack(background_map, fignum=1, projection=H.cartview, coord_plot='C', rot=[0,-32.5], 
#           partialmap=True, latra=[-20,20], lonra=[-20,20], config='../footprint.cfg')
#   fp = footprint.SurveyStack(background_map, fignum=1, projection=H.mollview, coord_plot='C', rot=[0,0],
#                              config='../footprint.cfg') 
    fp = footprint.SurveyStack('PLANCK-DUSTPOL-BACKGROUND', fignum=1, projection=H.mollview, coord_plot='C', rot=[0,0],
                               config='../footprint.cfg') 

    fp.superimpose_experiment('ACT',color='green')
    fp.superimpose_experiment('SPT',color='red')
    fp.superimpose_experiment('POLARBEAR',color='magenta')
    fp.superimpose_experiment('BICEP2',color='blue')
    fp.superimpose_experiment('QUIET',color='orange')

    pl.show()
