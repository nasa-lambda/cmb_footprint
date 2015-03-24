#!/usr/bin/env python

import numpy as np
import pylab as pl
import healpy as H

import footprint

if __name__ == '__main__':

#   fn_background = 'maps/Planck/HFI_SkyMap_353_2048_R2.00_full.fits'
    fn_background = 'maps/Planck/COM_CompMap_DustPol-commander_1024_R2.00.fits'

    data = H.read_map(fn_background,field=(0,1))
    background_map = np.sqrt(data[0]**2 + data[1]**2)

    fp = footprint.SurveyStack(background_map, fignum=1, projection=H.cartview, coord_plot='C')

#   fn_148_south_s2 = 'maps/ACT/ACT_148_south_season_2_1way_hits_v3.fits'
    fn_148_equ_s3 = 'maps/ACT/ACT_148_equ_season_3_1way_hits_v3.fits'
#   fn_148_south_s3 = 'maps/ACT/ACT_148_south_season_3_1way_hits_v3.fits'
#   fn_148_equ_s4 = 'maps/ACT/ACT_148_equ_season_4_1way_hits_v3.fits'
#   fn_148_south_s4 = 'maps/ACT/ACT_148_south_season_4_1way_hits_v3.fits'

#   fn_220_equ_s4 = 'maps/ACT/ACT_220_equ_season_4_1way_hits_v3.fits'
#   fn_220_south_s4 = 'maps/ACT/ACT_220_south_season_4_1way_hits_v3.fits'

    fn_148_spt = 'maps/SPT/spt_wgt_ra5h30dec-55_2008_150ghz_zea_dr1.fits'

#   fns_ACT = [fn_148_equ_s4, fn_148_south_s4, fn_148_equ_s3, fn_148_south_s3,
#              fn_148_south_s2]
    fns_ACT = [fn_148_equ_s3]
    fns_SPT = [fn_148_spt]

#   fp.superimpose_fits(fns_SPT,color='red',maptype='WCS')
#   footprint.superimpose_boundary_cen((0,0),(40,5),color='blue')
    fp.superimpose_experiment('ACT', color='green')
    fp.superimpose_experiment('SPT', color='red')

    pl.show()
