import footprint

if __name__ == '__main__':

    fn_background = 'maps/Planck/HFI_SkyMap_353_2048_R2.00_full.fits'

    fp = footprint.SurveyStack(fn_background, fignum=1)

#   fn_148_south_s2 = 'maps/ACT/ACT_148_south_season_2_1way_hits_v3.fits'
#   fn_148_equ_s3 = 'maps/ACT/ACT_148_equ_season_3_1way_hits_v3.fits'
#   fn_148_south_s3 = 'maps/ACT/ACT_148_south_season_3_1way_hits_v3.fits'
#   fn_148_equ_s4 = 'maps/ACT/ACT_148_equ_season_4_1way_hits_v3.fits'
#   fn_148_south_s4 = 'maps/ACT/ACT_148_south_season_4_1way_hits_v3.fits'

#   fn_220_equ_s4 = 'maps/ACT/ACT_220_equ_season_4_1way_hits_v3.fits'
#   fn_220_south_s4 = 'maps/ACT/ACT_220_south_season_4_1way_hits_v3.fits'

#   fns_ACT = [fn_148_equ_s4, fn_148_south_s4, fn_148_equ_s3, fn_148_south_s3,
#              fn_148_south_s2]
#   fns_ACT = [fn_148_equ_s4]

#   footprint.superimpose_fits(fns_ACT,color='red',maptype='WCS')
#   footprint.superimpose_boundary_cen((0,0),(40,5),color='blue')
    fp.superimpose_experiment('ACT', color='red')

    pl.show()
