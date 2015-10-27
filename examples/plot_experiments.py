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

    #QUIET
    fp = footprint.SurveyStack('PLANCK-DUSTPOL', fignum=1, projection='mollweide', coord_plot='C', rot=[0,0],
                               config='../footprint.cfg', title='QUIET Survey Areas')
    fp.superimpose_survey('QUIET-Q-CMB',color='blue', label='CMB')
    fp.superimpose_survey('QUIET-Q-GAL',color='red', label='Galaxy')
    pl.savefig('QUIET_footprint.png')

    #POLARBEAR
    fp = footprint.SurveyStack('PLANCK-DUSTPOL', fignum=2, projection='mollweide', coord_plot='C', rot=[0,0],
                               config='../footprint.cfg', title='POLARBEAR Survey Areas')
    fp.superimpose_survey('POLARBEAR',color='blue', add_cb=False)
    pl.savefig('POLARBEAR_footprint.png')

    #BICEP
    fp = footprint.SurveyStack('PLANCK-DUSTPOL', fignum=3, projection='mollweide', coord_plot='C', rot=[0,0],
                               config='../footprint.cfg', title='BICEP Survey Area')
    fp.superimpose_survey('BICEP2',color='blue', add_cb=False)
    pl.savefig('BICEP2_footprint.png')

    #SPIDER
    fp = footprint.SurveyStack('PLANCK-DUSTPOL', fignum=4, projection='mollweide', coord_plot='C', rot=[0,0],
                               config='../footprint.cfg', title='SPIDER Survey Area')
    fp.superimpose_survey('SPIDER-90',color='blue', add_cb=False)
    pl.savefig('SPIDER_footprint.png')
    
    #SPT
    fp = footprint.SurveyStack('PLANCK-DUST', fignum=5, projection='mollweide', coord_plot='C', rot=[0,0],
                               config='../footprint.cfg', title='SPT-SZ Survey Areas', max=5.0e5)
    fp.superimpose_survey_outline('SPT-SZ',color='blue', add_cb=False)
    pl.savefig('SPT_footprint.png')
    
    #ACT
    fp = footprint.SurveyStack('PLANCK-DUST', fignum=6, projection='mollweide', coord_plot='C', rot=[0,0],
                               config='../footprint.cfg', title='ACT Survey Areas', max=5.0e5)
    fp.superimpose_survey('ACT-EQU',color='blue', label='Equatorial')
    fp.superimpose_survey('ACT-SOUTH',color='red', label='Southern')
    pl.savefig('ACT_footprint.png')

