import numpy as np
import pylab as pl
import healpy as H
import matplotlib

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
    fp = footprint.SurveyStack('PLANCK-DUSTPOL', fignum=1, projection='mollweide', coord_plot='C', rot=[0,0],
                               config='../footprint.cfg', cbar=True, log=True)
    fp.superimpose_survey('SPIDER-90',color='orange', cbar=False)
    fp.superimpose_survey('QUIET-Q',color='green', label='QUIET', cbar=False)
    pl.show()
    
    
