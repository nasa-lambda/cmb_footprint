import numpy as np
import pylab as pl
import healpy as H
import matplotlib.cm as cm
import pyfits
from astropy import wcs
import scipy.sparse

def plot_background(fn,fignum=None):

	maps = H.read_map(fn,field=(1,2))

	polamp = np.sqrt(maps[0]**2 + maps[1]**2)

	map2d = project_healpix(polamp,coord=['G','C'])

	pl.figure(fignum)
	pl.imshow(map2d,cmap=cm.Greys,vmin=0.0,vmax=0.001)

	#pl.show()

def project_healpix(mapinp,coord=None):

	map2d = H.mollview(mapinp,xsize=1600,coord=coord,return_projected_map=True)
	pl.close()

	return map2d

def plot_ACT(fignum=None):
	
	#ACT data are stored in FITS files in WCS format. 

	fn_148_south_s2 = 'maps/ACT/ACT_148_south_season_2_1way_hits_v3.fits'
	fn_148_equ_s3 = 'maps/ACT/ACT_148_equ_season_3_1way_hits_v3.fits'
	fn_148_south_s3 = 'maps/ACT/ACT_148_south_season_3_1way_hits_v3.fits'
	fn_148_equ_s4 = 'maps/ACT/ACT_148_equ_season_4_1way_hits_v3.fits'
	fn_148_south_s4 = 'maps/ACT/ACT_148_south_season_4_1way_hits_v3.fits'

	fn_220_equ_s4 = 'maps/ACT/ACT_220_equ_season_4_1way_hits_v3.fits'
	fn_220_south_s4 = 'maps/ACT/ACT_220_south_season_4_1way_hits_v3.fits'

	fns = [fn_148_equ_s4,fn_148_south_s4,fn_148_equ_s3,fn_148_south_s3,fn_148_south_s2]

	nside = 2048
	npix = H.nside2npix(nside)

	hpx_map = np.zeros(npix)
	for fn in fns:
		hdulist = pyfits.open(fn)
		hpx_map += wcs_to_healpix(hdulist,nside)
		hdulist.close()

	hpx_148 = hpx_map

	hpx_148_2d = project_healpix(hpx_148)

	#This is so these points are not plotted. NaNs are ignored in imshow.
	idx = np.where(hpx_148_2d == 0)
	hpx_148_2d[idx] = np.NaN

	pl.figure(fignum)
	pl.imshow(hpx_148_2d,alpha=0.5)

	#H.mollview(hpx_148_s4)
	#pl.show()

def wcs_to_healpix(hdulist,nside):
	'''Converts data in an opened FITS file from a generic WCS to Healpix'''
	
	#use wcs_pix2world to get ra/dec coordinates, then use ang2pix to get Healpix coordinates and 
	#sum over all pixels with the same Healpix pixel numer

	w = wcs.WCS(hdulist[0].header)

	data = hdulist[0].data

	xsize, ysize = data.shape
	print "x/ysize = ", xsize, ysize

	pixcrd = [(x,y) for x in xrange(xsize) for y in xrange(ysize)]

	print "np.array"
	pixcrd2 = np.array(pixcrd)

	#pix2world to get ra/dec values of each pixel
	print "wcs_pix2world"
	world = w.wcs_pix2world(pixcrd,1)

	idx = np.isfinite(world[:,0])

	pixcrd2 = pixcrd2[idx,:]
	world = world[idx,:]

	#Use Healpy to do ang2pix for Healpix pixel numbers
	theta = np.radians(world[:,0])
	phi = np.radians(world[:,1]) #need to check sign of this
	print "ang2pix"
	pixnum = H.ang2pix(nside,theta,phi)
	

	print "hpx_data"
	npix = H.nside2npix(nside)
	row = pixnum
	col = np.zeros_like(row)

	#coo_matrix will sum entries that have multiple elements
	hpx_data = scipy.sparse.coo_matrix((data[(pixcrd2[:,0],pixcrd2[:,1])], (row, col)), shape=np.array([npix,1])).A.transpose()
	hpx_data.shape = (npix)

	return hpx_data

if __name__=='__main__':
	
	fn = 'maps/Planck/HFI_SkyMap_353_2048_R2.00_full.fits'
	plot_background(fn,fignum=0)

	plot_ACT(fignum=0)

	pl.show()
