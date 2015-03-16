import numpy as np
import healpy as H
import astropy.wcs as wcs
import scipy.sparse

def wcs_to_healpix(hdulist,nside):
	'''Converts data in an opened FITS file from a generic WCS to Healpix'''
	
	#use wcs_pix2world to get ra/dec coordinates, then use ang2pix to get Healpix coordinates and 
	#sum over all pixels with the same Healpix pixel numer

	w = wcs.WCS(hdulist[0].header)

	data = hdulist[0].data

	xsize, ysize = data.shape
	print("x/ysize = ", xsize, ysize)

	pixcrd = [(y,x) for x in xrange(xsize) for y in xrange(ysize)]

	print("np.array")
	pixcrd2 = np.array(pixcrd)

	#pix2world to get ra/dec values of each pixel
	print("wcs_pix2world")
	world = w.wcs_pix2world(pixcrd,1)

	#needed since for some reason I was getting NaNs for some ACT data. Probably because of my swapping of x and y
	idx = np.isfinite(world[:,0])

	pixcrd2 = pixcrd2[idx,:]
	world = world[idx,:]

	#Use Healpy to do ang2pix for Healpix pixel numbers
	theta = np.pi/2.0 - np.radians(world[:,1])
	phi = np.radians(world[:,0]) #Sign seems to be correct based on looking at ACT hits map.
	print("ang2pix")
	pixnum = H.ang2pix(nside,theta,phi)
	

	print("hpx_data")
	npix = H.nside2npix(nside)
	row = pixnum
	col = np.zeros_like(row)

	#coo_matrix will sum entries that have multiple elements
	hpx_data = scipy.sparse.coo_matrix((data[(pixcrd2[:,1],pixcrd2[:,0])], (row, col)), shape=np.array([npix,1])).A.transpose()
	hpx_data.shape = (npix)

	return hpx_data

