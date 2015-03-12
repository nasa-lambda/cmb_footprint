import numpy as np
import pylab as pl
import healpy as H
import matplotlib.cm as cm
import pyfits
from astropy import wcs
import scipy.sparse

class SurveyStack(object):

	def __init__(self,background,xsize=800,nside=None,fignum=None):
		self.xsize = xsize
		self.fn_background = background
		self.fig = pl.figure(fignum)

		maps = H.read_map(background,field=(1,2)) #Q,U should be fields 1 and 2. Background should be a Healpix map

		if nside is None:
			nside = H.npix2nside(len(maps[0]))

		self.nside = nside

		polamp = np.sqrt(maps[0]**2 + maps[1]**2)

		cm.Greys.set_under(alpha=0.0)
		H.mollview(polamp,title='Experiment Footprints',xsize=1600,coord=['G','C'],fig=self.fig.number,cmap=cm.Greys,max=0.001,min=0.0,notext=True,cbar=None)
		H.graticule(dpar=15.0,dmer=15.0,coord='C')

	def read_hpx_maps(self,fns):
		
		hpx_map = np.zeros(H.nside2npix(self.nside))
		for fn in fns:
			tmp_map = H.read_map(fn)
			nside = H.npix2nside(len(tmp_map))
			hpx_map += H.udgrade(tmp_map,self.nside)
	
		return hpx_map

	def read_wcs_maps(self,fns):
		
		hpx_map = np.zeros(H.nside2npix(self.nside))
		for fn in fns:
			hdulist = pyfits.open(fn)
			hpx_map += wcs_to_healpix(hdulist,self.nside)
			hdulist.close()

		return hpx_map

	def superimpose_hpxmap(self,hpx_map,color='red'):
		
		idx_nan = (hpx_map == 0)
		idx_good = (hpx_map != 0)

		hpx_map[idx_nan] = np.NaN

		if color == 'red':
			color_val = 255.0
		elif color == 'blue':
			color_val = 0.0

		#hpx_map[idx_good] = color_val

		#Can't use mollview because it doesn't allow alpha so must create our own HpxMollweideAxes object and add it to the figure. 
		#H.mollview(hpx_map,title='',xsize=1600,coord='C',fig=self.fig.number,min=0.0,max=255.0)
		extent = (0.02,0.05,0.96,0.9) #this is the extent used in H.mollview(). Needed so maps are overlaid correctly
		#cm.jet.set_bad(alpha=0.0) #so parts of the Nhits map that are NaN are completely transparent
		#cm.jet.set_under(alpha=0.0) #so parts in the rectangle outside the mollweide oval are completely transparent
		cm1 = get_color_map(color)

		ax = H.projaxes.HpxMollweideAxes(self.fig,extent,coord='C',rot=None,format='%g',flipconv='astro')
		map2d = ax.projmap(hpx_map,xsize=1600,coord='C',cmap=cm1)
		#map2d = ax.projmap(hpx_map,xsize=1600,coord='C',vmin=0.0,vmax=255.0,alpha=0.5,cmap=cm1)
		self.fig.add_axes(ax)

	def superimpose_fits(self,fns,color='red',maptype='WCS'):
		'''Superimpose the footprint of an experiment on the background image. Can be a single fits file or a list of them that will 
		be added together.'''

		if maptype == 'WCS':
			hpx_map = self.read_wcs_maps(fns)
		elif maptype == 'HPX':
			hpx_map = self.read_hpx_maps(fns)

		self.superimpose_hpxmap(hpx_map,color=color)

	def superimpose_boundary_cen(self,radec_cen,radec_size,color='red'):
		'''Superimpose the footprint of an experiment on the background image by giving input radec boundaries for the map. Boundaries are 
		defined as the center and "radius" in ra/dec. Much easier to input than corners of the boundary'''

		hpx_map = np.zeros(H.nside2npix(self.nside))

		self.superimpose_hpxmap(hpx_map,color=color)

	def superimpose_boundary_corners(self,radec_corners,color='red'):
		
		hpx_map = np.zeros(H.nside2npix(self.nside))

		self.superimpose_hpxmap(hpx_map,color=color)

def get_color_map(color):
	'''Generate a LinearSegmentedColormap with a single color and varying transparency'''

	from matplotlib.colors import LinearSegmentedColormap
	from matplotlib.colors import colorConverter

	if type(color) is str:
		rgb = colorConverter.to_rgb(color)
	elif len(color) == 3:
		#Not a string and length is 3, then assuming it is an rgb value
		rgb = color

	cdict = {'red':   [(0,rgb[0],rgb[0]),
					   (1,rgb[0],rgb[0])],
			'green':  [(0,rgb[1],rgb[1]),
					   (1,rgb[1],rgb[1])],
			'blue':   [(0,rgb[2],rgb[2]),
					   (1,rgb[2],rgb[2])],
			'alpha':  [(0,0,0),
					   (1,1,1)]}

	colormap1 = LinearSegmentedColormap('testcm',cdict)
	colormap1.set_bad(alpha=0.0)
	colormap1.set_under(alpha=0.0)

	return colormap1


def wcs_to_healpix(hdulist,nside):
	'''Converts data in an opened FITS file from a generic WCS to Healpix'''
	
	#use wcs_pix2world to get ra/dec coordinates, then use ang2pix to get Healpix coordinates and 
	#sum over all pixels with the same Healpix pixel numer

	w = wcs.WCS(hdulist[0].header)

	data = hdulist[0].data

	xsize, ysize = data.shape
	print "x/ysize = ", xsize, ysize

	pixcrd = [(y,x) for x in xrange(xsize) for y in xrange(ysize)]

	print "np.array"
	pixcrd2 = np.array(pixcrd)

	#pix2world to get ra/dec values of each pixel
	print "wcs_pix2world"
	world = w.wcs_pix2world(pixcrd,1)

	#needed since for some reason I was getting NaNs for some ACT data. Probably because of my swapping of x and y
	idx = np.isfinite(world[:,0])

	pixcrd2 = pixcrd2[idx,:]
	world = world[idx,:]

	#Use Healpy to do ang2pix for Healpix pixel numbers
	theta = np.pi/2.0 - np.radians(world[:,1])
	phi = np.radians(world[:,0]) #CHECK SIGN
	print "ang2pix"
	pixnum = H.ang2pix(nside,theta,phi)
	

	print "hpx_data"
	npix = H.nside2npix(nside)
	row = pixnum
	col = np.zeros_like(row)

	#coo_matrix will sum entries that have multiple elements
	hpx_data = scipy.sparse.coo_matrix((data[(pixcrd2[:,1],pixcrd2[:,0])], (row, col)), shape=np.array([npix,1])).A.transpose()
	hpx_data.shape = (npix)

	return hpx_data

if __name__=='__main__':

	fn_background = 'maps/Planck/HFI_SkyMap_353_2048_R2.00_full.fits'
	
	footprint = SurveyStack(fn_background,fignum=1)
	
	fn_148_south_s2 = 'maps/ACT/ACT_148_south_season_2_1way_hits_v3.fits'
	fn_148_equ_s3 = 'maps/ACT/ACT_148_equ_season_3_1way_hits_v3.fits'
	fn_148_south_s3 = 'maps/ACT/ACT_148_south_season_3_1way_hits_v3.fits'
	fn_148_equ_s4 = 'maps/ACT/ACT_148_equ_season_4_1way_hits_v3.fits'
	fn_148_south_s4 = 'maps/ACT/ACT_148_south_season_4_1way_hits_v3.fits'

	fn_220_equ_s4 = 'maps/ACT/ACT_220_equ_season_4_1way_hits_v3.fits'
	fn_220_south_s4 = 'maps/ACT/ACT_220_south_season_4_1way_hits_v3.fits'

	fns_ACT = [fn_148_equ_s4,fn_148_south_s4,fn_148_equ_s3,fn_148_south_s3,fn_148_south_s2]
	#fns_ACT = [fn_148_equ_s4]

	footprint.superimpose_fits(fns_ACT,color='red',maptype='WCS')
	
	pl.show()

