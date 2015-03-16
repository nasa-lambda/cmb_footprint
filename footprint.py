import numpy as np
import pylab as pl
import healpy as H
import matplotlib.cm as cm
from astropy.io import fits

import util

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
		'''Read in one or more healpix maps and add them together. Must input an array of strings even if only 
		inputting a single map'''
		
		hpx_map = np.zeros(H.nside2npix(self.nside))
		for fn in fns:
			tmp_map = H.read_map(fn)
			nside = H.npix2nside(len(tmp_map))
			hpx_map += H.ud_grade(tmp_map,self.nside)
	
		return hpx_map

	def read_wcs_maps(self,fns):
		
		hpx_map = np.zeros(H.nside2npix(self.nside))
		for fn in fns:
			hdulist = fits.open(fn)
			hpx_map += util.wcs_to_healpix(hdulist,self.nside)
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

		corner1 = (radec_cen[0]+radec_size[0],radec_cen[1]+radec_size[1])
		corner2 = (radec_cen[0]+radec_size[0],radec_cen[1]-radec_size[1])
		corner3 = (radec_cen[0]-radec_size[0],radec_cen[1]-radec_size[1])
		corner4 = (radec_cen[0]-radec_size[0],radec_cen[1]+radec_size[1])

		corners = (corner1,corner2,corner3,corner4)

		self.superimpose_boundary_corners(corners,color=color)

	def superimpose_boundary_corners(self,radec_corners,color='red'):
		
		radec_corners = np.array(radec_corners)

		thetas = np.pi/2 - np.radians(radec_corners[:,1])
		phis = np.radians(radec_corners[:,0])

		vecs = H.ang2vec(thetas,phis)

		ipix = H.query_polygon(self.nside,vecs)

		hpx_map = np.zeros(H.nside2npix(self.nside))
		hpx_map[ipix] = 1.0

		self.superimpose_hpxmap(hpx_map,color=color)

	def superimpose_experiment(self,experiment_name,color='red'):
		'''Superimpose a specific experiment whose Healpix footprints we have pregenerated'''

		if experiment_name == 'ACT':
			fns = ['maps/ACT_148_equ_hits_hpx.fits','maps/ACT_148_south_hits_hpx.fits']
		else:
			print('We do not have Healpix maps for this experiment')

		self.superimpose_fits(fns,color=color,maptype='HPX')

def get_color_map(color):
	'''Generate a LinearSegmentedColormap with a single color and varying transparency. Bad values and values below the 
	lower limit are set to be completely transparent.'''

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
			'alpha':  [(0,0.5,0.5),
					   (1,1,1)]}

	colormap1 = LinearSegmentedColormap('FootprintCM',cdict)
	colormap1.set_bad(alpha=0.0)
	colormap1.set_under(alpha=0.0)

	return colormap1

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

	#footprint.superimpose_fits(fns_ACT,color='red',maptype='WCS')
	#footprint.superimpose_boundary_cen((0,0),(40,5),color='blue')
	footprint.superimpose_experiment('ACT',color='red')

	pl.show()

