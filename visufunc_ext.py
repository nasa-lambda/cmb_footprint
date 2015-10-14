# 
#  This file copies parts of Healpy in order to add a specific function 
#  to a class that is inherited. As such, we are reproducing the license header
# 
#  This is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
# 
#  This is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with Healpy; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# 
#  For more information about Healpy, see http://code.google.com/p/healpy
# 

'''Code to plot contours of images on projections of Healpix maps.'''

import numpy as np
import healpy.projaxes as PA
import healpy.projector as P
import healpy.pixelfunc as pixelfunc
import healpy as hp
import matplotlib.axes
import pylab as pl

def mollcontour(map=None,NV=None,fig=None,rot=None,coord=None,unit='',
             xsize=800,title='Mollweide view',nest=False,
             min=None,max=None,flip='astro',
             remove_dip=False,remove_mono=False,
             gal_cut=0,
             format='%g',format2='%g',
             cbar=True,cmap=None, notext=False,
             norm=None,hold=False,margins=None,sub=None,
             return_projected_map=False):
    """Plot an healpix map (given as an array) as contours in Mollweide projection.
    
    Parameters
    ----------
    map : float, array-like or None
      An array containing the map, supports masked maps, see the `ma` function.
      If None, will display a blank map, useful for overplotting.
    NV : int, array-like, or None
        if int, the number of contour levels
        if an array, the values of the contour levels
        if None, contour levels will be automatically determined
    fig : int or None, optional
      The figure number to use. Default: create a new figure
    rot : scalar or sequence, optional
      Describe the rotation to apply.
      In the form (lon, lat, psi) (unit: degrees) : the point at
      longitude *lon* and latitude *lat* will be at the center. An additional rotation
      of angle *psi* around this direction is applied.
    coord : sequence of character, optional
      Either one of 'G', 'E' or 'C' to describe the coordinate
      system of the map, or a sequence of 2 of these to rotate
      the map from the first to the second coordinate system.
    unit : str, optional
      A text describing the unit of the data. Default: ''
    xsize : int, optional
      The size of the image. Default: 800
    title : str, optional
      The title of the plot. Default: 'Mollweide view'
    nest : bool, optional
      If True, ordering scheme is NESTED. Default: False (RING)
    min : float, optional
      The minimum range value
    max : float, optional
      The maximum range value
    flip : {'astro', 'geo'}, optional
      Defines the convention of projection : 'astro' (default, east towards left, west towards right)
      or 'geo' (east towards roght, west towards left)
    remove_dip : bool, optional
      If :const:`True`, remove the dipole+monopole
    remove_mono : bool, optional
      If :const:`True`, remove the monopole
    gal_cut : float, scalar, optional
      Symmetric galactic cut for the dipole/monopole fit.
      Removes points in latitude range [-gal_cut, +gal_cut]
    format : str, optional
      The format of the scale label. Default: '%g'
    format2 : str, optional
      Format of the pixel value under mouse. Default: '%g'
    cbar : bool, optional
      Display the colorbar. Default: True
    notext : bool, optional
      If True, no text is printed around the map
    norm : {'hist', 'log', None}
      Color normalization, hist= histogram equalized color mapping,
      log= logarithmic color mapping, default: None (linear color mapping)
    hold : bool, optional
      If True, replace the current Axes by a MollweideAxes.
      use this if you want to have multiple maps on the same
      figure. Default: False
    sub : int, scalar or sequence, optional
      Use only a zone of the current figure (same syntax as subplot).
      Default: None
    margins : None or sequence, optional
      Either None, or a sequence (left,bottom,right,top)
      giving the margins on left,bottom,right and top
      of the axes. Values are relative to figure (0-1).
      Default: None
    return_projected_map : bool
      if True returns the projected map in a 2d numpy array

    See Also
    --------
    gnomcontour, cartcontour, orthcontour
    """
    
    # Create the figure
    import pylab
    if not (hold or sub):
        f=pylab.figure(fig,figsize=(8.5,5.4))
        extent = (0.02,0.05,0.96,0.9)
    elif hold:
        f=pylab.gcf()
        left,bottom,right,top = np.array(f.gca().get_position()).ravel()
        extent = (left,bottom,right-left,top-bottom)
        f.delaxes(f.gca())
    else: # using subplot syntax
        f=pylab.gcf()
        if hasattr(sub,'__len__'):
            nrows, ncols, idx = sub
        else:
            nrows, ncols, idx = sub/100, (sub%100)/10, (sub%10)
        if idx < 1 or idx > ncols*nrows:
            raise ValueError('Wrong values for sub: %d, %d, %d'%(nrows,
                                                                 ncols,
                                                                 idx))
        c,r = (idx-1)%ncols,(idx-1)//ncols
        if not margins:
            margins = (0.01,0.0,0.0,0.02)
        extent = (c*1./ncols+margins[0], 
                  1.-(r+1)*1./nrows+margins[1],
                  1./ncols-margins[2]-margins[0],
                  1./nrows-margins[3]-margins[1])
        extent = (extent[0]+margins[0],
                  extent[1]+margins[1],
                  extent[2]-margins[2]-margins[0],
                  extent[3]-margins[3]-margins[1])
        #extent = (c*1./ncols, 1.-(r+1)*1./nrows,1./ncols,1./nrows)
    #f=pylab.figure(fig,figsize=(8.5,5.4))

    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = np.zeros(12)+np.inf
            cbar=False
        map = pixelfunc.ma_to_array(map)
        ax = HpxMollweideAxes(f,extent,coord=coord,rot=rot,
                               format=format2,flipconv=flip)
        f.add_axes(ax)
        if remove_dip:
            map=pixelfunc.remove_dipole(map,gal_cut=gal_cut,
                                        nest=nest,copy=True,
                                        verbose=True)
        elif remove_mono:
            map=pixelfunc.remove_monopole(map,gal_cut=gal_cut,nest=nest,
                                          copy=True,verbose=True)
        if NV is not None:
            img, cs = ax.projcontour(map,NV,nest=nest,xsize=xsize,coord=coord,vmin=min,vmax=max,
                                 cmap=cmap,norm=norm)
        else:
            img, cs = ax.projcontour(map,nest=nest,xsize=xsize,coord=coord,vmin=min,vmax=max,
                                 cmap=cmap,norm=norm)

        if cbar:
            b = cs.norm.inverse(np.linspace(0,1,cs.cmap.N+1))
            v = np.linspace(cs.norm.vmin,cs.norm.vmax,cs.cmap.N)
            if matplotlib.__version__ >= '0.91.0':
                cb=f.colorbar(cs,ax=ax,
                              orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
            else:
                # for older matplotlib versions, no ax kwarg
                cb=f.colorbar(cs,orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
            #cb.solids.set_rasterized(True)
        ax.set_title(title)
        if not notext:
            ax.text(0.86,0.05,ax.proj.coordsysstr,fontsize=14,
                    fontweight='bold',transform=ax.transAxes)
        if cbar:
            cb.ax.text(0.5,-1.0,unit,fontsize=14,
                       transform=cb.ax.transAxes,ha='center',va='center')
        f.sca(ax)
        #graticule(dpar=90, dmer=180.0, verbose=False)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
            #pylab.show()
    if return_projected_map:
        return img

def mollcontourf(map=None,NV=None,fig=None,rot=None,coord=None,unit='',
             xsize=800,title='Mollweide view',nest=False,
             min=None,max=None,flip='astro',
             remove_dip=False,remove_mono=False,
             gal_cut=0,
             format='%g',format2='%g',
             cbar=True,cmap=None, notext=False,
             norm=None,hold=False,margins=None,sub=None,
             return_projected_map=False):
    
    # Create the figure
    import pylab
    if not (hold or sub):
        f=pylab.figure(fig,figsize=(8.5,5.4))
        extent = (0.02,0.05,0.96,0.9)
    elif hold:
        f=pylab.gcf()
        left,bottom,right,top = np.array(f.gca().get_position()).ravel()
        extent = (left,bottom,right-left,top-bottom)
        f.delaxes(f.gca())
    else: # using subplot syntax
        f=pylab.gcf()
        if hasattr(sub,'__len__'):
            nrows, ncols, idx = sub
        else:
            nrows, ncols, idx = sub/100, (sub%100)/10, (sub%10)
        if idx < 1 or idx > ncols*nrows:
            raise ValueError('Wrong values for sub: %d, %d, %d'%(nrows,
                                                                 ncols,
                                                                 idx))
        c,r = (idx-1)%ncols,(idx-1)//ncols
        if not margins:
            margins = (0.01,0.0,0.0,0.02)
        extent = (c*1./ncols+margins[0], 
                  1.-(r+1)*1./nrows+margins[1],
                  1./ncols-margins[2]-margins[0],
                  1./nrows-margins[3]-margins[1])
        extent = (extent[0]+margins[0],
                  extent[1]+margins[1],
                  extent[2]-margins[2]-margins[0],
                  extent[3]-margins[3]-margins[1])
        #extent = (c*1./ncols, 1.-(r+1)*1./nrows,1./ncols,1./nrows)
    #f=pylab.figure(fig,figsize=(8.5,5.4))

    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = np.zeros(12)+np.inf
            cbar=False
        map = pixelfunc.ma_to_array(map)
        ax=HpxMollweideAxes(f,extent,coord=coord,rot=rot,
                               format=format2,flipconv=flip)
        f.add_axes(ax)
        if remove_dip:
            map=pixelfunc.remove_dipole(map,gal_cut=gal_cut,
                                        nest=nest,copy=True,
                                        verbose=True)
        elif remove_mono:
            map=pixelfunc.remove_monopole(map,gal_cut=gal_cut,nest=nest,
                                          copy=True,verbose=True)
        if NV is not None:
            img, cs = ax.projcontourf(map,NV,nest=nest,xsize=xsize,coord=coord,vmin=min,vmax=max,
                                 cmap=cmap,norm=norm)
        else:
            img, cs = ax.projcontourf(map,nest=nest,xsize=xsize,coord=coord,vmin=min,vmax=max,
                                 cmap=cmap,norm=norm)
        if cbar:
            b = cs.norm.inverse(np.linspace(0,1,cs.cmap.N+1))
            v = np.linspace(cs.norm.vmin,cs.norm.vmax,cs.cmap.N)
            if matplotlib.__version__ >= '0.91.0':
                cb=f.colorbar(cs,ax=ax,
                              orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
            else:
                # for older matplotlib versions, no ax kwarg
                cb=f.colorbar(cs,orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
            #cb.solids.set_rasterized(True)
        ax.set_title(title)
        if not notext:
            ax.text(0.86,0.05,ax.proj.coordsysstr,fontsize=14,
                    fontweight='bold',transform=ax.transAxes)
        if cbar:
            cb.ax.text(0.5,-1.0,unit,fontsize=14,
                       transform=cb.ax.transAxes,ha='center',va='center')
        f.sca(ax)
        #graticule(dpar=90, dmer=180.0, verbose=False)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
            #pylab.show()
    if return_projected_map:
        return img

def gnomcontour(map=None,NV=None,fig=None,rot=None,coord=None,unit='',
                xsize=200,ysize=None,reso=1.5,degree=False,
                title='Gnomonic view',nest=False,remove_dip=False,
                remove_mono=False,gal_cut=0,
                min=None,max=None,flip='astro',
                format='%.3g',cbar=True,
                cmap=None, norm=None,
                hold=False,sub=None,margins=None,notext=False,
                return_projected_map=False):
    """Plot an healpix map (given as an array) as contours in Gnomonic projection.

    Parameters
    ----------
    map : array-like
      The map to project, supports masked maps, see the `ma` function.
      If None, use a blank map, useful for
      overplotting.
    NV : int, array-like, or None
        if int, the number of contour levels
        if an array, the values of the contour levels
        if None, contour levels will be automatically determined
    fig : None or int, optional
      A figure number. Default: None= create a new figure
    rot : scalar or sequence, optional
      Describe the rotation to apply.
      In the form (lon, lat, psi) (unit: degrees) : the point at
      longitude *lon* and latitude *lat* will be at the center. An additional rotation
      of angle *psi* around this direction is applied.
    coord : sequence of character, optional
      Either one of 'G', 'E' or 'C' to describe the coordinate
      system of the map, or a sequence of 2 of these to rotate
      the map from the first to the second coordinate system.
    unit : str, optional
      A text describing the unit of the data. Default: ''
    xsize : int, optional
      The size of the image. Default: 200
    ysize : None or int, optional
      The size of the image. Default: None= xsize
    reso : float, optional
      Resolution (in arcmin if degree is False). Default: 1.5 arcmin
    degree : bool, optional
      if True, reso is in degree. Default: False
    title : str, optional
      The title of the plot. Default: 'Gnomonic view'
    nest : bool, optional
      If True, ordering scheme is NESTED. Default: False (RING)
    min : float, scalar, optional
      The minimum range value
    max : float, scalar, optional
      The maximum range value
    flip : {'astro', 'geo'}, optional
      Defines the convention of projection : 'astro' (default, east towards left, west towards right)
      or 'geo' (east towards roght, west towards left)
    remove_dip : bool, optional
      If :const:`True`, remove the dipole+monopole
    remove_mono : bool, optional
      If :const:`True`, remove the monopole
    gal_cut : float, scalar, optional
      Symmetric galactic cut for the dipole/monopole fit.
      Removes points in latitude range [-gal_cut, +gal_cut]
    format : str, optional
      The format of the scale label. Default: '%g'
    hold : bool, optional
      If True, replace the current Axes by a MollweideAxes.
      use this if you want to have multiple maps on the same
      figure. Default: False
    sub : int or sequence, optional
      Use only a zone of the current figure (same syntax as subplot).
      Default: None
    margins : None or sequence, optional
      Either None, or a sequence (left,bottom,right,top)
      giving the margins on left,bottom,right and top
      of the axes. Values are relative to figure (0-1).
      Default: None
    notext: bool, optional
      If True: do not add resolution info text. Default=False
    return_projected_map : bool
      if True returns the projected map in a 2d numpy array

    See Also
    --------
    mollcontour, cartcontour, orthcontour
    """
    
    import pylab
    if not (hold or sub):
        f=pylab.figure(fig,figsize=(5.8,6.4))
        if not margins:
                margins = (0.075,0.05,0.075,0.05)
        extent = (0.0,0.0,1.0,1.0)
    elif hold:
        f=pylab.gcf()
        left,bottom,right,top = np.array(pylab.gca().get_position()).ravel()
        if not margins:
            margins = (0.0,0.0,0.0,0.0)
        extent = (left,bottom,right-left,top-bottom)
        f.delaxes(pylab.gca())
    else: # using subplot syntax
        f=pylab.gcf()
        if hasattr(sub,'__len__'):
            nrows, ncols, idx = sub
        else:
            nrows, ncols, idx = sub/100, (sub%100)/10, (sub%10)
        if idx < 1 or idx > ncols*nrows:
            raise ValueError('Wrong values for sub: %d, %d, %d'%(nrows,
                                                                 ncols,
                                                                 idx))
        c,r = (idx-1)%ncols,(idx-1)//ncols
        if not margins:
            margins = (0.01,0.0,0.0,0.02)
        extent = (c*1./ncols+margins[0], 
                  1.-(r+1)*1./nrows+margins[1],
                  1./ncols-margins[2]-margins[0],
                  1./nrows-margins[3]-margins[1])
    extent = (extent[0]+margins[0],
              extent[1]+margins[1],
              extent[2]-margins[2]-margins[0],
              extent[3]-margins[3]-margins[1])
    #f=pylab.figure(fig,figsize=(5.5,6))

    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = np.zeros(12)+np.inf
            cbar=False
        map = pixelfunc.ma_to_array(map)
        ax=HpxGnomonicAxes(f,extent,coord=coord,rot=rot,
                              format=format,flipconv=flip)
        f.add_axes(ax)
        if remove_dip:
            map=pixelfunc.remove_dipole(map,gal_cut=gal_cut,nest=nest,copy=True)
        elif remove_mono:
            map=pixelfunc.remove_monopole(map,gal_cut=gal_cut,nest=nest,copy=True)

        if NV is not None:
            img, cs = ax.projcontour(map,NV,nest=nest,coord=coord,vmin=min,
                                     vmax=max,xsize=xsize,ysize=ysize,reso=reso,
                                     cmap=cmap,norm=norm)
        else:
            img, cs = ax.projcontour(map,nest=nest,coord=coord,vmin=min,
                                     vmax=max,xsize=xsize,ysize=ysize,reso=reso,
                                     cmap=cmap,norm=norm)
        
        if cbar:
            b = cs.norm.inverse(np.linspace(0,1,cs.cmap.N+1))
            v = np.linspace(cs.norm.vmin,cs.norm.vmax,cs.cmap.N)
            if matplotlib.__version__ >= '0.91.0':
                cb=f.colorbar(cs,ax=ax,
                              orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.08,fraction=0.1,boundaries=b,values=v,
                              format=format)
            else:
                cb=f.colorbar(cs,orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.08,fraction=0.1,boundaries=b,values=v,
                              format=format)
            cb.solids.set_rasterized(True)                  
        ax.set_title(title)
        if not notext:
            ax.text(-0.07,0.02,
                     "%g '/pix,   %dx%d pix"%(ax.proj.arrayinfo['reso'],
                                              ax.proj.arrayinfo['xsize'],
                                              ax.proj.arrayinfo['ysize']),
                     fontsize=12,verticalalignment='bottom',
                     transform=ax.transAxes,rotation=90)
            ax.text(-0.07,0.6,ax.proj.coordsysstr,fontsize=14,
                     fontweight='bold',rotation=90,transform=ax.transAxes)
            lon,lat = np.around(ax.proj.get_center(lonlat=True),ax._coordprec)
            ax.text(0.5,-0.03,'(%g,%g)'%(lon,lat),
                    verticalalignment='center', horizontalalignment='center',
                    transform=ax.transAxes)
        if cbar:
            cb.ax.text(1.05,0.30,unit,fontsize=14,fontweight='bold',
                       transform=cb.ax.transAxes,ha='left',va='center')
        f.sca(ax)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
            #pylab.show()
    if return_projected_map:
        return img

def cartcontour(map=None,NV=None,fig=None,rot=None,zat=None,coord=None,unit='',
                xsize=800,ysize=None,lonra=None,latra=None,
                title='Cartesian view',nest=False,remove_dip=False,
                remove_mono=False,gal_cut=0,
                min=None,max=None,flip='astro',
                format='%.3g',cbar=True,
                cmap=None, norm=None,aspect=None,
                hold=False,sub=None,margins=None,notext=False,
                return_projected_map=False):
    """Plot an healpix map (given as an array) as contours in Cartesian projection.

    Parameters
    ----------
    map : float, array-like or None
      An array containing the map, 
      supports masked maps, see the `ma` function.
      If None, will display a blank map, useful for overplotting.
    NV : int, array-like, or None
        if int, the number of contour levels
        if an array, the values of the contour levels
        if None, contour levels will be automatically determined
    fig : int or None, optional
      The figure number to use. Default: create a new figure
    rot : scalar or sequence, optional
      Describe the rotation to apply.
      In the form (lon, lat, psi) (unit: degrees) : the point at
      longitude *lon* and latitude *lat* will be at the center. An additional rotation
      of angle *psi* around this direction is applied.
    coord : sequence of character, optional
      Either one of 'G', 'E' or 'C' to describe the coordinate
      system of the map, or a sequence of 2 of these to rotate
      the map from the first to the second coordinate system.
    unit : str, optional
      A text describing the unit of the data. Default: ''
    xsize : int, optional
      The size of the image. Default: 800
    lonra : sequence, optional
      Range in longitude. Default: [-180,180]
    latra : sequence, optional
      Range in latitude. Default: [-90,90]
    title : str, optional
      The title of the plot. Default: 'Mollweide view'
    nest : bool, optional
      If True, ordering scheme is NESTED. Default: False (RING)
    min : float, optional
      The minimum range value
    max : float, optional
      The maximum range value
    flip : {'astro', 'geo'}, optional
      Defines the convention of projection : 'astro' (default, east towards left, west towards right)
      or 'geo' (east towards roght, west towards left)
    remove_dip : bool, optional
      If :const:`True`, remove the dipole+monopole
    remove_mono : bool, optional
      If :const:`True`, remove the monopole
    gal_cut : float, scalar, optional
      Symmetric galactic cut for the dipole/monopole fit.
      Removes points in latitude range [-gal_cut, +gal_cut]
    format : str, optional
      The format of the scale label. Default: '%g'
    cbar : bool, optional
      Display the colorbar. Default: True
    notext : bool, optional
      If True, no text is printed around the map
    norm : {'hist', 'log', None}, optional
      Color normalization, hist= histogram equalized color mapping,
      log= logarithmic color mapping, default: None (linear color mapping)
    hold : bool, optional
      If True, replace the current Axes by a CartesianAxes.
      use this if you want to have multiple maps on the same
      figure. Default: False
    sub : int, scalar or sequence, optional
      Use only a zone of the current figure (same syntax as subplot).
      Default: None
    margins : None or sequence, optional
      Either None, or a sequence (left,bottom,right,top)
      giving the margins on left,bottom,right and top
      of the axes. Values are relative to figure (0-1).
      Default: None
    return_projected_map : bool
      if True returns the projected map in a 2d numpy array

    See Also
    --------
    mollcontour, gnomcontour, orthcontour
    """
     
    import pylab
    if not (hold or sub):
        f=pylab.figure(fig,figsize=(8.5,5.4))
        if not margins:
                margins = (0.075,0.05,0.075,0.05)
        extent = (0.0,0.0,1.0,1.0)
    elif hold:
        f=pylab.gcf()
        left,bottom,right,top = np.array(pylab.gca().get_position()).ravel()
        if not margins:
            margins = (0.0,0.0,0.0,0.0)
        extent = (left,bottom,right-left,top-bottom)
        f.delaxes(pylab.gca())
    else: # using subplot syntax
        f=pylab.gcf()
        if hasattr(sub,'__len__'):
            nrows, ncols, idx = sub
        else:
            nrows, ncols, idx = sub/100, (sub%100)/10, (sub%10)
        if idx < 1 or idx > ncols*nrows:
            raise ValueError('Wrong values for sub: %d, %d, %d'%(nrows,
                                                                 ncols,
                                                                 idx))
        c,r = (idx-1)%ncols,(idx-1)//ncols
        if not margins:
            margins = (0.01,0.0,0.0,0.02)
        extent = (c*1./ncols+margins[0], 
                  1.-(r+1)*1./nrows+margins[1],
                  1./ncols-margins[2]-margins[0],
                  1./nrows-margins[3]-margins[1])
    extent = (extent[0]+margins[0],
              extent[1]+margins[1],
              extent[2]-margins[2]-margins[0],
              extent[3]-margins[3]-margins[1])

    #f=pylab.figure(fig,figsize=(5.5,6))
    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = np.zeros(12)+np.inf
            cbar=False
        map = pixelfunc.ma_to_array(map)
        if zat and rot:
            raise ValueError('Only give rot or zat, not both')
        if zat:
            rot = np.array(zat,dtype=np.float64)
            rot.resize(3)
            rot[1] -= 90
        ax=HpxCartesianAxes(f,extent,coord=coord,rot=rot,
                            format=format,flipconv=flip)
        f.add_axes(ax)
        if remove_dip:
            map=pixelfunc.remove_dipole(map,gal_cut=gal_cut,nest=nest,copy=True)
        elif remove_mono:
            map=pixelfunc.remove_monopole(map,gal_cut=gal_cut,nest=nest,copy=True)
        
        if NV is not None:
            img, cs = ax.projcontour(map,NV,nest=nest,coord=coord,vmin=min,
                                     vmax=max,xsize=xsize,ysize=ysize,lonra=lonra,
                                     latra=latra,cmap=cmap,norm=norm,aspect=aspect)
        else:
            img, cs = ax.projcontour(map,nest=nest,coord=coord,vmin=min,
                                     vmax=max,xsize=xsize,ysize=ysize,lonra=lonra,
                                     latra=latra,cmap=cmap,norm=norm,aspect=aspect)

        if cbar:
            b = cs.norm.inverse(np.linspace(0,1,cs.cmap.N+1))
            v = np.linspace(cs.norm.vmin,cs.norm.vmax,cs.cmap.N)
            if matplotlib.__version__ >= '0.91.0':
                cb=f.colorbar(cs,ax=ax,
                              orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.08,fraction=0.1,boundaries=b,values=v,
                              format=format)
            else:
                cb=f.colorbar(cs,orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.08,fraction=0.1,boundaries=b,values=v,
                              format=format)
            cb.solids.set_rasterized(True)                 
        ax.set_title(title)
        if not notext:
            ax.text(-0.07,0.6,ax.proj.coordsysstr,fontsize=14,
                     fontweight='bold',rotation=90,transform=ax.transAxes)
        if cbar:
            cb.ax.text(1.05,0.30,unit,fontsize=14,fontweight='bold',
                       transform=cb.ax.transAxes,ha='left',va='center')
        f.sca(ax)
    finally:
        if wasinteractive:
            pylab.ion()
            pylab.draw()
            #pylab.show()
    if return_projected_map:
        return img

def orthcontour(map=None,NV=None,fig=None,rot=None,coord=None,unit='',
             xsize=800,half_sky=False,
             title='Orthographic view',nest=False,
             min=None,max=None,flip='astro',
             remove_dip=False,remove_mono=False,
             gal_cut=0,
             format='%g',format2='%g',
             cbar=True,cmap=None, notext=False,
             norm=None,hold=False,margins=None,sub=None,
             return_projected_map=False):
    """Plot an healpix map (given as an array) in Orthographic projection.
    
    Parameters
    ----------
    map : float, array-like or None
      An array containing the map.
      If None, will display a blank map, useful for overplotting.
    NV : int, array-like, or None
        if int, the number of contour levels
        if an array, the values of the contour levels
        if None, contour levels will be automatically determined
    fig : int or None, optional
      The figure number to use. Default: create a new figure
    rot : scalar or sequence, optional
      Describe the rotation to apply.
      In the form (lon, lat, psi) (unit: degrees) : the point at
      longitude *lon* and latitude *lat* will be at the center. An additional rotation
      of angle *psi* around this direction is applied.
    coord : sequence of character, optional
      Either one of 'G', 'E' or 'C' to describe the coordinate
      system of the map, or a sequence of 2 of these to rotate
      the map from the first to the second coordinate system.
    half_sky : bool, optional
      Plot only one side of the sphere. Default: False
    unit : str, optional
      A text describing the unit of the data. Default: ''
    xsize : int, optional
      The size of the image. Default: 800
    title : str, optional
      The title of the plot. Default: 'Orthographic view'
    nest : bool, optional
      If True, ordering scheme is NESTED. Default: False (RING)
    min : float, optional
      The minimum range value
    max : float, optional
      The maximum range value
    flip : {'astro', 'geo'}, optional
      Defines the convention of projection : 'astro' (default, east towards left, west towards right)
      or 'geo' (east towards roght, west towards left)
    remove_dip : bool, optional
      If :const:`True`, remove the dipole+monopole
    remove_mono : bool, optional
      If :const:`True`, remove the monopole
    gal_cut : float, scalar, optional
      Symmetric galactic cut for the dipole/monopole fit.
      Removes points in latitude range [-gal_cut, +gal_cut]
    format : str, optional
      The format of the scale label. Default: '%g'
    format2 : str, optional
      Format of the pixel value under mouse. Default: '%g'
    cbar : bool, optional
      Display the colorbar. Default: True
    notext : bool, optional
      If True, no text is printed around the map
    norm : {'hist', 'log', None}
      Color normalization, hist= histogram equalized color mapping,
      log= logarithmic color mapping, default: None (linear color mapping)
    hold : bool, optional
      If True, replace the current Axes by an OrthographicAxes.
      use this if you want to have multiple maps on the same
      figure. Default: False
    sub : int, scalar or sequence, optional
      Use only a zone of the current figure (same syntax as subplot).
      Default: None
    margins : None or sequence, optional
      Either None, or a sequence (left,bottom,right,top)
      giving the margins on left,bottom,right and top
      of the axes. Values are relative to figure (0-1).
      Default: None
    return_projected_map : bool
      if True returns the projected map in a 2d numpy array
    
    See Also
    --------
    mollcontour, gnomcontour, cartcontour
    """
    # Create the figure
    import pylab
    if not (hold or sub):
        f=pylab.figure(fig,figsize=(8.5,5.4))
        extent = (0.02,0.05,0.96,0.9)
    elif hold:
        f=pylab.gcf()
        left,bottom,right,top = np.array(f.gca().get_position()).ravel()
        extent = (left,bottom,right-left,top-bottom)
        f.delaxes(f.gca())
    else: # using subplot syntax
        f=pylab.gcf()
        if hasattr(sub,'__len__'):
            nrows, ncols, idx = sub
        else:
            nrows, ncols, idx = sub/100, (sub%100)/10, (sub%10)
        if idx < 1 or idx > ncols*nrows:
            raise ValueError('Wrong values for sub: %d, %d, %d'%(nrows,
                                                                 ncols,
                                                                 idx))
        c,r = (idx-1)%ncols,(idx-1)//ncols
        if not margins:
            margins = (0.01,0.0,0.0,0.02)
        extent = (c*1./ncols+margins[0],
                  1.-(r+1)*1./nrows+margins[1],
                  1./ncols-margins[2]-margins[0],
                  1./nrows-margins[3]-margins[1])
        extent = (extent[0]+margins[0],
                  extent[1]+margins[1],
                  extent[2]-margins[2]-margins[0],
                  extent[3]-margins[3]-margins[1])
        #extent = (c*1./ncols, 1.-(r+1)*1./nrows,1./ncols,1./nrows)
    #f=pylab.figure(fig,figsize=(8.5,5.4))
    
    # Starting to draw : turn interactive off
    wasinteractive = pylab.isinteractive()
    pylab.ioff()
    try:
        if map is None:
            map = np.zeros(12)+np.inf
            cbar=False
        ax=HpxOrthographicAxes(f,extent,coord=coord,rot=rot,
                               format=format2,flipconv=flip)
        f.add_axes(ax)
        if remove_dip:
            map=pixelfunc.remove_dipole(map,gal_cut=gal_cut,
                                        nest=nest,copy=True,
                                        verbose=True)
        elif remove_mono:
            map=pixelfunc.remove_monopole(map,gal_cut=gal_cut,nest=nest,
                                          copy=True,verbose=True)
        
        if NV is not None:
            img, cs = ax.projcontour(map,NV,nest=nest,half_sky=half_sky,
                                     coord=coord,vmin=min,vmax=max,cmap=cmap,
                                     norm=norm)
        else:
            img, cs = ax.projcontour(map,nest=nest,half_sky=half_sky,
                                     coord=coord,vmin=min,vmax=max,cmap=cmap,
                                     norm=norm)
        if cbar:
            b = cs.norm.inverse(np.linspace(0,1,cs.cmap.N+1))
            v = np.linspace(cs.norm.vmin,cs.norm.vmax,cs.cmap.N)
            if matplotlib.__version__ >= '0.91.0':
                cb=f.colorbar(cs,ax=ax,
                              orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
            else:
                # for older matplotlib versions, no ax kwarg
                cb=f.colorbar(cs,orientation='horizontal',
                              shrink=0.5,aspect=25,ticks=PA.BoundaryLocator(),
                              pad=0.05,fraction=0.1,boundaries=b,values=v,
                              format=format)
            cb.solids.set_rasterized(True)                  
        ax.set_title(title)
        if not notext:
            ax.text(0.86,0.05,ax.proj.coordsysstr,fontsize=14,
                    fontweight='bold',transform=ax.transAxes)
        if cbar:
            cb.ax.text(0.5,-1.0,unit,fontsize=14,
                       transform=cb.ax.transAxes,ha='center',va='center')
        f.sca(ax)
    finally:
        pylab.draw()
        if wasinteractive:
            pylab.ion()
            #pylab.show()
    if return_projected_map:
        return img

class SphericalProjAxes(PA.SphericalProjAxes):

    def projcontour(self, map, vec2pix_func, *args, **kwds):
        """Projcontour is a wrapper around :func:`matplotlib.Axes.contour` to 
        take into account the spherical projection.

        Notes
        -----
        Other keywords are transmitted to :func:`matplotlib.Axes.contour`
        """
        vmin = kwds.pop('vmin', None)
        vmax = kwds.pop('vmax', None)
        norm = kwds.pop('norm', None)
        rot = kwds.pop('rot', None)
        coord = kwds.pop('coord', None)
        badval = kwds.pop('badval', hp.UNSEEN)
        cmap = kwds.pop('cmap', None)
        img = self.proj.projmap(map,vec2pix_func,rot=rot,coord=coord)
        w = ~( np.isnan(img) | 
               np.isinf(img) | 
               pixelfunc.mask_bad(img, badval = badval) )
        try:
            if vmin is None: vmin = img[w].min()
        except ValueError:
            vmin = 0.
        try:
            if vmax is None: vmax = img[w].max()
        except ValueError:
            vmax = 0.
        if vmin > vmax:
            vmin = vmax
        if vmin == vmax:
            vmin -= 1.
            vmax += 1.
        cm,nn = PA.get_color_table(vmin,vmax,img[w],cmap=cmap,norm=norm)
        ext = self.proj.get_extent()
        img = np.ma.masked_values(img, badval)
        axcont = self.contour(img, *args, extent = ext, cmap=cm,
                              norm=nn, origin='lower', vmin=vmin, 
                              vmax=vmax, **kwds)
        xmin,xmax,ymin,ymax = self.proj.get_extent()
        self.set_xlim(xmin,xmax)
        self.set_ylim(ymin,ymax)

        return img, axcont
    
    def projcontourf(self, map, vec2pix_func, *args, **kwds):
        """Projcontour is a wrapper around :func:`matplotlib.Axes.contour` to 
        take into account the spherical projection.

        Notes
        -----
        Other keywords are transmitted to :func:`matplotlib.Axes.contour`
        """
        vmin = kwds.pop('vmin', None)
        vmax = kwds.pop('vmax', None)
        norm = kwds.pop('norm', None)
        rot = kwds.pop('rot', None)
        coord = kwds.pop('coord', None)
        badval = kwds.pop('badval', hp.UNSEEN)
        cmap = kwds.pop('cmap', None)
        img = self.proj.projmap(map,vec2pix_func,rot=rot,coord=coord)
        w = ~( np.isnan(img) | 
               np.isinf(img) | 
               pixelfunc.mask_bad(img, badval = badval) )
        try:
            if vmin is None: vmin = img[w].min()
        except ValueError:
            vmin = 0.
        try:
            if vmax is None: vmax = img[w].max()
        except ValueError:
            vmax = 0.
        if vmin > vmax:
            vmin = vmax
        if vmin == vmax:
            vmin -= 1.
            vmax += 1.
        cm,nn = PA.get_color_table(vmin,vmax,img[w],cmap=cmap,norm=norm)
        ext = self.proj.get_extent()
        img = np.ma.masked_values(img, badval)
        axcont = self.contourf(img, *args, extent = ext, cmap=cm,
                              norm=nn, origin='lower', vmin=vmin, 
                              vmax=vmax, **kwds)
        xmin,xmax,ymin,ymax = self.proj.get_extent()
        self.set_xlim(xmin,xmax)
        self.set_ylim(ymin,ymax)

        return img, axcont

class GnomonicAxes(SphericalProjAxes):
    """Define a gnomonic Axes to handle gnomonic projection.

    Input:
      - rot=, coord= : define rotation and coordinate system. See rotator.
      - coordprec= : number of digit after floating point for coordinates display.
      - format= : format string for value display.
      
      Other keywords from Axes (see Axes).
    """
    def __init__(self,*args,**kwds):
        kwds.setdefault('coordprec',3)
        super(GnomonicAxes,self).__init__(P.GnomonicProj, *args,**kwds)
        self._do_border = False
        self._gratdef['local'] = True
        self._gratdef['dpar'] = 1.

    def projmap(self,map,vec2pix_func,xsize=200,ysize=None,reso=1.5,**kwds):
        self.proj.set_proj_plane_info(xsize=xsize,ysize=ysize,reso=reso)
        return super(GnomonicAxes,self).projmap(map,vec2pix_func,**kwds)
    
    def projcontour(self,map,vec2pix_func,*args,**kwds):
        xsize = kwds.pop('xsize', 200)
        ysize = kwds.pop('ysize', None)
        reso = kwds.pop('reso', 1.5)
        self.proj.set_proj_plane_info(xsize=xsize,ysize=ysize,reso=reso)
        return super(GnomonicAxes,self).projcontour(map,vec2pix_func,*args,**kwds)
    
    def projcontourf(self,map,vec2pix_func,*args,**kwds):
        xsize = kwds.pop('xsize', 200)
        ysize = kwds.pop('ysize', None)
        reso = kwds.pop('reso', 1.5)
        self.proj.set_proj_plane_info(xsize=xsize,ysize=ysize,reso=reso)
        return super(GnomonicAxes,self).projcontourf(map,vec2pix_func,*args,**kwds)
        
class HpxGnomonicAxes(GnomonicAxes):
    def projmap(self,map,nest=False,**kwds):
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        xsize = kwds.pop('xsize',200)
        ysize = kwds.pop('ysize',None)
        reso = kwds.pop('reso',1.5)
        return super(HpxGnomonicAxes,self).projmap(map,f,xsize=xsize,
                                            ysize=ysize,reso=reso,**kwds)
    
    def projcontour(self,map,*args,**kwds):
        nest = kwds.pop('nest', False)
        xsize = kwds.pop('xsize',200)
        ysize = kwds.pop('ysize',None)
        reso = kwds.pop('reso',1.5)
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxGnomonicAxes,self).projcontour(map,f,*args,xsize=xsize,
                                            ysize=ysize,reso=reso,**kwds)
    
    def projcontourf(self,map,*args,**kwds):
        nest = kwds.pop('nest', False)
        xsize = kwds.pop('xsize',200)
        ysize = kwds.pop('ysize',None)
        reso = kwds.pop('reso',1.5)
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxGnomonicAxes,self).projcontourf(map,f,*args,xsize=xsize,
                                            ysize=ysize,reso=reso,**kwds)

class MollweideAxes(SphericalProjAxes):
    """Define a mollweide Axes to handle mollweide projection.

    Input:
      - rot=, coord= : define rotation and coordinate system. See rotator.
      - coordprec= : number of digit after floating point for coordinates display.
      - format= : format string for value display.
      
      Other keywords from Axes (see Axes).
    """
    def __init__(self,*args,**kwds):
        kwds.setdefault('coordprec',2)
        super(MollweideAxes,self).__init__(P.MollweideProj, *args,**kwds)
        self.set_xlim(-2.01,2.01)
        self.set_ylim(-1.01,1.01)

    def projmap(self,map,vec2pix_func,xsize=800,**kwds):
        self.proj.set_proj_plane_info(xsize=xsize)
        img = super(MollweideAxes,self).projmap(map,vec2pix_func,**kwds)
        self.set_xlim(-2.01,2.01)
        self.set_ylim(-1.01,1.01)
        return img

    def projcontour(self,map,vec2pix_func,*args,**kwds):
        xsize = kwds.pop('xsize', 800)
        self.proj.set_proj_plane_info(xsize=xsize)
        img = super(MollweideAxes,self).projcontour(map,vec2pix_func,*args,**kwds)
        self.set_xlim(-2.01,2.01)
        self.set_ylim(-1.01,1.01)
        return img
    
    def projcontourf(self,map,vec2pix_func,*args,**kwds):
        xsize = kwds.pop('xsize', 800)
        self.proj.set_proj_plane_info(xsize=xsize)
        img = super(MollweideAxes,self).projcontourf(map,vec2pix_func,*args,**kwds)
        self.set_xlim(-2.01,2.01)
        self.set_ylim(-1.01,1.01)
        return img
        
class HpxMollweideAxes(MollweideAxes):
    def projmap(self,map,nest=False,**kwds):
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxMollweideAxes,self).projmap(map,f,**kwds)

    def projcontour(self,map,*args,**kwds):
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        nest = kwds.pop('nest', False)
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxMollweideAxes,self).projcontour(map,f,*args,**kwds)
    
    def projcontourf(self,map,*args,**kwds):
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        nest = kwds.pop('nest', False)
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxMollweideAxes,self).projcontourf(map,f,*args,**kwds)

class CartesianAxes(SphericalProjAxes):
    """Define a cylindrical Axes to handle cylindrical projection.
    """
    def __init__(self,*args,**kwds):
        kwds.setdefault('coordprec',2)
        super(CartesianAxes,self).__init__(P.CartesianProj, *args, **kwds)
        self._segment_threshold = 180
        self._segment_step_rad = 0.1*np.pi/180
        self._do_border = True
        
    def projmap(self,map,vec2pix_func,xsize=800,ysize=None,lonra=None,latra=None,**kwds):
        self.proj.set_proj_plane_info(xsize=xsize,ysize=ysize,lonra=lonra,latra=latra)
        return super(CartesianAxes,self).projmap(map,vec2pix_func,**kwds)
    
    def projcontour(self,map,vec2pix_func,*args,**kwds):
        xsize = kwds.pop('xsize', 800)
        ysize = kwds.pop('ysize', None)
        lonra = kwds.pop('lonra', None)
        latra = kwds.pop('latra', None)
        self.proj.set_proj_plane_info(xsize=xsize, ysize=ysize, lonra=lonra,
                                      latra=latra)
        img = super(CartesianAxes,self).projcontour(map,vec2pix_func,*args,**kwds)
        return img
    
    def projcontourf(self,map,vec2pix_func,*args,**kwds):
        xsize = kwds.pop('xsize', 800)
        ysize = kwds.pop('ysize', None)
        lonra = kwds.pop('lonra', None)
        latra = kwds.pop('latra', None)
        self.proj.set_proj_plane_info(xsize=xsize, ysize=ysize, lonra=lonra,
                                      latra=latra, **kwds)
        img = super(CartesianAxes,self).projcontourf(map,vec2pix_func,*args,**kwds)
        return img
        
class HpxCartesianAxes(CartesianAxes):
    
    def projmap(self,map,nest=False,**kwds):
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxCartesianAxes,self).projmap(map,f,**kwds)
    
    def projcontour(self,map,*args,**kwds):
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        nest = kwds.pop('nest', False)
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxCartesianAxes,self).projcontour(map,f,*args,**kwds)
    
    def projcontourf(self,map,*args,**kwds):
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        nest = kwds.pop('nest', False)
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxCartesianAxes,self).projcontourf(map,f,*args,**kwds)

        
class OrthographicAxes(SphericalProjAxes):
    """Define an orthographic Axes to handle orthographic projection.
    
    Input:
    - rot=, coord= : define rotation and coordinate system. See rotator.
    - coordprec= : num of digits after floating point for coordinates display.
    - format= : format string for value display.
    
    Other keywords from Axes (see Axes).
    """
    def __init__(self,*args,**kwds):
        kwds.setdefault('coordprec',2)
        super(OrthographicAxes,self).__init__(P.OrthographicProj, *args,**kwds)
        self._segment_threshold = 0.01
        self._do_border = False
        
    def projmap(self,map,vec2pix_func,xsize=800,half_sky=False,**kwds):
        self.proj.set_proj_plane_info(xsize=xsize,half_sky=half_sky)
        img = super(OrthographicAxes,self).projmap(map,vec2pix_func,**kwds)
        if half_sky: ratio = 1.01
        else: ratio = 2.01
        self.set_xlim(-ratio,ratio)
        self.set_ylim(-1.01,1.01)
        return img
    
    def projcontour(self,map,vec2pix_func,*args,**kwds):
        xsize = kwds.pop('xsize', 800)
        half_sky = kwds.pop('half_sky', False)
        self.proj.set_proj_plane_info(xsize=xsize,half_sky=half_sky)
        img = super(OrthographicAxes,self).projcontour(map,vec2pix_func,*args,**kwds)
        if half_sky: ratio = 1.01
        else: ratio = 2.01
        self.set_xlim(-ratio,ratio)
        self.set_ylim(-1.01,1.01)
        return img
    
    def projcontourf(self,map,vec2pix_func,*args,**kwds):
        xsize = kwds.pop('xsize', 800)
        half_sky = kwds.pop('half_sky', False)
        self.proj.set_proj_plane_info(xsize=xsize,half_sky=half_sky)
        img = super(OrthographicAxes,self).projcontourf(map,vec2pix_func,*args,**kwds)
        if half_sky: ratio = 1.01
        else: ratio = 2.01
        self.set_xlim(-ratio,ratio)
        self.set_ylim(-1.01,1.01)
        return img
        
class HpxOrthographicAxes(OrthographicAxes):
    
    def projmap(self,map,nest=False,**kwds):
        nside = pixelfunc.npix2nside(len(map))
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxOrthographicAxes,self).projmap(map,f,**kwds)
    
    def projcontour(self,map,*args,**kwds):
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        nest = kwds.pop('nest', False)
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxOrthographicAxes,self).projcontour(map,f,*args,**kwds)
    
    def projcontourf(self,map,*args,**kwds):
        nside = pixelfunc.npix2nside(pixelfunc.get_map_size(map))
        nest = kwds.pop('nest', False)
        f = lambda x,y,z: pixelfunc.vec2pix(nside,x,y,z,nest=nest)
        return super(HpxOrthographicAxes,self).projcontourf(map,f,*args,**kwds)

