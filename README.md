# README #

This is the README for the CMB footprint module provided by LAMBDA.

## What is this repository for? ##

The CMB footprint module is designed to show where different surveys are
observing on the sky. All surveys for which we provide footprint regions
are listed in the footprint.cfg file.

## How do I get set up? ##

The code runs on both Python 2 (checken on 2.7.10) and Python 3 (checked on 
3.4.3). The Python modules that are required for the code to run are

* Numpy
* Matplotlib
* Healpy
* Astropy
* Scipy

Numpy, Matplotlib, and Healpy are fairly obvious requirements. Astropy is 
needed so that human readable values can be input into the configuration
file. Additionally, it can be used to read in WCS FITS files. Scipy is only
used for the code that converts a WCS map to a Healpix map.

As long as your PYTHONPATH environment variable points to the directory 
where this module is located, the code should be able to run.

## Overview of the library ##

There are a couple examples on how to run the code shown in the examples/
subdirectory. Additionally, a Jupyter notebook is provided that shows the
results from both of the example scripts. A background for the footprint
is needed for the code to run. Backgrounds can be added to the configuration
file, but will need to be downloaded the first time and the files can be large.
Alternatively, any Healpix map can be input. 

There are many different options when plotting the survey footprints. You
can plot using any projection that Healpy provides and use most of the options
that these projections provide. The coordinate system of the background map
and the coordinate system of the plot should be provided though they do have
default values ('G' and 'C' respectively). The coordinate system of each
survey in the configuration file is assumed to be in Equatorial
coordinates, though if an input map is provided instead of using the
definitions in the configuration file, a coordinate system should be provided.

The input colors can be any string understood by matplotlib or an rgb triplet.
The plotted regions for single survey will have the same color but a
varying (0.5 to 1) alpha depending on the pixel value.

A copy of the latest configuration file is additionally stored on LAMBDA and
an option can be set to download the latest version of the configuration file
every time the code is run.

## Types of Configuration File Entries ##

There are several different ways to input surveys to the configuration
file. Here we will go over all the different types. For each types I will
list the different keys that must go in the configuration file for that entry. The 
regions are specified as strings in the standard astropy [coordinate format](http://astropy.readthedocs.org/en/latest/coordinates/).

### Common Keys ###

These are keys that are common to all the different types of configuration file entries.

* **coord** : 'C', 'E', or 'G' describing the coordinate system of the input map or 
the coordinate system of the input values.
* **instrument** : Instrument with which this data comes from. Multiple entries can
have the same instrument.
* **survey_type** : Either 'background', 'cmb', 'cmbpol', or 'lss'. Used to separate 
entries into different tables when we generate a markdown file listing all the 
entries in the configuration file.
* **citation** : Citation describing how we got our survey footprint for the given 
experiment. Possible entries include a link to a file that we got our data form or 
and link to a paper from which we pulled values.


### Healpix File ###

This reads in a footprint stored as a Healpix file. If the file does not 
exist in the path specified when you initialize the class, it will attempt
to download the file from LAMBDA

* **handler** : hpx_file
* **file** : filename of the file
* **checksum** : The MD5 checksum of the file. It the checksum does not match the
checksum of the local file, it will attempt to download it from LAMBDA.

### Disc ###

This generates a disc footprint given a location and radius for the disc.

* **handler** : disc
* **center** : The center of the disc in lon,lat in human readable form
(i.e. 4h12m,-12d4m). Values must be separated by a comma
* **radius** : The radius of the disc in human readable form

### Polygon ###

This generates a footprint given vertices of a polygon. There is no limit to
the number of vertices, but the resulting polygon must be convex.

* **handler** : polygon
* **vertex1**, **vertex2**, ... : Lon,lat location of the vertex. The lon,lat values
must be separated by a comma

### Rectangle ###

This type generates a rectangle for the footprint. The rectangle is generated
from a center point and a length of each side. The rectangle is always
oriented along ra/dec lines. This will only look like a rectangle in a
cartesian like projection where the x- and y-axes are ra/dec.

* **handler** : rectangle
* **center** : The center lon,lat of the rectangle.
* **size** : The length of the sides of the rectangle.

### Combination ###

The last option is a combination option which combines the footprints of
multiple other surveys listed in the configuration file. This is used so
that a single survey can have multiple entries for different patches, so
we can choose to plot a single patch or multiple patches without the code
thinking they are different surveys. All component maps must have the same
coordinate system because we just sum the Healpix maps together.

* **handler** : combination
* **components** : The names of the other entries in the configuration file that
will be combined. Names must be separated by a comma.

## Examples and Config File ##

Examples are shown in the examples directory. The Jupyter (.ipynb) notebook 
(https://github.com/nasa-lambda/cmb_footprint/blob/master/examples/plot_footprints.ipynb) 
can be visualized by github. The output is what I get on my computer.

The entries in the config file can be seen at 
https://github.com/nasa-lambda/cmb_footprint/blob/master/footprint_config.md or
http://lambda.gsfc.nasa.gov/toolbox/footprint/configfile.cfm. These pages 
are automatically generated by gen_footprint_md.py based on the entries in 
the configuration file and is an easier way of seeing what surveys are 
predefined.

## Who do I talk to? ##

* Nathan Miller (nathan.j.miller@nasa.gov)
