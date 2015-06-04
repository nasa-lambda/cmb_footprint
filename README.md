# README #

This is the README for the CMB footprint module provided by LAMBDA.

### What is this repository for? ###

The CMB footprint module is designed to show where different experiments are
observing on the sky. All experiments for which we provide survey footprints
are listed in the footprint.cfg file.

### How do I get set up? ###

The code runs on both Python 2 and Python 3. The Python modules that are
required are

* Numpy
* Matplotlib
* Healpy
* Astropy
* Scipy

As long as your PYTHONPATH environment variable points to the directory 
where this module is located, the code should be able to run.

### Overview of the library ###

There are a couple examples on how to run the code shown in the examples/
subdirectory. A background for the footprint is needed for the code to run.
Any Healpix map can be used. The map used in the example scripts is not
included, but can be downloaded from the Planck archive.

There are many different options when plotting the survey footprints. You
can plot using any projection that Healpy provides and use most of the options
that these projections provide. The coordinate system of the background map
and the coordinate system of the plot should be provided though they do have
default values ('G' and 'C' respectively). The coordinate system of each
experiment in the configuration file is assumed to be in Equatorial
coordinates, though if an input map is provided instead of using the
definitions in the configuration file, a coordinate system should be provided.

The input colors can be any string understood by matplotlib or an rgb triplet.
The plotted regions for single experiment will have the same color but a
varying (0.5 to 1) alpha depending on the pixel value.

A copy of the latest configuration file is additionally stored on LAMBDA and
an option can be set to download the latest version of the configuration file
every time the code is run.

### Who do I talk to? ###

* Nathan Miller
