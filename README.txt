DELSIMI - A SIMULATION CODE FOR ASTRONOMY IMAGES WITH DELPHINI-1
****************************************************************
Codeauthor: Jonas Svenstrup Hansen, jonas.svenstrup@gmail.com, 2018


DESCRIPTION
===========
This code serves to simulate astronomy images from the Delphini-1 satellite.
For information about this satellite see http://projects.au.dk/ausat/delphini1/.
This file explains the structure of the code and how to run it. For information
on the development of the various method see the file protocol.pdf.

.. note:: The code does not yet support input catalog stars!
.. note:: Angle definitions are still in development!


STRUCTURE
---------
The code is split in a run file run_delsimi.py which can be executed from the 
commandline as described in `RUNNING`_. This file runs the Delsimi class located
in delsimi.py. This class defines basic parameters, but uses the PSF class 
located in psf.py to define and integrate the point spread function. Various
functions used by either class are available in utilities.py.

UML diagrams are shown in packages.pdf and classes.pdf. They are created using
pyreverse as guided by 
https://twigstechtips.blogspot.dk/2014/09/python-visualise-your-class-hierarchy.html


I/O
---
The input and output are loaded and written to directories ''infiles'' and 
''outfiles'', which should be located in the parent directory of this code.


RUNNING
-------
The code is designed to be run on Ubuntu 16.04 with Python 3.5 using the 
packages listed in the file requirements.txt. The following procedure describes 
how to run the simulation.

1. Create input and output folders described in `I/O`_.
2. In a virtual Python 3.5 environment install the requirements.
3. In a terminal in the simulation directory run the file run_delsimi.py.

HOW TO run_delsimi
++++++++++++++++++
The file run_delsimi.py can be run from a bash terminal using this command::
	
	python run_delsimi.py

This call will generate a FITS file image.fits in the output directory along
with an ASCII file catalog.txt containing a description of the catalog of stars
used to create the image. The stars in the image will be two default stars.
Several parameters can also be passed, e.g.::
	
	python run_delsimi.py -c [45,60] -t 5. -v 0.1 -s 0.2

This call will place the "viewfinder" at an angle of 0.2 radians at right
ascension 45 and declination 60. The two default stars will be integrated for
five seconds, and they will move at an angle of 0.1 in the image.

Show all available parameters by executing::

	python run_delsimi.py --help


FUTURE EXTENSIONS
=================
Several extensions are proposed in the following which is largely inspired by 
the SPyFFI manual (<https://github.com/TESScience/SPyFFI>`_).
- Pixel-integrated Gaussian at high time resolution to replace convolution?
- Catalog stars
- Position variable PSF (not compatible with convolution)
- Realistic constant used for magnitude to flux conversion (color dependent)
- Realistic Johnson-Cousins BVR to RGB magnitude conversion constants
- Background
- Focus change
- Jitter
- Pixel sensitivity variations
- Photon noise
- Readout smear
- Rolling shutter
- Saturated pixels
- Cosmic rays


DISTRIBUTION
============
The code is available on the online repository 
<https://github.com/jonasshansen/delsimi/>`_
