#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utility functions for the delsimi simulation code.

@author: Jonas Svenstrup Hansen, jonas.svenstrup@gmail.com
"""

import numpy as np
from astroquery.simbad import Simbad
from astropy import coordinates


def rvb2rgb(rvb):
	return rvb


def star_CCD_speed(pixel_scale, H_ISS=400*1e3, V_ISS=7.66*1e3, 
				R_Earth=6371*1e3):
	"""
	Estimate the speed of a star on the Delphini-1 CCD.
	
	The speed of a star on the CCD of Delphini-1 is estimated from knowledge of
	the pixel size, and the orbit height and speed of the ISS, which Delphini-1
	is assumed to follow. The orbit is assumed circular.
	
	Parameters
	----------
	pixel_scale (float):
		Size of a pixel on the CCD in arcsec per pixel.
	H_ISS (float):
		Approximate height of the ISS orbit in meters. Default is ``400*1e3``.
	V_ISS (float):
		Approximate speed of the ISS in orbit in meters per second. Default is 
		``7.66*1e3``.
	R_Earth (float):
		Approximate radius of Earth in meters. Default is ``6371*1e3``.
	
	Returns
	-------
	(float):
		Approximate speed of a star on the CCD of Delphini-1 in pixels per 
		second.
	"""

	return 360*3600/(2*np.pi*(R_Earth+H_ISS)) * V_ISS / pixel_scale


def bvr2rgb_discarded(bvr, A_inv=None, C=None):
	"""
	WARNING: This function does not yield believable output. It is kept here for
	reference, but is not to be used.
	
	Convert Johnson-Cousins B, V and R magnitudes to RGB magnitudes using the
	inversed conversion algortihm from Park (2016).

	Equations (1)-(3) in Park (2016) provide the conversion from RGB to BVR.
	Here, these equations are used in their matrix representation
		$$Ab+c=x,$$
	where `A` is the matrix of constants, `b` is the RGB vector, `c` is a vector 
	of constants and `x` is the BVR vector. The matrix `A` is then inverted and 
	the solution vector for the conversion from BVR to RGB is given as
		$$b = A^{-1} (x-c)$$

	The constants used here are taken from the Park (2016), Table 2, 2nd sigma-
	clipping. We use these values in the knowledge that they are found using 
	both a different camera system and with ground based photometry besides 
	being based on observations of a single open cluster, IC4665. However,
	they serve as a rough estimate, and if desired realistic constant could be
	optained from measurements.
	
	The article Park (2016) with the title "Photometric transformation from RGB 
	Bayer filter system to Johnson–Cousins BVR filter system" can be found here:
	https://www.sciencedirect.com/science/article/pii/S0273117715005694?via%3Dihub

	Parameters
	----------
	bvr (numpy array):
		Array with Johnson-Cousins B, V and R magnitudes: [B,V,R].
	A_inv (numpy array):
		Inverted constants matrix. Default is ``None``, which will create 
		``A_inv`` from scratch.
	C (numpy array):
		Constants vector. Default is ``None``, which will create ``C`` from 
		scratch.

	Returns
	-------
	rgb (numpy array):
		Array with converted RGB magnitudes: [R,G,B]

	Code Authors
	------------
	Andreas Kjær Dideriksen
	Jonas Svenstrup Hansen, jonas.svenstrup@gmail.com
	"""

	if A_inv is None:
		# Set constant values from the Park (2016):
		C_BBG = 0.280
		C_BGR = 0.600
		C_VBG = 0.542
		C_VGR = -0.064
		C_RBG = 0.051
		C_RGR = 0.468

		# Define matrix:
		A = np.array([[1+C_BBG, -C_BBG+C_BGR, -C_BGR],
					[C_VBG, 1-C_VBG+C_VGR, -C_VGR],
					[C_RBG, -C_RBG+C_RGR, 1-C_RGR]])

		# Invert matrix:
		A_inv = np.linalg.inv(A)

	if C is None:
		# Set constant values from the Park (2016):
		B_BZP = -0.291
		G_BZP = -0.252
		R_BZP = -0.226

		# Define vector:
		C = np.array([B_BZP, G_BZP, R_BZP])

	# Calculate and return rgb values:
	return np.flip((A_inv.dot(bvr - C)))


def make_bayer_filter(img_shape):
	"""
	Apply Bayer filter scaling of a normalised image with the level of flux.
	
	The Bayer pattern used here is BGGR, which corresponds to the following
	pixel filter setup for a 4x6 pixel sensor::

		B G B G B G
		G R G R G R
		B G B G B G
		G R G R G R

	Parameters
	----------
	img_shape (tuple):
		Shape of the image.

	Returns
	-------
	Bayer_flux (numpy array):
		2D image the same shape as the input, but with Bayer filter flags, 
		the structure of which is	visualized using the filter example above as 
		an example::
		
			2 1 2 1 2 1
			1 0 1 0 1 0
			2 1 2 1 2 1
			1 0 1 0 1 0
	"""

	# Preallocate:
	bayer_filter = np.empty(img_shape, dtype=int)

	# Red:
	bayer_filter[1::2,1::2] = 0
	# Green:
	bayer_filter[0::2,1::2] = 1
	bayer_filter[1::2,0::2] = 1
	# Blue:
	bayer_filter[0::2,0::2] = 2

	return bayer_filter



def mag2flux(mag, mag_type=None):
	"""
	Convert from magnitude to flux using scaling relation from
	aperture photometry for the TESS satellite. This is an estimate.

	Parameters
	----------
	mag (float): 
		Magnitude in TESS band.
	mag_type (string):
		Magnitude type. Can be either of ``R``, ``G`` and ``B``. 
		Accurate conversion constants not implemented.

	Returns
	-------
	float: 
		Corresponding flux value.
	"""

	# Set magnitude to flux conversion constants for RGB:
	# TODO: use measured constants for magnitude to flux conversion in RGB
	mag_consts_RGB = {
		"R": 28.24,
		"G": 28.24,
		"B": 28.24
	}

	return 10**(-0.4*(mag - mag_consts_RGB.get(mag_type, 28.24)))


def make_astroquery(ra=None, dec=None, radius=None, maxVmag=6.):
	"""
	Make a Simbad query for stars with a given maximum V magnitude in a radius 
	around a given position.
	
	Parameters
	----------
	ra (float):
		Right ascension in degrees (ICRSd).
	dec (float):
		Declination in degrees (ICRSd).
	radius (int):
		Maximum radius in degrees where to search around the specified (ra,dec)
		position.
	maxVmag (float):
		Maximum V magnitude of objects to include.
	
	Returns
	-------
	(numpy array):
		Array containing the columns (ra, dec, R, V, B), where the last three
		values are the Johnson-Cousins magnitudes of the given stars.
	"""
	# Use Pleiades cluster (ra,dec) in ICRSd found using Aladin as default:
	if ra or dec is None:
		ra, dec = (56.75,24.11670)
	if radius is None:
		radius = 3
	
	# Make astroquery:
	customSimbad = Simbad()
	customSimbad.add_votable_fields('ra(d)','dec(d)','flux(R)','flux(V)','flux(B)')
	customSimbad.remove_votable_fields('coordinates')
	
	C = coordinates.SkyCoord(ra,dec,unit=('deg','deg'), frame='icrs')
	result = customSimbad.query_region(C, radius=np.str(radius)+' degrees')
	
	# Select only objects brighter than some V magnitude:
	result = result[result['FLUX_V']<maxVmag]
	
	# Convert output to numpy array:
	return np.transpose(np.array([result['RA_d'],result['DEC_d'],
					result['FLUX_R'],result['FLUX_V'],result['FLUX_B']]))

