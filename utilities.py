#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities for the delsimi simulation code.

@author: jonas
"""

import numpy as np



def rvb2rgb(rvb):
	return rvb



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
	return np.flip((A_inv.dot(bvr - C))


def bayer_scaling(img, flux):
	"""
	Apply Bayer filter scaling of a normalised image with the level of flux.
	
	The Bayer pattern used here is BGGR.
	
	Parameters
	----------
	img (numpy array):
		2D image with an even number of rows and columns.
	flux (list of floats):
		List with the three values for red, green and blue to scale image.
	
	Returns
	-------
	img (numpy array):
		2D image like the input, but with Bayer filter scaled values.
	"""
	# Red:
	img[1::2,1::2] *= flux[0]
	# Green:
	img[0::2,1::2] *= flux[1]
	img[1::2,0::2] *= flux[1]
	# Blue:
	img[0::2,0::2] *= flux[2]
	return img



def mag2flux(mag):
	"""
	Convert from magnitude to flux using scaling relation from
	aperture photometry for the TESS satellite. This is an estimate.

	Parameters
	----------
	mag (float): 
		Magnitude in TESS band.

	Returns
	-------
	float: 
		Corresponding flux value.
	"""
	return 10**(-0.4*(mag - 28.24))

