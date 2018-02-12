#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Delsimi image simulation code.

@author: Jonas Svenstrup Hansen, jonas.svenstrup@gmail.com
"""

#import matplotlib.pyplot as plt
#import math as mt
#import matplotlib
#from scipy import misc
from PIL import Image
import os
import numpy as np
import scipy
from astropy.io import fits
from astropy.table import Table, Column
from astropy.wcs import WCS

from utilities import rvb2rgb, star_CCD_speed, mag2flux
from psf import PSF

class delsimi(object):
	def __init__(self, input_dir='../infiles',
					output_dir='../outfiles',
					coord_cen=[0.,0.],
					integration_time=30.,
					angle_vel=0.,
					angle_sat=0.,):
		"""
		Simulate stellar images from Delphini-1.

		It is assumed that the camera pointing is orthogonal to the tangent of
		the orbit.

		The output images will be in digital units.

		Parameters
		----------
		input_dir (string):
			Input file directory. Default is ``'../infiles'``.
		output_dir (string):
			Output file directory. Default is ``'../outfiles'``.
		coord_cen (list of floats):
			Right ascension and declination coordinate center at the midtime of
			exposure.
		integration_time (float):
			CCD integration time in seconds. Default is ``30.``.
		angle_vel (float):
			Angle of velocity in radians with respect to the ecliptic 
			coordinates. Default is ``0.``.
		angle_sat (float):
			Angle of satellite in radians with respect to the ecliptic 
			coordinates. Default is ``0.``.
		

		Future Extensions
		-----------------
		This is a list of possible future extensions to the code. The list is
		largely inspired by the SPyFFI manual which is linked to here:
			
			https://github.com/TESScience/SPyFFI
			
		- Additional stars
		- PSF length inferred from estimated satellite speed and exposure time
		- Catalog stars
		- Simple noise
		- Focus change
		- Jitter
		- Pixel sensitivity variations
		- Realistic constant used for magnitude to flux conversion (color dependent)
		- Realistic Johnson-Cousins BVR to RGB magnitude conversion constants
		- Background
		- Photon noise
		- Readout smear
		- Rolling shutter
		- Saturated pixels
		- Cosmic rays
		- Position variable PSF (not compatible with convolution)

		PSF
		---
		For the point spread function (PSF) a Gaussian is applied. This is an
		approximation.

		Sensor
		------
		The Delphini-1 sensor is an Aptina MT9T031 1/2" (4:3) CMOS Bayer 
		pattern sensor. Thus, each pixel has a colour filter. Binning is 
		applied to convert the Bayer pattern pixels to greyscale. Four pixels,
		red, blue and two green, are summed to yield one greyscale pixel. This
		binning is implemented in the code, assuming no cromatic abberation, 
		which changes the PSF based on wavelength.


		Code Authors
		------------
		Carolina von Essen, cessen@phys.au.dk (base code, see first commit on
		Github)

		Jonas Svenstrup Hansen, jonas.svenstrup@gmail.com (extensions, see all
		but the first commit on Github)

		Online Repository
		-----------------
		https://github.com/jonasshansen/delsimi
		"""



		""" Set constants """
		self.input_dir = input_dir
		self.output_dir = output_dir

		self.coord_cen = coord_cen
		self.integration_time = integration_time
		self.angle_vel = angle_vel
		self.angle_sat = angle_sat

		# Set gain:
		self.gain = 25 # electrons per LSB or ADU

		# Set size in pixels to Delphini's CCD size, Aptina MT9T031 1/2" (4:3):
		self.ccd_shape = np.array([1536,2048])

		# Set pixel scale in arcseconds:
		pixel_size = 3.2 # micrometer
		focal_length = 35 # millimeter
		conv_const = 206.265 # convert micrometer to mm and radians to arcsec
		self.pixel_scale = pixel_size / focal_length * conv_const # arcsec/pixel


		""" Load Catalog """
		# TODO: add input catalog loading here
		catalog, w = self.make_catalog(cat_input=None)


		""" Make image """
		# Calculate average speed of a star on the CCD in pixel units:
		speed = star_CCD_speed(self.pixel_scale)

		# Instantiate PSF class:
		dpsf = PSF(imshape=self.ccd_shape, superres=10)

		# Evaluate PSF class:
		star_row = 20
		star_col = 20
		star_flux = 1e3
		star = [star_row, star_col, star_flux]
		img, smearKernel, PSFhighres, highresConvPSF, highresImageInterp = \
			dpsf.integrate_to_image(star=star,
				integration_time=self.integration_time,
				angle_vel=self.angle_vel, speed=speed, fwhm=1.)


		""" Make noise """
		# Set random number generator seed to a random state:
		np.random.seed(seed=None)

		# Generate dark current:
		dark_current = self.integration_time * np.random.normal(
						loc=1e2, scale=1e1, size=img.shape)

		# Generate read noise:
		read_noise = self.integration_time * np.random.normal(
						loc=0., scale=6., size=img.shape)

		# White noise to simulate background, faint stars and other sources:
		other_noise = np.random.normal(loc=1e3, scale=1e5, size=img.shape)

		# Apply noise to image:
		img += dark_current + read_noise + other_noise


		""" Apply binning """
		# Apply sum binning of Bayer pixels using the method linked to below:
		# https://stackoverflow.com/questions/14916545/numpy-rebinning-a-2d-array
		row_bin = 2
		col_bin = 2
		img_view = img.reshape(img.shape[0] // row_bin, row_bin,
							img.shape[1] // col_bin, col_bin)
		img_binned = img_view.sum(axis=3).sum(axis=1)

		# Update the WCS solution:
		


		""" Convert from electrons to digital units """
		img_binned /= self.gain


		""" Export to fits """
		# Instantiate primary header data unit:
		header = w.to_header()
		hdu = fits.PrimaryHDU(data=img_binned, header=header)
		
		# Save data to fits file:
		#img + noise


		""" Export catalog to ASCII file """
		# Save catalog to file:
		# TODO



	def make_catalog(self, cat_input=None):
		"""
		Generate an astropy Table with the catalog from stellar position and
		magnitudes.

		Parameters
		----------
		cat_input (list of arrays):
			List with numpy arrays containing the following:
			 - starid (int): Star id for internal use, e.g. ``0``, ``1``, ...
			 - RA (float): Right ascension of stars.
			 - DEC (float): Declination of stars.
			 - R (float): Cousins filters R magnitude of stars.
			 - V (float): Johnson filters V magnitude of stars.
			 - B (float): Johnson filters B magnitude of stars.
			Default is ``None`` which generates an ``cat_input`` with two test
			stars with the following parameters:
			 - starid (int): np.array([0, 1])
			 - RA (float): np.array([20., 20.01])
			 - DEC (float): np.array([20., 20.015])
			 - R (float): np.array([4., 6.5])
			 - V (float): np.array([5., 5.5])
			 - B (float): np.array([6., 4.5])

		Returns
		-------
		catalog (astropy Table): 
			Table with catalog information. The columns contain the following:
			 - starid (int): Star id for internal use, e.g. ``0``, ``1``, ...
			 - RA (float): Right ascension of stars.
			 - DEC (float): Declination of stars.
			 - row (float): Pixel row position on CCD.
			 - col (float): Pixel column position on CCD.
			 - R_Cousins (float): Cousins filters R magnitude of stars.
			 - V_Johnson (float): Johnson filters V magnitude of stars.
			 - B_Johnson (float): Johnson filters B magnitude of stars.
			 - R_Bayer (float): Bayer filter R magnitude of stars.
			 - G_Bayer (float): Bayer filter G magnitude of stars.
			 - B_Bayer (float): Bayer filter B magnitude of stars.
			 - flux_R (float): Bayer filter R flux of stars.
			 - flux_G (float): Bayer filter G flux of stars.
			 - fluxB (float): Bayer filter B flux of stars.
		w (astropy WCS solution):
			World Coordinate System solution for the current position. Can be
			used to transform from ``(ra,dec)`` to pixel ``(row,col)`` with 
			``row_col = w.wcs_world2pix([RA, DEC], 0)``.
		"""

		if cat_input is None:
			# Generate two test stars:
			cat_input = 	[np.array([0, 1]),   # starid
				np.array([20., 20.01]),      # rigth ascension
				np.array([20., 20.015]),     # declination
				np.array([4., 6.5]),         # R magnitude (Cousins)
				np.array([5., 5.5]),         # V magnitude (Johnson)
				np.array([6., 4.5])]         # B magnitude (Johnson)

		# Extract parameters from catalog input:
		starid, RA, DEC, R_Cousins, V_Johnson, B_Johnson = cat_input

		# WCS initialisation:
		# FIXME: binning and wcs solution problem
		w = WCS(naxis=2)
		w.wcs.crpix = self.ccd_shape/2
		w.wcs.cdelt = [self.pixel_scale/3600, self.pixel_scale/3600]
		w.wcs.crval = self.coord_cen
		w.wcs.ctype = ["RA---AIR", "DEC--AIR"]

		# WCS conversion from (ra,dec) to pixel (row,col):
		row_col = w.wcs_world2pix([RA, DEC], 0)
		CCD_row = np.array([coord[0] for coord in row_col])
		CCD_col = np.array([coord[1] for coord in row_col])

		# Collect star parameters in list for catalog:
		cat = [starid, RA, DEC, CCD_row, CCD_col,
			 R_Cousins, V_Johnson, B_Johnson]

		# Convert Johnson-Cousins filter magnitudes to RGB magnitudes:
		# TODO: apply accurate filter transformation
		R_Bayer, G_Bayer, B_Bayer = rvb2rgb([R_Cousins, V_Johnson, B_Johnson])
		for M_Bayer in [R_Bayer, G_Bayer, B_Bayer]:
			cat.append(M_Bayer)

		# Convert magnitudes to flux:
		# TODO: get constants for accurate magnitude to flux transformation
		for M_Bayer, M_type in zip([R_Bayer, G_Bayer, B_Bayer],['R', 'G', 'B']):
			cat.append(mag2flux(M_Bayer, M_type))

		# Make astropy table with catalog:
		catalog = Table(
			cat,
			names=('starid', 'RA', 'DEC', 'row', 'col', 
					'R_Cousins', 'V_Johnson', 'B_Johnson',
					'R_Bayer', 'G_Bayer', 'B_Bayer',
					'flux_R', 'flux_G', 'fluxB'),
			dtype=('int64', 'float64', 'float64', 'float64', 'float64',
					'float32', 'float32', 'float32',
					'float32', 'float32', 'float32',
					'float64', 'float64', 'float64')
		)
		
		return catalog, w



if __name__ == '__main__':
	# Make delsimi class instance:
	simtest = delsimi()
#	
#	# Make PSF class instance:
#	# TODO: do not make a PSF instance here, do it in the delsimi class __init__
#	dpsf = PSF(imshape=simtest.jpgfile.size, superres=5)
#	
#	# Evaluate PSF with specified parameters:
#	img, smearKernel, PSFhighres, highresImage, highresImageInterp = dpsf.evaluate(
#			star=[10,15], integrationTime=10, angle=np.pi/7, speed=3, fwhm=1)
	
	