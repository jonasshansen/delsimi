#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

from utilities import uvb2rgb
from psf import PSF

class delsimi(object):
	def __init__(self, input_dir='../infiles',
					output_dir='../outfiles',
					coord_cen=[0.,0.],
					integrationTime=30.0):
		"""
		Simulate stellar images from Delphini-1.

		It is assumed that the camera pointing is orthogonal to the tangent of
		the orbit.

		Parameters
		----------
		input_dir (string):
			Input file directory. Default is ``'../infiles'``.
		output_dir (string):
			Output file directory. Default is ``'../outfiles'``.
		coord_cen (list of floats):
			Right ascension and declination coordinate center at the midtime of
			exposure.
		integrationTime (float):
			CCD integration time in seconds. Default is 30 seconds.

		Future Extensions
		-----------------
		This is a list of possible future extensions to this code. The list is
		largely inspired by the SPyFFI manual which is linked to here:
			
			https://github.com/TESScience/SPyFFI
			
		- Additional stars
		- PSF length inferred from estimated satellite speed and exposure time
		- True pixel sizes
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

		self.integrationTime = integrationTime

		self.input_dir = input_dir
		self.output_dir = output_dir

		# Set size in pixels to Delphini's CCD size, Aptina MT9T031 1/2" (4:3):
		self.ccd_shape = np.array([1536,2048])

		# Set pixel scale in arcseconds:
		pixel_size = 3.2 # micrometer
		focal_length = 35 # millimeter
		conv_const = 206.265 # convert micrometer to mm and radians to arcsec
		self.pixel_scale = pixel_size / focal_length * conv_const # arcsec/pixel

		# TODO: Load catalog here:
		

		# TODO: WCS conversion from (ra,dec) to pixel (row,col):

		# Collect star parameters in list for catalog:
		cat = [starids, starrows, starcols, mag_b, mag_v, mag_r]

		# Make astropy table with catalog:
		return Table(
			cat,
			names=('starid', 'row', 'col', 'mag_b', 'mag_v', 'mag_r'),
			dtype=('int64', 'float64', 'float64', 'float32', 'float32', 'float32')
		)

		for star in cat:
			# Convert Johnson filters to RGB colors:
			rgb = bvr2rgb(bvr)
	#		flux_r, flux_b, flux_g = bvr2rgb(np.array([flux_b, flux_v, flux_r]))
			star_flux = rgb
	
			# Convert magnitudes to flux:
			mag2flux(star_flux)
			# TODO: requires absolute scaling only obtainable from photograph
			# TODO: add rgb fluxes to catalog

		# Instantiate PSF class:
		dpsf = PSF(imshape=self.ccd_shape, superres=10)

		# Evaluate PSF class:
		star_row = 20
		star_col = 20
		star_flux = 1e3
		star = [star_row, star_col, star_flux]
		img, smearKernel, PSFhighres, highresConvPSF, highresImageInterp = \
			dpsf.evaluate(star=star, integrationTime=self.integrationTime,
				angle=angle, speed=speed, fwhm=fwhm)

		# TODO: Add white noise to image:

		# Apply sum binning of Bayer pixels using the method linked to below:
		# https://stackoverflow.com/questions/14916545/numpy-rebinning-a-2d-array
		row_bin = 2
		col_bin = 2
		img_view = img.reshape(img.shape[0] // row_bin, row_bin,
							img.shape[1] // col_bin, col_bin)
		img_binned = img_view.sum(axis=3).sum(axis=1)

		# Save data to fits file:
#		# Make WCS solution parameters:
#		w = WCS(naxis=2)
#		w.wcs.crpix = [0,0]
#		w.wcs.cdelt = [self.pixel_scale/3600, self.pixel_scale/3600]
#		w.wcs.crval = self.coord_zero_point # [0.,0.]
#		w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
#		header = w.to_header()
#		
#		# Instantiate primary header data unit:
#		hdu = fits.PrimaryHDU(data=img, header=header)
#		
#		# Add timestamp to header with a unit of days:
#		hdu.header['BJD'] = (timestamp/3600/24, 
#			'time in days (arb. starting point)')
#		hdu.header['NAXIS'] = (2, 'Number of data dimension')
#		hdu.header['NAXIS1'] = (self.Ncols, 'Number of pixel columns')
#		hdu.header['NAXIS2'] = (self.Nrows, 'Number of pixel rows')
#		# TODO: write more info to header
#		
#		
#		# Specify output directory:
#		if outdir is None:
#			outdir = os.path.join(self.output_folder, 'images')
#		
#		# Remove any previous hdf5 file made by prepare_photometry:
#		try:
#			hdf5filename = 'camera1_ccd1.hdf5'
#			os.remove(os.path.join(self.output_folder,hdf5filename))
#		except:
#			pass
#		
#		# Write FITS file to output directory:
#		hdu.writeto(os.path.join(outdir, 'test%02d.fits' % i),
#					overwrite=self.overwrite_images)

		# Save catalog to file:
		# TODO


if __name__ == '__main__':
	# Make delsimi class instance:
	simtest = delsimi()
	
	# Make PSF class instance:
	# TODO: do not make a PSF instance here, do it in the delsimi class __init__
	dpsf = PSF(imshape=simtest.jpgfile.size, superres=5)
	
	# Evaluate PSF with specified parameters:
	img, smearKernel, PSFhighres, highresImage, highresImageInterp = dpsf.evaluate(
			star=[10,15], integrationTime=10, angle=np.pi/7, speed=3, fwhm=1)
	
	