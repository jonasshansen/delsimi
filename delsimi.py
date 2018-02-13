#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Delsimi image simulation code.

@author: Jonas Svenstrup Hansen, jonas.svenstrup@gmail.com
"""

import os
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

from utilities import rvb2rgb, star_CCD_speed, mag2flux
from psf import PSF

class delsimi(object):
	def __init__(self, input_dir='../infiles',
					output_dir='../outfiles',
					overwrite=True,
					coord_cen=[0.,0.],
					integration_time=0.1,
					angle_vel=0.,
					angle_sat=0.):
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
		overwrite (boolean):
			``True`` if to overwrite FITS images in output.
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
		Jonas Svenstrup Hansen, jonas.svenstrup@gmail.com (extensions, see all
		but the first commit on Github)

		Carolina von Essen, cessen@phys.au.dk (initial idea, see first commit on
		Github for code)

		Online Repository
		-----------------
		https://github.com/jonasshansen/delsimi
		"""


		""" Set constants """
		self.input_dir = input_dir
		self.output_dir = output_dir
		self.overwrite = overwrite

		self.coord_cen = coord_cen
		self.integration_time = float(integration_time)
		self.angle_vel = float(angle_vel)
		self.angle_sat = float(angle_sat)

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
		# TODO: add limitation of catalog stars from position and velocity

		# Generate astropy table and WCS solution:
		catalog, w = self.make_catalog(cat_input=None)


		""" Make image """
		# Calculate average speed of a star on the CCD in pixel units:
		speed = star_CCD_speed(self.pixel_scale)

		# Instantiate PSF class:
		dpsf = PSF(imshape=self.ccd_shape, superres=10)

		# Prepare list with star information:
		stars = [[row, col, [flux_R, flux_G, flux_B]] for 
				(row, col, flux_R, flux_G, flux_B) in 
				zip(catalog['row'], catalog['col'], \
					catalog['flux_R'], catalog['flux_G'], catalog['flux_B'])]

		# Integrate stars to image:
		img, smearKernel, PSFhighres, highresConvPSF, highresImageInterp = \
			dpsf.integrate_to_image(stars=stars,
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

		# Apply noise to image:
		img += dark_current + read_noise


		""" Convert from electrons to digital units """
		img /= self.gain


		""" Apply 2x2 pixel binning """
		# Apply sum binning of Bayer pixels using the method linked to below:
		# https://stackoverflow.com/questions/14916545/numpy-rebinning-a-2d-array
		row_bin = 2
		col_bin = 2
		row_col_bin = np.array([row_bin, col_bin])
		img_view = img.reshape(img.shape[0] // row_bin, row_bin,
							img.shape[1] // col_bin, col_bin)
		img_binned = img_view.sum(axis=3).sum(axis=1)

		# Update the WCS solution:
		w.wcs.crpix = np.divide(w.wcs.crpix, row_col_bin)
		w.wcs.cdelt = np.divide(w.wcs.cdelt, row_col_bin)

		# Update catalog positions
		catalog['row'] /= row_bin
		catalog['col'] /= col_bin


		""" Export binned image to fits """
		# Make primary header data unit:
		header = w.to_header()
		hdu = fits.PrimaryHDU(data=img_binned, header=header)
		hdu.header['NAXIS'] = (2, 'Number of data dimension')
		hdu.header['NAXIS1'] = (np.shape(img_binned)[1], 'Number of pixel columns')
		hdu.header['NAXIS2'] = (np.shape(img_binned)[0], 'Number of pixel rows')
		
		# Write image to fits file:
		image_output_fname = os.path.join(self.output_dir, 'image.fits')
		print('Writing image to '+image_output_fname)
		hdu.writeto(image_output_fname, overwrite=self.overwrite)


		""" Export catalog to ASCII file """
		# Save catalog to file:
		catalog_output_fname = os.path.join(self.output_dir, 'catalog.txt')
		print('Writing catalog to '+catalog_output_fname)
		print(catalog)
		np.savetxt(catalog_output_fname,
			np.asarray(catalog),
			delimiter='\t',
			header='    '.join(catalog.colnames))



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
			 - flux_R (int): Bayer filter R flux of stars in electrons.
			 - flux_G (int): Bayer filter G flux of stars in electrons.
			 - flux_B (int): Bayer filter B flux of stars in electrons.
		w (astropy WCS solution):
			World Coordinate System solution for the current position. Can be
			used to transform from ``(ra,dec)`` to pixel ``(row,col)`` with 
			``row_col = w.wcs_world2pix([RA, DEC], 0)``.
		"""

		if cat_input is None:
			# Generate two test stars:
			cat_input = 	[np.array([0, 1]),     # starid
				np.array([self.coord_cen[0], 
						self.coord_cen[0] + 
						10/self.pixel_scale]),# rigth ascension
				np.array([self.coord_cen[1], 
						self.coord_cen[0] + 
						5/self.pixel_scale]), # declination
				np.array([4., 6.5]),           # R magnitude (Cousins)
				np.array([5., 5.5]),           # V magnitude (Johnson)
				np.array([6., 4.5])]           # B magnitude (Johnson)

		# Extract parameters from catalog input:
		starid, RA, DEC, R_Cousins, V_Johnson, B_Johnson = cat_input

		# WCS initialisation:
		w = WCS(naxis=2)
		w.wcs.crpix = self.ccd_shape[::-1]/2
		w.wcs.cdelt = [self.pixel_scale/3600, self.pixel_scale/3600]
		w.wcs.crval = self.coord_cen
		w.wcs.ctype = ["RA---AIR", "DEC--AIR"]

		# WCS conversion from (ra,dec) to pixel (row,col):
		col_row = w.wcs_world2pix(np.transpose(np.array([RA, DEC])), 0)
		CCD_row = np.array([coord[1] for coord in col_row])
		CCD_col = np.array([coord[0] for coord in col_row])

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
					'flux_R', 'flux_G', 'flux_B'),
			dtype=('int64', 'float64', 'float64', 'float64', 'float64',
					'float32', 'float32', 'float32',
					'float32', 'float32', 'float32',
					'int64', 'int64', 'int64')
		)
		
		return catalog, w



if __name__ == '__main__':
	# Make delsimi class instance:
	simtest = delsimi(angle_vel=np.pi/7)
