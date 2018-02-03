#import matplotlib.pyplot as plt
#import math as mt
#import matplotlib
#from scipy import misc
from PIL import Image
import os
import numpy as np
import scipy
from psf import PSF


class delsimi(object):
	def __init__(self, integrationTime = 30.0):
		"""
		Simulate stellar images from Delphini-1.
		
		Input
		-----
		integrationTime (float)
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
		
		MTF
		---
		# TODO: is this correct?
		The Modulation transfer function (MTF) is a description of both jitter
		(high temporal frequency compared to exposure time) and smear (low
		temporal frequency compared to exposure time). It ferciliates the 
		simulation of movement during CCD integration. Thus the smear caused
		by the orbital motion of the satellite as well as possible jitter can
		be included in the simulation.
		NOTE: This is currently not implemented.
		
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
		
		self.infiledir = "../infiles"
		self.outfiledir = "../outfiles"
		
		# Set size in pixels to Delphini's CCD size, Aptina MT9T031 1/2" (4:3):
		self.ccdshape = np.array([1536,2048])
		
		# Set pixel size in arcseconds:
		self.pixel_scale = 42 # TODO: get correct pixel scale
		
		
		# Load catalog here:
		
		
		# Instantiate PSF class:
		dpsf = PSF(imshape=self.ccdshape, superres=10)
		
		# Evaluate PSF class:
		star_row = 20
		star_col = 20
		star_flux = 1e3
		star = [star_row, star_col, star_flux]
		img, smearKernel, PSFhighres, highresConvPSF, highresImageInterp = \
			dpsf.evaluate(star=star, integrationTime=self.integrationTime,
				angle=angle, speed=speed, fwhm=fwhm)
		
		# Apply binning of Bayer pixels following the method linked to below:
		# https://stackoverflow.com/questions/14916545/numpy-rebinning-a-2d-array
		row_bin = 2
		col_bin = 2
		img_view = img.reshape(img.shape[0] // row_bin, row_bin,
							img.shape[1] // col_bin, col_bin)
		img_binned = img_view.sum(axis=3).sum(axis=1)



if __name__ == '__main__':
	# Make delsimi class instance:
	simtest = delsimi()
	
	# Make PSF class instance:
	# TODO: do not make a PSF instance here, do it in the delsimi class __init__
	dpsf = PSF(imshape=simtest.jpgfile.size, superres=5)
	
	# Evaluate PSF with specified parameters:
	img, smearKernel, PSFhighres, highresImage, highresImageInterp = dpsf.evaluate(
			star=[10,15], integrationTime=10, angle=np.pi/7, speed=3, fwhm=1)
	
	