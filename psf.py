#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PSF methods for the delsimi simulation code.

@author: Jonas Svenstrup Hansen, jonas.svenstrup@gmail.com
"""

import numpy as np
from scipy.special import erf
from scipy.signal import convolve2d
from scipy.interpolate import RectBivariateSpline
from skimage.draw import line_aa

from utilities import make_bayer_filter

class PSF():
	def __init__(self, imshape, superres = 10, rgb_filter = 'bggr'):
		"""
		Approximate the PSF in stellar images from Delphini-1.
		
		Future Extensions
		-----------------
		Linear filters for convolution with smear and jitter:
		They posses additivity (L[f+g] = L[f] + L[g]) and homogeneity (L[cf] =
		cL[f], where c is a constant) according to Bogges and Narcowich 2009, 
		p. 110.
		
		Use the following sources for ideas on applying smear and jitter:
		https://www.osapublishing.org/ao/fulltext.cfm?uri=ao-32-32-6503&id=40363
		https://github.com/TESScience/SPyFFI/blob/master/PSF.py
		https://github.com/TESScience/SPyFFI/blob/master/Jitter.py
		
		Code Author
		-----------
		Jonas Svenstrup Hansen, jonas.svenstrup@gmail.com
		"""
		self.imshape = np.array(imshape) # shape of image in pixels (row, col)
		self.superres = superres # subpixel resolution
		self.rgb_filter = rgb_filter



	def integrate_to_image(self, stars, integration_time, angle_vel, 
			speed = None, fwhm = 1., jitter = False, focus = False):
		"""
		Integrate a PSF that is smeared in one direction to an image.

		Parameters
		----------
		stars (list):
			List with an element for each star. Each element contains the
			elements ``[row, col, [flux_R, flux_G, flux_B]]`` which are used
			to generate the star. The row and column position corresponds
			to the star position at the midtime of exposure.
		integration_time (float):
			CCD integration time.
		angle_vel (float):
			Angle in radians of star CCD movement.
		speed (float):
			Speed of a star in pixels. Default is ``Ǹone```which yields the 
			standard speed estimated from the speed of the ISS.
		fwhm (float):
			Full width at half maximum of PSF in pixels. Default is ``1.``.
		jitter (string):
			``True`` if jitter is to be applied. Default is ``False``. Not
			implemented.
		focus (string):
			``True`` if focus is to be applied. Default is ``False``. Not
			implemented.

		Returns
		-------
		img (2D array, float):
			Smeared, pixel-integrated and Bayer filter scaled PSF.
		smearKernel (2D array, float):
			Kernel with a line that specifies the smear of a single PSF.
		PSFhighres (2D array, float):
			Subpixel resolution PSF that is convolved with the smear kernel.
		highresConvPSF (2D array, float):
			Normalised convolution of smear kernel and subpixel resolution PSF.
		highresImageInterp (interpolation object):
			Interpolated smeared PSF. Can be integrated efficiently.
		"""
		# Set speed to standard speed if not given as parameter:
		if speed is None:
			speed = self.speed

		# Define subpixel buffer. Needs to be large for correct interpolation:
		self.buffer = np.int(3*fwhm*self.superres)

		# Create smear kernel:
		smearKernel, r0, c0, r1, c1 = self.makeSmearKernel(
				integration_time, angle_vel, speed, fwhm)
		self.kernelShape = smearKernel.shape

		# Get highres PSF:
		PSFhighres = self.highresPSF(fwhm)

		# TODO: convolve highres PSF with focus and jitter here

		# Convolve the PSF with the smear kernel:
		highresConvPSF = self.convolvePSF(PSFhighres, smearKernel)

		# Normalise the PRF:
		highresConvPSF /= np.nansum(highresConvPSF) * self.superres**2

		# Define pixel centered index arrays for the interpolater:
		PRFrow = np.arange(0.5, self.kernelShape[0] + 0.5)
		PRFcol = np.arange(0.5, self.kernelShape[1] + 0.5)

		# Center around 0:
		PRFrow = PRFrow - PRFrow.size / 2
		PRFcol = PRFcol - PRFcol.size / 2

		# Convert from subpixel to pixel resolution:
		PRFrow /= self.superres
		PRFcol /= self.superres

		# Interpolate highresImage:
		highresImageInterp = RectBivariateSpline(PRFrow, PRFcol, highresConvPSF)

		# Preallocate image array:
		img = np.zeros(self.imshape, dtype=np.float64)

		# Prepare Bayer filter scaling:
		Bayer_filter = make_bayer_filter(img.shape)

		# Integrate the interpolation object in each pixel:
		for star in stars:
			for row in range(self.imshape[0]):
				for col in range(self.imshape[1]):
					# Get star position in PSF(t=0)-based coordinates:
					row_cen = row - star[0]
					col_cen = col - star[1]
					# Integrate only significant contributions to avoid artefacts:
					withinBoundary = highresImageInterp(row_cen, col_cen) > 1e-9
					if withinBoundary:
						# Get flux value:
						# Red:
						if Bayer_filter[row,col] == 0:
							Bayer_flux = star[2][0]
						# Green:
						elif Bayer_filter[row,col] == 1:
							Bayer_flux = star[2][1]
						# Blue:
						elif Bayer_filter[row,col] == 2:
							Bayer_flux = star[2][2]
						else:
							raise ValueError(
									'Bayer filter flag must be 0, 1 or 2.')
						# Integrate normalised interpolation in the current pixel:
						img[row,col] += Bayer_flux * integration_time * \
								highresImageInterp.integral(
									row_cen-0.5, row_cen+0.5,
									col_cen-0.5, col_cen+0.5)

		return img, smearKernel, PSFhighres, highresConvPSF, highresImageInterp



	def makeSmearKernel(self, integration_time, angle_vel, speed, fwhm):
		"""
		Make a smear kernel that describes the large-scale movement of a star
		on the CCD. Do this by making a high resolution line that approximates
		the movement of a star with a certain constant angle and speed across 
		the CCD during the integrationTime.
		
		The buffer is the area around the kernel corresponding to the extent of
		the PSF. Since the smear kernel array must have the same shape as the 
		PSF array, this buffer is applied here too.
		
		Note that the current implementation only includes line start at 
		integer subpixel positions. With sufficient oversampling this effect is
		not expected to change the final star positions significantly.
		
		Parameters
		----------
		integration_time (float)
			CCD integration time.
		angle_vel (float)
			Angle in radians of star CCD movement.
		speed (float)
			Speed of star CCD movement.
		fwhm (float)
			Full width at half maximum of PSF. Used to determine buffer pixels
			around the line.
		
		Returns
		-------
		(array, float)
			Interpolated pixel positions of a star.
		"""
		
		# Set buffer around the line in subpixel units:
		buffer = self.buffer
		
		if speed is 0:
			out = np.zeros([2*buffer+1, 2*buffer+1])
			out[buffer,buffer] = 1.
			r0 = buffer
			c0 = buffer
			r1 = buffer
			c1 = buffer
		else:
			# Get integer pixel final position coordinates:
			finalposRow = np.int(
					self.superres*speed*integration_time*np.sin(angle_vel))
			finalposCol = np.int(
					self.superres*speed*integration_time*np.cos(angle_vel))
			
			# Set starting point (r0,c0):
			r0 = 0 + buffer -1
			c0 = 0 + buffer -1
			
			# Set ending point (r1,c1):
			r1 = r0 + finalposRow
			c1 = c0 + finalposCol
			
			# Preallocate array in which to draw line:
			out = np.zeros([r1+buffer, c1+buffer])
			
			# Draw line:
			# TODO: change line implementation to subsubpixel line definition
			# Non-anti-aliased:
	#		rr, cc = line(r0, c0, r1, c1)
	#		out[rr, cc] = 1
			# Anti-aliased:
			rr, cc, val = line_aa(r0, c0, r1, c1)
			out[rr, cc] = val
			
		return out, r0, c0, r1, c1



	def makeJitterKernel(self):
		# TODO: dependent on satellite pitch, yaw, roll
		# Make spline of parametric curve:
#		splprep([positions[:,1],positions[:,0]], s=0)
		pass



	def makeFocusKernel(self):
		# TODO: start out by making this a constant blur
		pass



	def convolvePSF(self, PSFunconvolved, kernel):
		"""
		Inspired by populateJitteredPSFLibrary in 
		
		https://github.com/TESScience/SPyFFI/blob/master/PSF.py
		
		Use only the subpixel-resolution PSF for this!
		
		Parameters
		----------
		PSFunconvolved (numpy array): 
			2D array with the subpixel PSF.
		kernel (numpy array): 
			Kernel with which to convolute.
		"""
		return convolve2d(PSFunconvolved, kernel, 'same', 'fill', 0)



	def highresPSF(self, fwhm):
		"""
		Make a subpixel resolution PSF.
		
		The PSF is placed at the center of the high resolution array.
		
		Parameters
		----------
		fwhm (float)
			Full width at half maximum of Gaussian PSF in pixels.
		
		Returns
		-------
		(2D array, float)
			Subpixel-resolution PSF.
		"""
		# Make subpixel resolution grid with superres oversampling:
		X, Y = np.meshgrid(np.arange(0,self.kernelShape[1]), 
						np.arange(0,self.kernelShape[0]))
		
		# Define centroid position:
		x_0 = self.kernelShape[1]/2-0.5
		y_0 = self.kernelShape[0]/2-0.5
		
		# Evaluate PSF on grid with superres increased width:
		return self.integratedGaussian(X, Y, 1, x_0, y_0, 
					sigma = self.superres*fwhm/(2*np.sqrt(2*np.log(2))))



	def integratedGaussian(self, x, y, flux, x_0, y_0, sigma=1):
		'''
		Evaluate a 2D symmetrical Gaussian integrated in pixels.
		
		Parameters
		----------
		x (numpy array): 
			x coordinates at which to evaluate the PSF.
		y (numpy array): 
			y coordinates at which to evaluate the PSF.
		flux (float): 
			Integrated value.
		x_0 (float): 
			Centroid position.
		y_0 (float): 
			Centroid position.
		sigma (float, optional): 
			Standard deviation of Gaussian. Default=1.
		
		Returns
		-------
		numpy array: 
			2D Gaussian integrated pixel values at (x,y).
		
		Example
		-------
		
		>>> import numpy as np
		>>> X, Y = np.meshgrid(np.arange(-1,2), np.arange(-1,2))
		>>> integratedGaussian(X, Y, 10, 0, 0)
		array([[ 0.58433556,  0.92564571,  0.58433556],
			[ 0.92564571,  1.46631496,  0.92564571],
			[ 0.58433556,  0.92564571,  0.58433556]])
		
		Code Author
		-----------
		Documentation: Jonas Svenstrup Hansen, jonas.svenstrup@gmail.com
		
		Code: 
		https://github.com/astropy/photutils/blob/master/photutils/psf/models.py
		(IntegratedGaussianPRF)
		'''
		return (flux / 4 *
			((erf((x - x_0 + 0.5) / (np.sqrt(2) * sigma)) -
			  erf((x - x_0 - 0.5) / (np.sqrt(2) * sigma))) *
			 (erf((y - y_0 + 0.5) / (np.sqrt(2) * sigma)) -
			  erf((y - y_0 - 0.5) / (np.sqrt(2) * sigma)))))




if __name__ == '__main__':
	import matplotlib.pyplot as plt
	from utilities import mag2flux
	import string

	# Make PSF class instance:
	dpsf = PSF(imshape=[50,80], superres=10)

	# Evaluate PSF with specified parameters:
	stars = [[20.,20.,[mag2flux(3.), mag2flux(4.), mag2flux(5.)]],
			[27.,50.,[mag2flux(5.5), mag2flux(4.5), mag2flux(3.5)]]]
	integration_time = 3 # seconds
	angle_vel = np.pi/5.7 # radians
	speed = 10 # pixel per second
	fwhm = 1.5 # pixel
	img, smearKernel, PSFhighres, highresImage, highresImageInterp = \
		dpsf.integrate_to_image(
			stars=stars, 
			integration_time=integration_time, 
			angle_vel=angle_vel, 
			speed=speed, 
			fwhm=fwhm)

	# Plot:
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)

	ax1.imshow(img, origin='lower')
	ax1.set_xlabel('Pixel column')
	ax1.set_ylabel('Pixel row')
	ax1.set_title('Pixel-integrated image')

	ax2.imshow(smearKernel, origin='lower')
	ax2.set_title('Smear kernel')

	ax3.imshow(highresImage, origin='lower')
	ax3.set_title('High res. convolved PSF')

	ax4.imshow(PSFhighres, origin='lower')
	ax4.set_title('High resolution PSF')

	for ax in (ax2, ax3, ax4):
		ax.set_xlabel('Subpixel column')
		ax.set_ylabel('Subpixel row')

	for n,ax in enumerate((ax1, ax2, ax3, ax4)):
		ax.text(0.05, 0.75, string.ascii_uppercase[n], transform=ax.transAxes, 
			size=20, weight='bold', color='w')

	# Make space for the subplot titles:
	fig.subplots_adjust(hspace=0.5, wspace=0.5)
