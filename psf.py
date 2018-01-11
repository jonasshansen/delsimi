import numpy as np
from scipy.special import erf
from scipy.signal import convolve2d
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import splprep
from skimage.draw import line

class PSF():
	def __init__(self, imshape, superres = 10):
		"""
		Approximate the PSF in stellar images from Delphini-1.
		
		Future Extensions
		-----------------
		Linear filters:
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

	def evaluate(self, starpos, integrationTime, angle, speed, fwhm = 2, 
			jitter = False, focus = False):
		"""
		Evaluate a PSF that is smeared in one direction.
		
		Input
		-----
		starpos (array, float)
			Row and column position in pixels of the star.
		integrationTime (float)
			CCD integration time.
		angle (float)
			Angle in radians of star CCD movement.
		speed (float)
			Speed of star CCD movement.
		fwhm (float)
			Full width at half maximum of PSF in pixels.
		
		Output
		------
		(2D array, float)
			Smeared and pixel-integrated PSF.
		"""
		# Interpolate positions to subpixel grid:
		positionsKernel = self.makeSmearKernel(starpos, integrationTime, 
										angle, speed)
		
		# Get highres PSF:
		PSFhighres = self.highresPSF(fwhm)
		
		# TODO: convolve highres PSF with focus and jitter here
		
		# Convolve the PSF with the interpolated positions:
		highresImage = self.convolvePSF(PSFhighres, positionsKernel)
		
		# Define pixel centered index arrays for the interpolater:
		PRFrow = np.arange(0.5, self.imshape[0] + 0.5)
		PRFcol = np.arange(0.5, self.imshape[1] + 0.5)
		
		# Center around 0 and convert to PSF subpixel resolution:
		PRFrow = (PRFrow - np.size(PRFrow) / 2) * self.superres
		PRFcol = (PRFcol - np.size(PRFcol) / 2) * self.superres
		
		# Interpolate highresImage:
		highresImageInterp = RectBivariateSpline(PRFcol, PRFrow, highresImage)
		
		# Integrate the interpolation to pixels:
		out = np.zeros(self.imshape)
		for row in range(self.imshape[0]):
			for col in range(self.imshape[1]):
				row_cen = - starpos[0]
				col_cen = - starpos[1]
				out[row,col] = highresImageInterp.integral(
					col_cen-0.5, col_cen+0.5, row_cen-0.5, row_cen+0.5)
		
		# Normalise the PSF:
		out /= np.nansum(out)
		
		return out

	def makeSmearKernel(self, starpos, integrationTime, angle, speed):
		"""
		Make a smear kernel that describes the large-scale movement of a star
		on the CCD. Do this by making a high resolution line that approximates
		the movement of a star with a certain constant angle and speed across 
		the CCD during the integrationTime.
		
		Input
		-----
		starpos (array, float)
			Row and column position in pixels of the star.
		integrationTime (float)
			CCD integration time.
		angle (float)
			Angle in radians of star CCD movement.
		speed (float)
			Speed of star CCD movement.
		
		Output
		------
		(array, float)
			Interpolated pixel positions of a star.
		"""
		imshapeHR = self.superres*self.imshape
		out = np.zeros(imshapeHR)
		r0 = starpos[0]
		c0 = starpos[1]
		r1 = r0 + np.int(self.superres*speed*integrationTime*np.sin(angle))
		c1 = c0 + np.int(self.superres*speed*integrationTime*np.cos(angle))
		# TODO: change implementation to subsubpixel line definition
		rr, cc = line(r0, c0, r1, c1)
		out[rr, cc] = 1
		return out

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
		"""
		return convolve2d(PSFunconvolved, kernel, 'same', 'fill', 0)

	def highresPSF(self, fwhm):
		"""
		Make a subpixel resolution PSF.
		
		Input
		-----
		fwhm (float)
			Full width at half maximum of Gaussian PSF in pixels.
		
		Output
		------
		(2D array, float)
			Subpixel-resolution PSF.
		"""
		# Make subpixel resolution grid with superres oversampling:
		PSFshapeHR = self.superres*self.imshape
		X, Y = np.meshgrid(np.arange(0,PSFshapeHR[1]), 
						np.arange(0,PSFshapeHR[0]))
		
		# Get coordinates of center of grid:
		x_0 = PSFshapeHR[1]/2 + 0.5
		y_0 = PSFshapeHR[0]/2 + 0.5
		
		# Evaluate PSF on grid with superres increased width:
		return self.integratedGaussian(X, Y, 1, x_0, y_0, 
					sigma = self.superres*fwhm/(2*np.sqrt(2*np.log(2))))

	def integratedGaussian(self, x, y, flux, x_0, y_0, sigma=1):
		'''
		Evaluate a 2D symmetrical Gaussian integrated in pixels.
		
		Parameters:
			x (numpy array) : x coordinates at which to evaluate the PSF.
			y (numpy array) : y coordinates at which to evaluate the PSF.
			flux (float) : Integrated value.
			x_0 (float) : Centroid position.
			y_0 (float) : Centroid position.
			sigma (float, optional) : Standard deviation of Gaussian. Default=1.
		
		Returns:
			numpy array : 2D Gaussian integrated pixel values at (x,y).
		
		Example:
		
		>>> import numpy as np
		>>> X, Y = np.meshgrid(np.arange(-1,2), np.arange(-1,2))
		>>> integratedGaussian(X, Y, 10, 0, 0)
		array([[ 0.58433556,  0.92564571,  0.58433556],
			[ 0.92564571,  1.46631496,  0.92564571],
			[ 0.58433556,  0.92564571,  0.58433556]])
		
		Code Author:
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
	
	bkg = np.zeros([1500,2000],dtype=float)
	dpsf = PSF(imshape=bkg.shape, superres=2)
	img = dpsf.evaluate(starpos=[100,150], integrationTime=30, angle=np.pi/4, 
			speed=40, fwhm=2)
	
	fig, ax = plt.subplots(1)
	ax.imshow(img)