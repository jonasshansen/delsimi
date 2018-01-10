import numpy as np
#import scipy
from scipy.special import erf


class PSF():
	def __init__(self):
		"""
		Approximate the PSF in stellar images from Delphini-1.
		
		Use the following sources to smear and jitter:
		https://www.osapublishing.org/ao/fulltext.cfm?uri=ao-32-32-6503&id=40363
		https://github.com/TESScience/SPyFFI/blob/master/PSF.py
		https://github.com/TESScience/SPyFFI/blob/master/Jitter.py
		"""

	def makeJitter(self):
		pass
	
	def makeSmear(self):
		pass

	def makeFocus(self):
		pass

	def convolvePSF(self):
#		scipy.signal.convolve2d
		pass


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
