import numpy as np
from scipy.special import erf
from scipy.signal import convolve2d
from scipy.interpolate import RectBivariateSpline
from scipy.interpolate import splprep
from skimage.draw import line_aa

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

	def evaluate(self, star, integrationTime, angle, speed, fwhm = 1, 
			jitter = False, focus = False):
		"""
		Evaluate a PSF that is smeared in one direction.
		
		Input
		-----
		star (array, float)
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
		# Define subpixel buffer:
		self.buffer = np.int(3*self.superres*fwhm)
		
		# Create smear kernel:
		smearKernel = self.makeSmearKernel(integrationTime, angle, speed, fwhm)
		self.kernelShape = smearKernel.shape
		
		# Get highres PSF:
		PSFhighres = self.highresPSF(fwhm)
		
		# TODO: convolve highres PSF with focus and jitter here
		
		# Convolve the PSF with the smear kernel:
		highresImage = self.convolvePSF(PSFhighres, smearKernel)
		
		# Normalise the PRF:
		highresImage /= np.nansum(highresImage) * self.superres**2
		
		# Define pixel centered index arrays for the interpolater:
		PRFrow = np.arange(0.5, self.kernelShape[0] + 0.5)
		PRFcol = np.arange(0.5, self.kernelShape[1] + 0.5)
		
		# Center around 0:
		PRFrow = PRFrow - PRFrow.size / 2
		PRFcol = PRFcol - PRFcol.size / 2
		
		# Move center to account for buffer zone:
		PRFrow += self.buffer
		PRFcol += self.buffer
		
		# Convert from subpixel to pixel resolution:
		PRFrow /= self.superres
		PRFcol /= self.superres
		
		# Interpolate highresImage:
		highresImageInterp = RectBivariateSpline(PRFrow, PRFcol, highresImage)
		
		# Integrate the interpolation in each pixel:
		# (integration of the spline outside the spline boundaries yield artefacts)
		img = np.zeros(self.imshape, dtype='float64')
		for row in range(self.imshape[0]):
			for col in range(self.imshape[1]):
				# Get star position in PSF(t=0)-based coordinates:
				row_cen = row - star[0]
				col_cen = col - star[1]
				# Determine whether to integrate current pixel value:
				withinBoundary = highresImageInterp(row_cen, col_cen) > 1e-9
				if withinBoundary:
					# Integrate intepolation in the current pixel:
					img[row,col] = highresImageInterp.integral(row_cen-0.5, row_cen+0.5, col_cen-0.5, col_cen+0.5)
		
		return img, smearKernel, PSFhighres, highresImage, highresImageInterp

	def makeSmearKernel(self, integrationTime, angle, speed, fwhm):
		"""
		Make a smear kernel that describes the large-scale movement of a star
		on the CCD. Do this by making a high resolution line that approximates
		the movement of a star with a certain constant angle and speed across 
		the CCD during the integrationTime.
		
		Input
		-----
		integrationTime (float)
			CCD integration time.
		angle (float)
			Angle in radians of star CCD movement.
		speed (float)
			Speed of star CCD movement.
		fwhm (float)
			Full width at half maximum of PSF. Used to determine buffer pixels
			around the line.
		
		Output
		------
		(array, float)
			Interpolated pixel positions of a star.
		"""
		buffer = self.buffer
		finalposRow = np.int(self.superres*speed*integrationTime*np.sin(angle))
		finalposCol = np.int(self.superres*speed*integrationTime*np.cos(angle))
		r0 = 0 + buffer
		c0 = 0 + buffer
		r1 = r0 + finalposRow
		c1 = c0 + finalposCol
		out = np.zeros([r1+buffer, c1+buffer])
		# TODO: change implementation to subsubpixel line definition
#		rr, cc = line(r0, c0, r1, c1)
#		out[rr, cc] = 1
		rr, cc, val = line_aa(r0, c0, r1, c1)
		out[rr, cc] = val
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
		X, Y = np.meshgrid(np.arange(0,self.kernelShape[1]), 
						np.arange(0,self.kernelShape[0]))
		
		# Define centroid position:
		# TODO: fix bug where the starting position changes when the integrationTime changes
		x_0 = self.kernelShape[1]/2 + 0.5
		y_0 = self.kernelShape[0]/2 + 0.5
		
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
	
	# Define background:
	bkg = np.zeros([30,40],dtype=float)
	
	# Make PSF class instance:
	dpsf = PSF(imshape=bkg.shape, superres=3)
	
	# Evaluate PSF with specified parameters:
	img, smearKernel, PSFhighres, highresImage, highresImageInterp = dpsf.evaluate(
			star=[5,10], integrationTime=10, angle=np.pi/3, speed=1, fwhm=1)
	
	# Plot:
	fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)
	
	ax1.imshow(img, origin='lower')
	ax1.set_xlabel('Pixel column')
	ax1.set_ylabel('Pixel row')
	ax1.set_title('Pixel-integrated image')
	
	ax2.imshow(smearKernel, origin='lower')
	ax2.set_title('Smear kernel')
	
	ax3.imshow(highresImage, origin='lower')
	ax3.set_title('High res. convolved PRF')
	
	ax4.imshow(PSFhighres, origin='lower')
	ax4.set_title('High resolution PRF')
	
	for ax in (ax2, ax3, ax4):
		ax.set_xlabel('Subpixel column')
		ax.set_ylabel('Subpixel row')
	
	fig.subplots_adjust(hspace=0.5)
