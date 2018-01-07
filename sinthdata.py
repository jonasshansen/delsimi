#import matplotlib.pyplot as plt
#import math as mt
#import matplotlib
#from scipy import misc
from PIL import Image
import os
import numpy as np
import scipy

##############################################################################
##############################################################################

class delsimi(object):
	def __init__(self):
		"""
		Simulate stellar images from Delphini-1.
		
		Future Extensions
		-----------------
		This is a list of possible future extensions to this code. The list is
		largely inspired by the SPyFFI manual which is linked to here:
			
			https://github.com/TESScience/SPyFFI
			
		- Smeared PSF with length inferred from satellite speed and exposure time
		- True pixel sizes / frame size
		- Catalog stars
		- Additional stars
		- Pixel sensitivity variations?
		- Background
		- Photon noise
		- Readout smear
		- Saturated pixels
		- Cosmic rays
		- Position variable PSF
		
		Code Authors
		------------
		Carolina von Essen, cessen@phys.au.dk
		Jonas Svenstrup Hansen, jonas.svenstrup@gmail.com
		"""
		
		self.infiledir = "../infiles"
		self.outfiledir = "../outfiles"
		
		# Load jpg file:
		# TODO: what is this used for? jpg encoding?
		self.jpgfile = Image.open(os.path.join("../infiles","Southern_Cross.jpg"))
#		print(jpgfile.bits, jpgfile.size, jpgfile.format)
		self.jpgfile.size = np.array([200,200])
		
		# Define number of output frames of each type:
		self.biasnr = 1
		self.flatsnr = 1
		self.sciencenr = 1
		#self.biasnr = 10
		#self.flatsnr = 3
		#self.sciencenr = 30


	def makebias(self):
		"""
		Make bias frames.
		"""
		
		biasmeanr = 3.75394283333 ; biasstdevr = 0.532235493781
		biasmeang = 2.34252829167 ; biasstdevg = 0.506928983646
		biasmeanb = 4.36334695833 ; biasstdevb = 0.993081162917
		
		self.biasmeanr = biasmeanr
		self.biasstdevr = biasstdevr
		
		for i in range(self.biasnr):
			bias_RGB = np.zeros((self.jpgfile.size[0], self.jpgfile.size[1], 3), "uint8") 
			biasr = np.random.normal(biasmeanr, biasstdevr, self.jpgfile.size[0]*self.jpgfile.size[1]).reshape((self.jpgfile.size[0], self.jpgfile.size[1]))
			biasg = np.random.normal(biasmeang, biasstdevg, self.jpgfile.size[0]*self.jpgfile.size[1]).reshape((self.jpgfile.size[0], self.jpgfile.size[1]))
			biasb = np.random.normal(biasmeanb, biasstdevb, self.jpgfile.size[0]*self.jpgfile.size[1]).reshape((self.jpgfile.size[0], self.jpgfile.size[1]))
		
			bias_RGB[:,:,0] = biasr
			bias_RGB[:,:,1] = biasg
			bias_RGB[:,:,2] = biasb
		
#			imout = 'bias_'+str(i)+'.jpg'
		
			scipy.misc.imsave(os.path.join(self.outfiledir,'bias_%02d.jpg' % i), bias_RGB)


	def makeflat(self):
		"""
		Make flat frames.
		"""
		
		master_flat = np.zeros((self.jpgfile.size[0],self.jpgfile.size[1], 3), "uint8") 
		master_flatr = np.zeros(self.jpgfile.size[0]*self.jpgfile.size[1]).reshape(self.jpgfile.size[0], self.jpgfile.size[1])
		master_flatg = np.zeros(self.jpgfile.size[0]*self.jpgfile.size[1]).reshape(self.jpgfile.size[0], self.jpgfile.size[1]) 
		master_flatb = np.zeros(self.jpgfile.size[0]*self.jpgfile.size[1]).reshape(self.jpgfile.size[0], self.jpgfile.size[1])
		
		flatmeanr = 133.7070775 ; flatstdr = 20.2505762111
		flatmeang = 133.0223080 ; flatstdg = 19.9602947291
		flatmeanb = 132.3117177 ; flatstdb = 20.2094106756
		
		for i in range(self.flatsnr):
			flats_RGB = np.zeros((self.jpgfile.size[0], self.jpgfile.size[1], 3), "uint8") 
			flatsr = np.random.normal(flatmeanr, flatstdr, self.jpgfile.size[0]*self.jpgfile.size[1]).reshape(self.jpgfile.size[0], self.jpgfile.size[1])
			flatsg = np.random.normal(flatmeang, flatstdg, self.jpgfile.size[0]*self.jpgfile.size[1]).reshape(self.jpgfile.size[0], self.jpgfile.size[1])
			flatsb = np.random.normal(flatmeanb, flatstdb, self.jpgfile.size[0]*self.jpgfile.size[1]).reshape(self.jpgfile.size[0], self.jpgfile.size[1])
		
			for l1 in range(self.jpgfile.size[0]):
				for l2 in range(self.jpgfile.size[1]):
		
					#central illumination/vignetting.
					flatsr[l1,l2] = flatsr[l1,l2] + 35.*np.exp( -( (float(l1)-self.jpgfile.size[0]/2.)**2/((2.*self.jpgfile.size[0]/5.)**2) + (float(l2)-self.jpgfile.size[1]/2.)**2/((2.*self.jpgfile.size[1]/5.)**2) ) )
					flatsg[l1,l2] = flatsr[l1,l2] + 35.*np.exp( -( (float(l1)-self.jpgfile.size[0]/2.)**2/((2.*self.jpgfile.size[0]/5.)**2) + (float(l2)-self.jpgfile.size[1]/2.)**2/((2.*self.jpgfile.size[1]/5.)**2) ) )
					flatsb[l1,l2] = flatsr[l1,l2] + 35.*np.exp( -( (float(l1)-self.jpgfile.size[0]/2.)**2/((2.*self.jpgfile.size[0]/5.)**2) + (float(l2)-self.jpgfile.size[1]/2.)**2/((2.*self.jpgfile.size[1]/5.)**2) ) )
		
			flats_RGB[:,:,0] = flatsr
			flats_RGB[:,:,1] = flatsg
			flats_RGB[:,:,2] = flatsb
		
			scipy.misc.imsave(os.path.join(self.outfiledir, 'flats_%02d.jpg' % i), flats_RGB)
		
			master_flatr = master_flatr + flatsr
			master_flatg = master_flatg + flatsg
			master_flatb = master_flatb + flatsb
		
		master_flatr = master_flatr/float(self.flatsnr)
		master_flatg = master_flatg/float(self.flatsnr)
		master_flatb = master_flatb/float(self.flatsnr)
		
		master_flat[:,:,0] = master_flatr
		master_flat[:,:,1] = master_flatg
		master_flat[:,:,2] = master_flatb
		
		scipy.misc.imsave(os.path.join(self.outfiledir,'master_flat.jpg'), master_flat)
		scipy.misc.imsave(os.path.join(self.outfiledir,'master_flat_norm.jpg'), master_flat/np.mean(master_flat))
		
		self.master_flat = master_flat


	def makescience(self):
		"""
		Make science frames using the master flat.
		"""
		
		sciencer = np.zeros(self.jpgfile.size[0]*self.jpgfile.size[1]).reshape((self.jpgfile.size[0], self.jpgfile.size[1]))
		scienceg = np.zeros(self.jpgfile.size[0]*self.jpgfile.size[1]).reshape((self.jpgfile.size[0], self.jpgfile.size[1]))
		scienceb = np.zeros(self.jpgfile.size[0]*self.jpgfile.size[1]).reshape((self.jpgfile.size[0], self.jpgfile.size[1]))
		
		for i in range(self.sciencenr):
			science_RGB = np.zeros((self.jpgfile.size[0], self.jpgfile.size[1], 3), "uint8") 
		
			randx = np.random.uniform(0.2*self.jpgfile.size[0], 0.8*self.jpgfile.size[0])
			randy = np.random.uniform(0.2*self.jpgfile.size[1], 0.8*self.jpgfile.size[1])
			sigmax_r = 13.0 + np.random.uniform(0.2,0.3)
			sigmay_r = 13.0 + np.random.uniform(0.1,0.2)
			sigmax_g = 9.5 + np.random.uniform(0.2,0.3)
			sigmay_g = 9.5 + np.random.uniform(0.1,0.2)
			sigmax_b = 9.0 + np.random.uniform(0.2,0.3)
			sigmay_b = 9.0 + np.random.uniform(0.1,0.2)
		
			for l1 in range(self.jpgfile.size[0]):
				for l2 in range(self.jpgfile.size[1]):
					
					biaslevel = np.random.normal(self.biasmeanr*4., self.biasstdevr*4.)
					
					sciencer[l1,l2] = (self.master_flat[l1,l2,0]/np.mean(self.master_flat))*biaslevel + 10. + 190.*(np.exp(- ((float(l1) - randx)**2/(2.*sigmax_r**2) + (float(l2) - randy)**2/(2.*sigmax_r**2) ))) 
					scienceg[l1,l2] = (self.master_flat[l1,l2,1]/np.mean(self.master_flat))*biaslevel + 5. + 190.*(np.exp(- ((float(l1) - randx)**2/(2.*sigmax_g**2) + (float(l2) - randy)**2/(2.*sigmax_g**2) )))
					scienceb[l1,l2] = (self.master_flat[l1,l2,2]/np.mean(self.master_flat))*biaslevel + 8. + 190.*(np.exp(- ((float(l1) - randx)**2/(2.*sigmax_b**2) + (float(l2) - randy)**2/(2.*sigmax_b**2) )))
		
			science_RGB[:,:,0] = sciencer
			science_RGB[:,:,1] = scienceg
			science_RGB[:,:,2] = scienceb
		
			scipy.misc.imsave(os.path.join(self.outfiledir,'science_%02d.jpg' % i), science_RGB)


if __name__ == '__main__':
	simtest = delsimi()
	print('Making bias...')
	simtest.makebias()
	print('...done!')
	print('Making flat')
	simtest.makeflat()
	print('...done!')
	print('Making science')
	simtest.makescience()
	print('...done!')