#import matplotlib.pyplot as plt
#import math as mt
#import matplotlib
from PIL import Image
import os
import numpy as np
from scipy import misc
import scipy

##############################################################################
##############################################################################

infiledir = "../infiles"
outfiledir = "../outfiles"

# Load jpg file:
# TODO: what is this used for?
jpgfile = Image.open(os.path.join("../infiles","Southern_Cross.jpg"))
print(jpgfile.bits, jpgfile.size, jpgfile.format)
jpgfile.size = np.array([200,200])

# Define number of output frames:
biasnr = 1
flatsnr = 1
sciencenr = 1
#biasnr = 10
#flatsnr = 3
#sciencenr = 30

biasmeanr = 3.75394283333 ; biasstdevr = 0.532235493781
biasmeang = 2.34252829167 ; biasstdevg = 0.506928983646
biasmeanb = 4.36334695833 ; biasstdevb = 0.993081162917

# Create bias images. #################################################################################

for i in range(biasnr):
    bias_RGB = np.zeros((jpgfile.size[0], jpgfile.size[1], 3), "uint8") 
    biasr = np.random.normal(biasmeanr, biasstdevr, jpgfile.size[0]*jpgfile.size[1]).reshape((jpgfile.size[0], jpgfile.size[1]))
    biasg = np.random.normal(biasmeang, biasstdevg, jpgfile.size[0]*jpgfile.size[1]).reshape((jpgfile.size[0], jpgfile.size[1]))
    biasb = np.random.normal(biasmeanb, biasstdevb, jpgfile.size[0]*jpgfile.size[1]).reshape((jpgfile.size[0], jpgfile.size[1]))

    bias_RGB[:,:,0] = biasr
    bias_RGB[:,:,1] = biasg
    bias_RGB[:,:,2] = biasb

    imout = 'bias_'+str(i)+'.jpg'
    
    scipy.misc.imsave(os.path.join(outfiledir,'bias_%02d.jpg' % i), bias_RGB)
    
# Create flats images. #################################################################################

master_flat = np.zeros((jpgfile.size[0], jpgfile.size[1], 3), "uint8") 
master_flatr = np.zeros(jpgfile.size[0]*jpgfile.size[1]).reshape(jpgfile.size[0], jpgfile.size[1])
master_flatg = np.zeros(jpgfile.size[0]*jpgfile.size[1]).reshape(jpgfile.size[0], jpgfile.size[1]) 
master_flatb = np.zeros(jpgfile.size[0]*jpgfile.size[1]).reshape(jpgfile.size[0], jpgfile.size[1])

flatmeanr = 133.7070775 ; flatstdr = 20.2505762111
flatmeang = 133.0223080 ; flatstdg = 19.9602947291
flatmeanb = 132.3117177 ; flatstdb = 20.2094106756

for i in range(flatsnr):
    flats_RGB = np.zeros((jpgfile.size[0], jpgfile.size[1], 3), "uint8") 
    flatsr = np.random.normal(flatmeanr, flatstdr, jpgfile.size[0]*jpgfile.size[1]).reshape(jpgfile.size[0], jpgfile.size[1])
    flatsg = np.random.normal(flatmeang, flatstdg, jpgfile.size[0]*jpgfile.size[1]).reshape(jpgfile.size[0], jpgfile.size[1])
    flatsb = np.random.normal(flatmeanb, flatstdb, jpgfile.size[0]*jpgfile.size[1]).reshape(jpgfile.size[0], jpgfile.size[1])
    
    for l1 in range(jpgfile.size[0]):
        for l2 in range(jpgfile.size[1]):
                
            #central illumination/vignetting.
            flatsr[l1,l2] = flatsr[l1,l2] + 35.*np.exp( -( (float(l1)-jpgfile.size[0]/2.)**2/((2.*jpgfile.size[0]/5.)**2) + (float(l2)-jpgfile.size[1]/2.)**2/((2.*jpgfile.size[1]/5.)**2) ) )
            flatsg[l1,l2] = flatsr[l1,l2] + 35.*np.exp( -( (float(l1)-jpgfile.size[0]/2.)**2/((2.*jpgfile.size[0]/5.)**2) + (float(l2)-jpgfile.size[1]/2.)**2/((2.*jpgfile.size[1]/5.)**2) ) )
            flatsb[l1,l2] = flatsr[l1,l2] + 35.*np.exp( -( (float(l1)-jpgfile.size[0]/2.)**2/((2.*jpgfile.size[0]/5.)**2) + (float(l2)-jpgfile.size[1]/2.)**2/((2.*jpgfile.size[1]/5.)**2) ) )
       
    flats_RGB[:,:,0] = flatsr
    flats_RGB[:,:,1] = flatsg
    flats_RGB[:,:,2] = flatsb

    misc.imsave(os.path.join(outfiledir, 'flats_%02d.jpg' % i), flats_RGB)

    master_flatr = master_flatr + flatsr
    master_flatg = master_flatg + flatsg
    master_flatb = master_flatb + flatsb

master_flatr = master_flatr/float(flatsnr)
master_flatg = master_flatg/float(flatsnr)
master_flatb = master_flatb/float(flatsnr)

master_flat[:,:,0] = master_flatr
master_flat[:,:,1] = master_flatg
master_flat[:,:,2] = master_flatb

misc.imsave(os.path.join(outfiledir,'master_flat.jpg'), master_flat)
misc.imsave(os.path.join(outfiledir,'master_flat_norm.jpg'), master_flat/np.mean(master_flat))
  
# Create science frames over master flat. ################################################################

sciencer = np.zeros(jpgfile.size[0]*jpgfile.size[1]).reshape((jpgfile.size[0], jpgfile.size[1]))
scienceg = np.zeros(jpgfile.size[0]*jpgfile.size[1]).reshape((jpgfile.size[0], jpgfile.size[1]))
scienceb = np.zeros(jpgfile.size[0]*jpgfile.size[1]).reshape((jpgfile.size[0], jpgfile.size[1]))

for i in range(sciencenr):
    science_RGB = np.zeros((jpgfile.size[0], jpgfile.size[1], 3), "uint8") 

    randx = np.random.uniform(0.2*jpgfile.size[0], 0.8*jpgfile.size[0])
    randy = np.random.uniform(0.2*jpgfile.size[1], 0.8*jpgfile.size[1])
    sigmax_r = 13.0 + np.random.uniform(0.2,0.3)
    sigmay_r = 13.0 + np.random.uniform(0.1,0.2)
    sigmax_g = 9.5 + np.random.uniform(0.2,0.3)
    sigmay_g = 9.5 + np.random.uniform(0.1,0.2)
    sigmax_b = 9.0 + np.random.uniform(0.2,0.3)
    sigmay_b = 9.0 + np.random.uniform(0.1,0.2)

    for l1 in range(jpgfile.size[0]):
        for l2 in range(jpgfile.size[1]):

            biaslevel = np.random.normal(biasmeanr*4., biasstdevr*4.)

            sciencer[l1,l2] = (master_flat[l1,l2,0]/np.mean(master_flat))*biaslevel + 10. + 190.*(np.exp(- ((float(l1) - randx)**2/(2.*sigmax_r**2) + (float(l2) - randy)**2/(2.*sigmax_r**2) ))) 
            scienceg[l1,l2] = (master_flat[l1,l2,1]/np.mean(master_flat))*biaslevel + 5. + 190.*(np.exp(- ((float(l1) - randx)**2/(2.*sigmax_g**2) + (float(l2) - randy)**2/(2.*sigmax_g**2) )))
            scienceb[l1,l2] = (master_flat[l1,l2,2]/np.mean(master_flat))*biaslevel + 8. + 190.*(np.exp(- ((float(l1) - randx)**2/(2.*sigmax_b**2) + (float(l2) - randy)**2/(2.*sigmax_b**2) )))

    science_RGB[:,:,0] = sciencer
    science_RGB[:,:,1] = scienceg
    science_RGB[:,:,2] = scienceb

    misc.imsave(os.path.join(outfiledir,'science_%02d.jpg' % i), science_RGB)
