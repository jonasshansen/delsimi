#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities for the delsimi simulation code.

@author: jonas
"""

# TODO: create function to convert RGB to Johnson filters and back
def uvb2rgb(uvb):
	# Assume the following conversion applies to space photometry
	#https://www.sciencedirect.com/science/article/pii/S0273117715005694?via%3Dihub
	rgb = uvb
	return rgb


def bayer_scaling(img, flux):
	"""
	Apply Bayer filter scaling of a normalised image with the level of flux.
	
	Parameters:
	-----------
	img (numpy array):
		2D image with an even number of rows and columns.
	flux (list of floats):
		List with the three values for red, green and blue to scale image.
	"""
	# Red:
	img[1::2,1::2] *= flux[0]
	# Green:
	img[0::2,1::2] *= flux[1]
	img[1::2,0::2] *= flux[1]
	# Blue:
	img[0::2,0::2] *= flux[2]
	return img