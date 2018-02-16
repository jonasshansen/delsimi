#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run function for the delsimi code for commandline calls.

@author: Jonas Svenstrup Hansen, jonas.svenstrup@gmail.com
"""

import argparse
from delsimi import delsimi

if __name__ == '__main__':
	# Parse command-line arguments:
	parser = argparse.ArgumentParser(description=
			'Run the delsimi code to simulate stellar images from Delphini-1.')
	parser.add_argument('-i', '--input_dir', 
			help='Input file directory.',
			default='../infiles')
	parser.add_argument('-o', '--output_dir', 
			help='Output file directory.',
			default='../outfiles')
	parser.add_argument('-w', '--overwrite', 
			help='True if to overwrite FITS images in output.',
			default=True)
	parser.add_argument('-c', '--coord_cen', 
			help='Coordinate center at the midtime of exposure.',
			default=[56.75,24.11670])
	parser.add_argument('-t', '--integration_time', 
			help='Integration (exposure) time in seconds.',
			default=1.)
	parser.add_argument('-v', '--angle_vel', 
			help='Velocity angle.',
			default=0.)
	parser.add_argument('-s', '--angle_sat', 
			help='Satellite angle.',
			default=0.)
	
	
	# Save to variables:
	args = parser.parse_args()
	
	input_dir = args.input_dir
	output_dir = args.output_dir
	overwrite = args.overwrite
	coord_cen = args.coord_cen
	integration_time = args.integration_time
	angle_vel = args.angle_vel
	angle_sat = args.angle_sat
	
	
	# Run delsimi.py with the specified parameters:
	delsimi(input_dir=input_dir,
			output_dir=output_dir,
			overwrite=overwrite,
			coord_cen=coord_cen,
			integration_time=integration_time,
			angle_vel=angle_vel,
			angle_sat=angle_sat)