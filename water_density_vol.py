#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog = 'water_density_vol', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
************************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/water_density_vol
DOI: 
************************************************

[ DESCRIPTION ]
 
This script calculates the water density in a specified volume of the simulation box.


[ REQUIREMENTS ]

The following python modules are needed :
 - MDAnalysis
 - matplotlib
 - numpy

 
[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: structure file [.gro] (required)
-x			: trajectory file [.xtc]
-o			: name of output folder
-b			: beginning time (ns) (the bilayer must exist by then!)
-e			: ending time (ns)	
-t 		[1]	: process every t-frames
-w		['W']	: resname of water molecules
 
Density profile options
-----------------------------------------------------
--sz		[200] 	: number of slizes along z
--sx		[200] 	: number of slizes along x and y

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#data options
parser.add_argument('-f', nargs=1, dest='grofilename', default=['no'], help=argparse.SUPPRESS, required=True)
parser.add_argument('-x', nargs=1, dest='xtcfilename', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-o', nargs=1, dest='output_folder', default=['no'], help=argparse.SUPPRESS)
parser.add_argument('-b', nargs=1, dest='t_start', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-e', nargs=1, dest='t_end', default=[-1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-t', nargs=1, dest='frames_dt', default=[1], type=int, help=argparse.SUPPRESS)
parser.add_argument('-w', nargs=1, dest='w_resname', default=['W'], help=argparse.SUPPRESS)

#density profile options
parser.add_argument('--sx', nargs=1, dest='sx', default=[200], type=int, help=argparse.SUPPRESS)
parser.add_argument('--sz', nargs=1, dest='sz', default=[200], type=int, help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
#data options
args.grofilename = args.grofilename[0]
args.xtcfilename = args.xtcfilename[0]
args.output_folder = args.output_folder[0]
args.t_start = args.t_start[0]
args.t_end = args.t_end[0]
args.frames_dt = args.frames_dt[0]
args.w_resname = args.w_resname[0]
#density profile options
args.sx = args.sx[0]
args.sz = args.sz[0]

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================
#generic science modules
try:
	import math
except:
	print "Error: you need to install the maths module."
	sys.exit(1)
try:
	import numpy as np
except:
	print "Error: you need to install the numpy module."
	sys.exit(1)
try:
	import scipy as sp
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)
try:
	import matplotlib as mpl
	mpl.use('Agg')
	import matplotlib.colors as mcolors
	mcolorconv = mcolors.ColorConverter()
	import matplotlib.cm as cm				#colours library
	import matplotlib.ticker
	from matplotlib.ticker import MaxNLocator
	from matplotlib.font_manager import FontProperties
	fontP=FontProperties()
except:
	print "Error: you need to install the matplotlib module."
	sys.exit(1)
try:
	import pylab as plt
except:
	print "Error: you need to install the pylab module."
	sys.exit(1)

#MDAnalysis module
try:
	import MDAnalysis
	from MDAnalysis import *
	import MDAnalysis.analysis
	import MDAnalysis.analysis.density
	import MDAnalysis.analysis.leaflet
	import MDAnalysis.analysis.distances
	#set MDAnalysis to use periodic boundary conditions
	MDAnalysis.core.flags['use_periodic_selections'] = True
	MDAnalysis.core.flags['use_KDTree_routines'] = False
except:
	print "Error: you need to install the MDAnalysis module first. See http://mdanalysis.googlecode.com"
	sys.exit(1)

#=========================================================================================
# sanity check
#=========================================================================================
if not os.path.isfile(args.grofilename):
	print "Error: file " + str(args.grofilename) + " not found."
	sys.exit(1)
if args.t_end != -1 and args.t_end < args.t_start:
	print "Error: the starting time (" + str(args.t_start) + "ns) for analysis is later than the ending time (" + str(args.t_end) + "ns)."
	sys.exit(1)

if args.xtcfilename == "no":
	if '-t' in sys.argv:
		print "Error: -t option specified but no xtc file specified."
		sys.exit(1)
	elif '-b' in sys.argv:
		print "Error: -b option specified but no xtc file specified."
		sys.exit(1)
	elif '-e' in sys.argv:
		print "Error: -e option specified but no xtc file specified."
		sys.exit(1)
elif not os.path.isfile(args.xtcfilename):
	print "Error: file " + str(args.xtcfilename) + " not found."
	sys.exit(1)

#=========================================================================================
# create folders and log file
#=========================================================================================
if args.output_folder == "no":
	if args.xtcfilename == "no":
		args.output_folder = "water_density_vol_" + args.grofilename[:-4]
	else:
		args.output_folder = "water_density_vol_" + args.xtcfilename[:-4]
if os.path.isdir(args.output_folder):
	print "Error: folder " + str(args.output_folder) + " already exists, choose a different output name via -o."
	sys.exit(1)
else:
	os.mkdir(args.output_folder)
	filename_log = os.getcwd() + '/' + str(args.output_folder) + '/water_density_vol.log'
	output_log = open(filename_log, 'w')		
	output_log.write("[water_density_vol v" + str(version_nb) + "]\n")
	output_log.write("\nThis folder and its content were created using the following command:\n\n")
	tmp_log="python water_density_vol.py"
	for c in sys.argv[1:]:
		tmp_log += " " + c
	output_log.write(tmp_log + "\n")
	output_log.close()
	
##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

global delta_x
global delta_y
global delta_z	
global coords_x
global coords_z
global z_middle
global upper_avg
global lower_avg
global w_density
global w_density_1D
global w_density_2D


w_density = np.zeros((args.sx, args.sx, args.sz))
w_density_1D = np.zeros(args.sz)
w_density_2D = np.zeros((args.sx,args.sz))

#=========================================================================================
# data loading
#=========================================================================================

def load_MDA_universe():
	
	global U
	global all_atoms
	global nb_atoms
	global nb_frames_xtc
	global frames_to_process
	global nb_frames_to_process
	global f_start
	global f_end
	global residues_list
	global water_nb
	global water_pres
	global water_sele
	global leaflet_sele
	global upper_sele
	global lower_sele
	f_start = 0
		
	#load universe
	#-------------
	if args.xtcfilename == "no":
		print "\nLoading file..."
		U = Universe(args.grofilename)
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = 1
		frames_to_process = [0]
		nb_frames_to_process = 1
	else:
		print "\nLoading trajectory..."
		U = Universe(args.grofilename, args.xtcfilename)
		U_timestep = U.trajectory.dt
		all_atoms = U.selectAtoms("all")
		nb_atoms = all_atoms.numberOfAtoms()
		nb_frames_xtc = U.trajectory.numframes

		U.trajectory.rewind()
		#sanity check
		if U.trajectory[nb_frames_xtc-1].time/float(1000) < args.t_start:
			print "Error: the trajectory duration (" + str(U.trajectory.time/float(1000)) + "ns) is shorted than the starting stime specified (" + str(args.t_start) + "ns)."
			sys.exit(1)
		if U.trajectory.numframes < args.frames_dt:
			print "Warning: the trajectory contains fewer frames (" + str(nb_frames_xtc) + ") than the frame step specified (" + str(args.frames_dt) + ")."

		#create list of index of frames to process
		if args.t_end != -1:
			f_end = int((args.t_end*1000 - U.trajectory[0].time) / float(U_timestep))
			if f_end < 0:
				print "Error: the starting time specified is before the beginning of the xtc."
				sys.exit(1)
		else:
			f_end = nb_frames_xtc - 1
		if args.t_start != -1:
			f_start = int((args.t_start*1000 - U.trajectory[0].time) / float(U_timestep))
			if f_start > f_end:
				print "Error: the starting time specified is after the end of the xtc."
				sys.exit(1)
		if (f_end - f_start)%args.frames_dt == 0:
			tmp_offset = 0
		else:
			tmp_offset = 1
		frames_to_process = map(lambda f:f_start + args.frames_dt*f, range(0,(f_end - f_start)//args.frames_dt+tmp_offset))
		nb_frames_to_process = len(frames_to_process)

	#identify leaflets the lazy way (in any case we need to assume to a plane bilayer as we plot against z so no need to be fancy)
	#------------------------------
	leaflet_sele = U.selectAtoms("name PO4")
	tmp_lipids_avg_z = leaflet_sele.centerOfGeometry()[2]
	upper_sele = leaflet_sele.selectAtoms("prop z > " + str(tmp_lipids_avg_z))
	lower_sele = leaflet_sele.selectAtoms("prop z < " + str(tmp_lipids_avg_z))

	#create water particles selection
	#--------------------------------
	water_pres = False
	water_sele = U.selectAtoms("resname " + args.w_resname)
	water_nb = water_sele.numberOfAtoms()
	if water_nb == 0:
		print "Error: particles with residue name '" + str(args.w_resname) + "'found, double check the -w option."
		sys.exit(1)
	else:
		water_pres = True

	return

#=========================================================================================
# core functions
#=========================================================================================

def calculate_density(f_index, box_dim):
	
	global delta_x
	global delta_y
	global delta_z	
	global upper
	global lower
	global w_density

	#define bins
	tmp_bins_x = np.linspace(0,box_dim[0],args.sx+1)
	tmp_bins_y = np.linspace(0,box_dim[1],args.sx+1)
	tmp_bins_z = np.linspace(0,box_dim[2],args.sz+1)
	delta_x = tmp_bins_x[1]-tmp_bins_x[0]
	delta_y = tmp_bins_y[1]-tmp_bins_y[0]	
	delta_z = tmp_bins_z[1]-tmp_bins_z[0]

	#store leaflets z position
	upper[f_index] = upper_sele.centerOfGeometry()[2]
	lower[f_index] = lower_sele.centerOfGeometry()[2]

	#get w coordinates
	w_coord = water_sele.coordinates()

	#transform coordinates into index of bins (the fanciness is to handle coords smaller/bigger than zero and U.dim)
	w_coord[:,0] = np.minimum(np.floor(w_coord[:,0]/float(delta_x)), np.remainder(np.floor(w_coord[:,0]/float(delta_x)),args.sx))
	w_coord[:,1] = np.minimum(np.floor(w_coord[:,1]/float(delta_y)), np.remainder(np.floor(w_coord[:,1]/float(delta_y)),args.sx))
	w_coord[:,2] = np.minimum(np.floor(w_coord[:,2]/float(delta_z)), np.remainder(np.floor(w_coord[:,2]/float(delta_z)),args.sz))
	w_coord = w_coord.astype(int)
	
	#bin water density
	w_density += np.histogramdd(w_coord,(range(0,args.sx+1), range(0,args.sx+1), range(0,args.sz+1)))[0]

	return
def calculate_stats():
	global coords_x
	global coords_z
	global z_middle
	global upper_avg
	global lower_avg
	global w_density
	global w_density_1D
	global w_density_2D
	
	#calculate coords and leaflets positions
	#---------------------------------------
	upper_avg = np.average(upper)
	lower_avg = np.average(lower)
	z_middle = (upper_avg+lower_avg) / float(2)
	upper_avg -= z_middle
	lower_avg -= z_middle	
	tmp_coords_x = np.linspace(0,U.dimensions[0],args.sx+1)
	tmp_coords_z = np.linspace(0,U.dimensions[2],args.sz+1)
	coords_x = tmp_coords_x[0:args.sx] + (tmp_coords_x[1]-tmp_coords_x[0])/float(2) - U.dimensions[0]/float(2)
	coords_z = tmp_coords_z[0:args.sz] + (tmp_coords_z[1]-tmp_coords_z[0])/float(2) - z_middle

	#calculate average volumic water density
	#---------------------------------------	
	#3D
	w_density = w_density / float(nb_frames_to_process)
	for nz in range(0, args.sz):
		#1D
		w_density_1D[nz] = np.average(w_density[:,:,nz])
		#2D
		for nx in range(0, args.sx):
			w_density_2D[nx,nz] = np.average(w_density[nx,:,nz])

	return

#=========================================================================================
# outputs
#=========================================================================================

#xvg file
def write_xvg_w_density():
	
	#open files
	filename_xvg = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_w_density_1D.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [water density profile - written by water_density_vol v" + str(version_nb) + "]\n")
	output_xvg.write("#  -> nb of slices x and y: " + str(args.sx) + "\n")
	output_xvg.write("#  -> nb of slices z: " + str(args.sz) + "\n")
	output_xvg.write("#  -> slices volume: " + str(round(delta_x*delta_y*delta_z,2)) + " (Angstrom3)\n")
	output_xvg.write("# nb of frames which contributed to this profile:\n")
	output_xvg.write("# -> weight = " + str(nb_frames_to_process) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Water density profile along z\"\n")
	output_xvg.write("@ xaxis label \"z distance to bilayer center (Angstrom)\"\n")
	output_xvg.write("@ yaxis label \"water particles relative frequency\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 1\n")
	output_xvg.write("@ s0 legend \"total\"\n")
	
	#data
	for nz in range(0,args.sz):
		results = str(coords_z[nz]) + "	" + "{:.6e}".format(w_density_1D[nz])
		output_xvg.write(results + "\n")	
	output_xvg.close()

		
	return

#graph
def graph_w_density():
			
	#1D profile
	#----------
	#filenames
	filename_svg = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_w_density_1D.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Water density profile along z")

	#plot data
	ax = fig.add_subplot(111)
	plt.plot(coords_z, w_density_1D, color = 'k', linewidth = 2)
	plt.vlines(lower_avg, min(w_density_1D), max(w_density_1D), linestyles = 'dashed')
	plt.vlines(upper_avg, min(w_density_1D), max(w_density_1D), linestyles = 'dashed')
	plt.vlines(0, min(w_density_1D), max(w_density_1D), linestyles = 'dashdot')
	plt.hlines(0, min(coords_z), max(coords_z))
	plt.hlines(0, min(coords_z), max(coords_z))
	plt.xlabel('z distance to bilayer center ($\AA$)')
	plt.ylabel('water partiles relative frequency')
	
	#save figure
	#ax.set_xlim(min(coords_z), max(coords_z))
	ax.set_xlim(-60, 60)
	#ax.set_ylim(min_density_charges, max_density_charges)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=9))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
	ax.xaxis.labelpad = 10
	ax.yaxis.labelpad = 10
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.15, right = 0.85)
	fig.savefig(filename_svg)
	plt.close()


	#2D profile
	#----------
	#filenames
	filename_svg = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_w_density_2D.svg'

	#create figure
	fig = plt.figure(figsize=(8, 6.2))
	fig.suptitle("Water density profile slice")

	#rotate data so that the x is horizontal and z vertical after imshow plotting
	w_density_2D_oriented = np.zeros((args.sz,args.sx))
	for nx in range(0, args.sx):
		for nz in range(0, args.sz):
			w_density_2D_oriented[nz,nx] = w_density_2D[nx,args.sz-1-nz]

	#plot data
	ax = fig.add_subplot(111)
	#im = plt.imshow(w_density_2D_oriented, extent = [min(coords_x),max(coords_x),min(coords_z),max(coords_z)], cmap = matplotlib.cm.jet_r, vmin = -0.08, vmax = 0.04)
	im = plt.imshow(w_density_2D_oriented[19:140,:], extent = [min(coords_x),max(coords_x),-60,60], cmap = matplotlib.cm.jet_r, vmin = 0, vmax = 0.01)
	plt.vlines(lower_avg, min(w_density_2D[:,0]), max(w_density_2D[:,0]), linestyles = 'dashed')
	plt.vlines(upper_avg, min(w_density_2D[:,0]), max(w_density_2D[:,0]), linestyles = 'dashed')
	plt.vlines(0, min(w_density_2D[:,0]), max(w_density_2D[:,0]), linestyles = 'dashdot')
	plt.xlabel('x ($\AA$)')
	plt.ylabel('z distance to bilayer center ($\AA$)')
	
	#color bar
	cax = fig.add_axes([0.85, 0.26, 0.025, 0.48])
	cbar = fig.colorbar(im, orientation='vertical', cax=cax)
	cbar.ax.tick_params(axis='y', direction='out')
	cbar.set_label(r'water particles relative frequency')
		
	#save figure
	ax.set_xlim(min(coords_x), max(coords_x))
	#ax.set_ylim(min(coords_z), max(coords_z))
	ax.set_ylim(-60,60)
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_major_locator(MaxNLocator(nbins=9))
	ax.yaxis.set_major_locator(MaxNLocator(nbins=7))
	ax.xaxis.labelpad = 10
	ax.yaxis.labelpad = 10
	ax.tick_params(axis='x', direction='out')
	ax.tick_params(axis='y', direction='out')
	plt.setp(ax.xaxis.get_majorticklabels(), fontsize = "small")
	plt.setp(ax.yaxis.get_majorticklabels(), fontsize = "small")
	plt.subplots_adjust(top = 0.9, bottom = 0.15, left = 0.1, right = 0.8)
	fig.savefig(filename_svg)
	plt.close()

	return

#Opendx file
def write_dx_w_density():

	filename_dx = os.getcwd() + '/' + args.output_folder + '/' + str(args.xtcfilename[:-4]) + '_w_density_3D.dx'
	D = MDAnalysis.analysis.density.Density(w_density,[np.arange(0,args.sx),np.arange(0,args.sx),np.arange(0,args.sz)])
	D.export(filename_dx)

	return

##########################################################################################
# ALGORITHM
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================
load_MDA_universe()
upper = np.zeros(nb_frames_to_process)
lower = np.zeros(nb_frames_to_process)

#=========================================================================================
# generate data
#=========================================================================================
print "\nCalculating water density..."
#case: structure only
#--------------------
if args.xtcfilename=="no":
	calculate_potential(0,U.trajectory.ts.dimensions)
#case: browse xtc frames
#-----------------------
else:
	for f_index in range(0,nb_frames_to_process):
		ts = U.trajectory[frames_to_process[f_index]]
		progress = '\r -processing frame ' + str(f_index+1) + '/' + str(nb_frames_to_process) + ' (every ' + str(args.frames_dt) + ' frame(s) from frame ' + str(f_start) + ' to frame ' + str(f_end) + ' out of ' + str(nb_frames_xtc) + ')      '  
		sys.stdout.flush()
		sys.stdout.write(progress)							
		calculate_density(f_index, U.trajectory.ts.dimensions)
	print ""

#=========================================================================================
# process data
#=========================================================================================
calculate_stats()

#=========================================================================================
# produce outputs
#=========================================================================================
print "\nWriting outputs..."
graph_w_density()
write_xvg_w_density()
write_dx_w_density()
	
#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check output in ./" + args.output_folder + "/"
print ""
sys.exit(0)
