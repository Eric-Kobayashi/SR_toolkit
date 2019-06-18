# @String jsonpath

from __future__ import with_statement
import os
import os.path as op
import shutil as sh
from ij import IJ
from ij import WindowManager as wm
import ij.plugin.frame.RoiManager as RoiManager
import json
IJ.redirectErrorMessages()

if __name__ in ['__builtin__', '__main__']:
	with open(jsonpath, 'r') as f:
		d = json.load(f)	# Configuration dictionary
	results_dir = d['results_dir']
	sr_scale = d['sr_scale']

	# Fix the configuration file non-existence problem
	imageJ_path = os.getcwd()
	gdsc_xml = op.join(imageJ_path, 'gdsc.smlm.settings.xml')
	# make a fake configuration file
	IJ.run("Fit Configuration", "configuration_file={0} config_file={1} calibration=107 ".format(gdsc_xml.replace('\\', '/'), gdsc_xml.replace('\\', '/'))+\
	"gain=84.4 exposure_time=50 initial_stddev0=2 initial_stddev1=2 initial_angle=0 "+\
	"spot_filter_type=Single spot_filter=Mean smoothing=0.5 search_width=3 border=1 "+\
	"fitting_width=3 fit_solver=[Least Squares Estimator (LSE)] fit_function=Circular "+\
	"fail_limit=10 include_neighbours neighbour_height=0.3 residuals_threshold=1 "+\
	"duplicate_distance=0.5 shift_factor=2 signal_strength=0 min_photons=0 min_width_factor=0.5 "+\
	"width_factor=2 precision_threshold=40 fit_criteria=[Least-squared error] "+\
	"significant_digits=5 coord_delta=0.0001 lambda=10.0000 max_iterations=20")
	for to_fit in ['burst', 'mol', 'cluster']:
		List_of_fits = []
		fit_name = 'All_{}_header.txt'.format(to_fit)
		for roots, dirs, files in os.walk(results_dir):
			for f in files:
				if f == fit_name:
					List_of_fits.append((roots, f))
		
		for r, f in List_of_fits:
			IJ.run("Close All")
			fitresults = op.join(r, f)
			try:
				IJ.redirectErrorMessages()
				IJ.run("Results Manager", "coordinate=["+
				fitresults+"] input=File input_file=["+fitresults+"] results_table=Uncalibrated "+
				"image=[Localisations (width=precision)] weighted equalised image_precision=1.50 "+
				"image_scale="+str(sr_scale)+" image_window=0")
			except:
				continue
			
			try:
				IJ.redirectErrorMessages()
				wm.getWindow("Fit Results").close()
			except:
				pass

			try:
				IJ.redirectErrorMessages()
				sr = wm.getImage("Image (LSE) SuperRes")
			except:
				continue
				
			if to_fit == 'cluster':
				try:
					roi_file = op.join(r, 'clusters_roi.txt')
					rm = RoiManager.getInstance()
					if not rm:
						rm = RoiManager()
					rm.reset()
					i = 0
					with open(roi_file, 'r') as roi_f:
						for lines in roi_f:
							if i == 0: 
								i += 1
								continue
							else:
								eval("sr.setRoi({})".format(lines[:-1]))
								rm.addRoi(sr.getRoi())
								sr.killRoi()
								
					rm.runCommand("Save", op.join(r, "roiset.zip"))
					rm.runCommand(sr, "Show All with labels")
					IJ.saveAs(sr, "Tiff", op.join(r, "cluster_SR_labelled.srf.tif"))
					rm.reset()
				except:
					pass
			try:
				IJ.redirectErrorMessages()
				IJ.saveAs(sr, "Tiff", op.join(r, "{}_SR.srf.tif".format(to_fit)))
			except:
				continue
				
			IJ.run("Close All")

	IJ.run("Quit")