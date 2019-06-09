# @String jsonpath

from __future__ import with_statement
import os
import os.path as op
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
	List_of_fits = []

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
				IJ.run("Results Manager", "input_file=["+fitresults+"] input=File input_file=["+fitresults+"] "+
				"table=ImageJ image=[Localisations (width=precision)] results_format=None results_file=[] "+
				"table_distance_unit=[pixel (px)] table_intensity_unit=photon table_angle_unit=[unknown (na)] "+
				"table_precision=4 equalised image_scale="+str(sr_scale))
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