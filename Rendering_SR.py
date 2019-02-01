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
	for roots, dirs, files in os.walk(results_dir):
		for f in files:
			if f == 'All_fits_with_header_py.txt':
				List_of_fits.append((roots, f))

	for r, f in List_of_fits:
		fitresults = op.join(r, f)
		roi_file = op.join(r, 'clusters_roi.txt')
		try:
			IJ.run("Results Manager", "coordinate=["+
			fitresults+"] input=File input_file=["+fitresults+"] results_table=Uncalibrated "+
			"image=[Localisations (width=precision)] weighted equalised image_precision=1.50 "+
			"image_scale="+str(sr_scale)+" image_window=0")
		except:
			continue
		wm.getWindow("Fit Results").close()
		sr = wm.getImage("Image (LSE) SuperRes")
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
		IJ.saveAs(sr, "Tiff", op.join(r, "Rendered_SR_labelled.tif"))
		rm.reset()
		IJ.saveAs(sr, "Tiff", op.join(r, "Rendered_SR.tif"))
		sr.close()

	IJ.run("Quit")