# @String jsonpath

from __future__ import with_statement
import os
from ij import IJ
from ij import WindowManager as wm
import ij.plugin.frame.RoiManager as RoiManager
import ij.plugin.ZProjector as ZProjector
import json
from ij.plugin import SubstackMaker

def GDSC_SMLM(directory, file_name, results_dir, dirnum, trim_track, gdsc_smlm_xml, pixel_size, frame_length, signal_strength, precision, min_photons, sr_scale):
	try:
		wm.getWindow("Fit Results").close()	# clean the previous fit results
	except:
		pass

	fit_name = "FitResults.txt"
	root_dir = os.path.split(results_dir)[0]
	mydir = os.path.join(directory.replace(root_dir, results_dir), str(dirnum))
	if not os.path.isdir(mydir):
		os.makedirs(mydir)

	pathfile = os.path.join(directory, file_name)
	imp = IJ.openImage(pathfile)
	if trim_track['run']: # trim the stack to the required frame number
		imp_frame = max(imp.getNSlices(), imp.getNFrames())
		if imp_frame > trim_track['frame_number']:
			s = SubstackMaker()
			trim_sting = "{}-{}".format(imp_frame-trim_track['frame_number']+1, imp_frame) if (trim_track['from_end']) else "{}-{}".format(1, trim_track['frame_number'])
			imp = s.makeSubstack(imp, trim_sting)
	imp.setTitle("Image")
	imp.show()
	IJ.run(imp, "Peak Fit", "template=[None] config_file=["+gdsc_smlm_xml+"] calibration="+str(pixel_size)+
	" gain=55.50 exposure_time="+str(frame_length)+" initial_stddev0=2.000 initial_stddev1=2.000 initial_angle=0.000 "+
	"smoothing=0.50 smoothing2=3 search_width=3 fit_solver=[Least Squares Estimator (LSE)] "+
	"fit_function=Circular local_background camera_bias=0 fit_criteria=[Least-squared error] significant_digits=5 coord_delta=0.0001 "+
	"lambda=10.0000 max_iterations=20 fail_limit=10 include_neighbours neighbour_height=0.30 +"
	"residuals_threshold=1 duplicate_distance=0.50 shift_factor=2 signal_strength="+
	str(signal_strength)+" width_factor=2 precision="+str(precision)+" min_photons="+str(min_photons)+" results_table=Uncalibrated "+
	"image=[Localisations (width=precision)] weighted equalised image_precision=5 image_scale="+str(sr_scale)+" results_dir=["+mydir+"]"+
	"local_background camera_bias=0 fit_criteria=[Least-squared error] significant_digits=5 coord_delta=0.0001 lambda=10.0000 "+
	"max_iterations=20 stack")
	
	sr = wm.getWindow("Image (LSE) SuperRes")
	IJ.selectWindow("Image (LSE) SuperRes")
	IJ.saveAs("Tiff", os.path.join(mydir, "SR_"+str(signal_strength)+"_"+str(precision)+"nm"+str(min_photons)+"photons.srf.tif"))
	IJ.run("Scale...", "x=0.125 y=0.125 width=512 height=512 interpolation=Bilinear average create")
	IJ.saveAs("Tiff", os.path.join(mydir, "SR_unscaled.srf.tif"))
	wm.getCurrentWindow().close()
	sr.close()
	imp.close()
	IJ.selectWindow("Fit Results")
	IJ.saveAs("text", os.path.join(mydir, fit_name))
	wm.getWindow("Fit Results").close()
	with open(os.path.join(mydir, "raw_img_{}.path.txt".format(os.path.splitext(file_name)[0])), 'w') as log:
		log.write(pathfile)

	return mydir

def cleanfi(f):
	f_set = set(f)
	L = []
	for coords in f_set:
		if f.count(coords) <= 10:
			L.append(coords)
	for _ in L: f_set.remove(_)
	return list(f_set)

def Find_fidicials(directory, fid_brightness, fid_size):
	to_load = os.path.join(directory, "FitResults.txt")
	dimensions = False
	fi = []
	with open(to_load, 'r') as fit_file:
		for lines in fit_file:
			if lines.startswith("Image:"):
				origX, origY, origValue = lines.split('\t')[2:5]
				if int(float(origValue)) > fid_brightness:
					fi.append((int(float(origX)), int(float(origY))))
				if not dimensions:
					x0 = lines.split('\t')[0].split('x0 y0 w')[1].split(' ')[0]
					y0 = lines.split('\t')[0].split('x0 y0 w')[1].split(' h')[1]
					dimensions = (int(x0), int(y0))
	if len(fi) < 100:
		return False
	else:
		fi = cleanfi(fi)
		if len(fi) == 0: return False
		imp = IJ.createImage("Fiducials", "8-bit black", dimensions[0], dimensions[1], 1)
		for x, y in fi:
			imp.setRoi(x, y, 1, 1)
			IJ.run(imp, "Add...", "value=255")
		imp.killRoi()
		IJ.run(imp, "Find Maxima...", "noise=10 output=[List]")
		result = IJ.getTextPanel()
		result.saveAs(os.path.join(directory, "Fiducials.txt"))
		result.clear()
		xylist = []
		with open(os.path.join(directory, "Fiducials.txt"), 'r') as fid:
			for lines in fid:
				if 'X' not in lines:
					xylist.append((int(lines.split('\t')[0]), int(lines.split('\t')[1])))
		
		rm = RoiManager.getInstance()
		if not rm:
			rm = RoiManager()
		rm.reset()
		for x, y in xylist:
			imp.setRoi(x, y, 1, 1)
			IJ.run(imp, "Enlarge...", "enlarge={}".format(fid_size)) 
			rm.addRoi(imp.getRoi())
			imp.killRoi()
		return True

def Correct_fidicials(directory, sr_scale):
	fit_name = "FitResults_Corrected.txt"
	IJ.run("Clear Memory Results", "All")
	fitresultfile = os.path.join(directory, "Image.results.xls")
	IJ.run("Results Manager", "coordinate=["+
	fitresultfile+"] input=File input_file=["+fitresultfile+"] results_table=Uncalibrated "+
	"image=[Localisations (width=precision)] weighted equalised image_precision=1.50 "+
	"image_scale="+str(sr_scale)+" image_window=0 results_in_memory")
	try:
		IJ.run("Drift Calculator", "input=[Image (LSE)] method=[Marked ROIs] "+
		"max_iterations=50 relative_error=0.010 smoothing=0.25 limit_smoothing "+
		"min_smoothing_points=10 max_smoothing_points=50 smoothing_iterations=1 "+
		"plot_drift update_method=[New dataset] save_drift "+
		"drift_file=["+os.path.join(directory, "Drift.txt")+"]")
	except:
		wm.getActiveWindow().close()
		return False
	wm.getWindow("Fit Results").close()
	IJ.run("Results Manager", "input=[Image (LSE) (Corrected)]  results_table=Uncalibrated "+
	"image=[Localisations (width=precision)] weighted equalised image_precision=1.50 "+
	"image_scale="+str(sr_scale)+" image_window=0 results_file=[] results_in_memory")
	IJ.selectWindow("Fit Results")
	IJ.saveAs("text", os.path.join(directory, fit_name))
	wm.getWindow("Fit Results").close()
	wm.getWindow("Drift X").close()
	wm.getWindow("Drift Y").close()
	wm.getWindow("Image (LSE) SuperRes").close()
	wm.getWindow("Image (LSE) (Corrected) SuperRes").close()
	
if __name__ in ['__builtin__', '__main__']:
	i = 0
	with open(jsonpath, 'r') as f:
		d = json.load(f)	# Configuration dictionary
	L = d['filelist']
	trim_track = d['trim_track']
	gdsc_smlm_xml = os.path.join(os.path.dirname(jsonpath), 'gdsc.smlm.settings.xml')
	results_dir = d['results_dir']
	pixel_size = d['pixel_size']
	frame_length = d['frame_length']
	signal_strength = d['signal_strength'] 
	precision = d['precision']
	min_photons = d['min_photons']
	sr_scale = d['sr_scale']
	run_fidicial_correction = d['fiducial_correction']['run']
	fid_brightness = d['fiducial_correction']['fid_brightness']
	fid_size = d['fiducial_correction']['fid_size']
	
#	with open(jsonpath, 'w') as f:
#		json.dump(d, f, indent=2) 	# Report back the results to python

	log_dict = {}
	for img in L:
		r = os.path.dirname(img)
		f = os.path.basename(img)
		try:
			mydir = GDSC_SMLM(r, f, results_dir, i, trim_track, gdsc_smlm_xml, pixel_size, frame_length, signal_strength, precision, min_photons, sr_scale)
		except Exception as e:
			log_dict[img] = {}
			log_dict[img]['GDSC_SMLM'] = 'Error: {}'.format(e) # GDSC_SMLM peakfit failed
		else:
			if run_fidicial_correction:
				try:
					flag = Find_fidicials(mydir, fid_brightness, fid_size)
					if flag:
						Correct_fidicials(mydir, sr_scale)
					else:
						log_dict[img] = {}
						log_dict[img]['Fiducials'] = 'none' # No fiducials detected
				except Exception as e:
					log_dict[img] = {}
					log_dict[img]['Fiducials'] = 'Error: {}'.format(e) # Correct fiducials failed
		finally:
			i += 1
	d['fit_name'] = "FitResults_Corrected.txt" if run_fidicial_correction else "FitResults.txt"
	d['ImageJ_log'] = log_dict
	with open(jsonpath, 'w') as f:
		json.dump(d, f, indent=2) 	# Report back the results to python

	IJ.run("Quit")
	
	
		
