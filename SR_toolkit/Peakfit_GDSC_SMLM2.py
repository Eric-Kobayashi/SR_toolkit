# @String jsonpath

from __future__ import with_statement
import os
import os.path as op
import shutil as sh
import math
from ij import IJ
from ij import WindowManager as wm
import ij.plugin.frame.RoiManager as RoiManager
import ij.plugin.ZProjector as ZProjector
import json
from ij.plugin import SubstackMaker
import subprocess
import time

def GDSC_SMLM(directory, file_name, mydir, trim_track, bg_measurement, spot_filter, pixel_size, camera_gain, camera_bias, frame_length, signal_strength, precision, min_photons, sr_scale):
	try:
		wm.getWindow("Fit Results").close()	# clean the previous fit results
	except:
		pass

	fit_name = "FitResults.txt"
	pathfile = op.join(directory, file_name)
	imp = IJ.openImage(pathfile)
	h, w = imp.getHeight(), imp.getWidth()
	if trim_track['run']: # trim the stack to the required frame number
		imp_frame = max(imp.getNSlices(), imp.getNFrames())
		if imp_frame > trim_track['frame_number']:
			s = SubstackMaker()
			trim_sting = "{}-{}".format(imp_frame-trim_track['frame_number']+1, imp_frame) if (trim_track['from_end']) else "{}-{}".format(1, trim_track['frame_number'])
			imp = s.makeSubstack(imp, trim_sting)
	if bg_measurement:
		IJ.run(imp, "Hyperstack to Stack", "")
		ave_imp = ZProjector.run(imp, "avg")
		stats = ave_imp.getStatistics(0x10000)
		with open(op.join(mydir, "median_intensity.txt"), 'w') as f:
			f.write(str(stats.median))

	imp.setTitle("Image")
	imp.show()
	IJ.redirectErrorMessages()
	IJ.run(imp, "Peak Fit", "template=[None] camera_type=EMCCD calibration="+str(pixel_size)+" exposure_time="+str(frame_length)+
	" psf=[Circular Gaussian 2D] spot_filter_type="+spot_filter['filter_type']+" spot_filter=Gaussian smoothing="+str(spot_filter['smoothing'])+" search_width=3 border_width=0.50 fitting_width=3 "+
	"fit_solver=[LVM LSE] fail_limit=10 pass_rate=0.50 include_neighbours neighbour_height=0.30 residuals_threshold=1 "+
	"duplicate_distance=0.50 shift_factor=2 signal_strength="+str(signal_strength)+" min_photons="+str(min_photons)+" min_width_factor=0.20 width_factor=2 precision="+str(precision)+" "+
	"show_results_table image=[Localisations (width=precision)] results_format=Text results_directory=["+mydir+"] "+
	"save_to_memory camera_bias="+str(camera_bias)+" gain="+str(camera_gain)+" read_noise=0.0000 psf_parameter_1=1.000 duplicate_distance_absolute "+
	"precision_method=Mortensen table_distance_unit=[pixel (px)] table_intensity_unit=[unknown (na)] table_angle_unit=[unknown (na)] "+
	"table_show_fitting_data table_show_noise_data table_show_precision table_precision=4 equalised image_scale="+str(sr_scale)+" "+
	"file_distance_unit=[pixel (px)] file_intensity_unit=[unknown (na)] file_angle_unit=[unknown (na)] "+
	"spot_filter2=Gaussian smoothing2="+str(spot_filter['smoothing2'])+" relative_threshold=1.0E-6 absolute_threshold=1.0E-16 parameter_relative_threshold=0.001 "+
	"parameter_absolute_threshold=1.0E-6 max_iterations=20 lambda=10.0000 stack");
	
	sr = wm.getWindow("Image (LVM LSE) SuperRes")
	IJ.selectWindow("Image (LVM LSE) SuperRes")
	IJ.saveAs("Tiff", op.join(mydir, "RawSR_"+str(signal_strength)+"_{:.2f}nm".format(precision)+str(min_photons)+"photons.srf.tif"))
	IJ.run("Scale...", "x={0} y={0} width={1} height={2} interpolation=Bilinear average create".format(1.0/sr_scale, w, h))
	IJ.saveAs("Tiff", op.join(mydir, "RawSR_unscaled.srf.tif"))
	wm.getCurrentWindow().close()
	sr.close()
	imp.close()
	IJ.selectWindow("Fit Results")
	IJ.saveAs("text", op.join(mydir, fit_name))
	wm.getWindow("Fit Results").close()
	with open(op.join(mydir, "raw_img_{}.path.txt".format(op.splitext(file_name)[0])), 'w') as log:
		log.write(pathfile)
	return w, h

def cleanfi(f, fid_last_time):
	f_set = set(f)
	L = []
	for coords in f_set:
		if f.count(coords) <= fid_last_time:  # this is a variable parameter
			L.append(coords)
	for _ in L: f_set.remove(_)
	return list(f_set)

def Find_fidicials(directory, fid_brightness, fid_size, fid_last_time):
	to_load = op.join(directory, "FitResults.txt")
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
		fi = cleanfi(fi, fid_last_time)
		if len(fi) == 0: return False
		imp = IJ.createImage("Fiducials", "8-bit black", dimensions[0], dimensions[1], 1)
		for x, y in fi:
			imp.setRoi(x, y, 1, 1)
			IJ.run(imp, "Add...", "value=255")
		imp.killRoi()
		IJ.run("Input/Output...", "jpeg=85 gif=-1 file=.txt use_file copy_row save_column")
		IJ.run(imp, "Find Maxima...", "noise=10 output=[List]")
		result = IJ.getTextPanel()
		result.saveAs(op.join(directory, "Fiducials.txt"))
		result.clear()
		xylist = []
		with open(op.join(directory, "Fiducials.txt"), 'r') as fid:
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

def Correct_fidicials_with_drift(img_dir, directory, sr_scale, w, h, fid_size, smoothing_para, limit_smoothing):	
	fit_name = "FitResults_Corrected.txt"
	IJ.run("Clear Memory Results", "All")
	fitresultfile = op.join(directory, "Image.results.xls")
	IJ.redirectErrorMessages()
	IJ.run("Results Manager", "input=File input_file=["+fitresultfile+"] table=ImageJ "+
	"image=[Localisations (width=precision)] results_format=None results_file=[] save_to_memory "+
	"table_distance_unit=[pixel (px)] table_intensity_unit=photon table_angle_unit=[unknown (na)] "+
	"table_show_fitting_data table_show_noise_data table_show_precision table_precision=4 "+
	"equalised image_scale="+str(sr_scale))
	drift_file = op.join(directory, "Drift.txt")
	try:
		sh.copyfile(op.join(img_dir, "Drift.txt"), drift_file)
	except:
		IJ.run("Close All", "")
		return False
	try:
		IJ.redirectErrorMessages()
		IJ.run("Drift Calculator", "input=[Image (LVM LSE)] method=[Drift File] "+
		"max_iterations=50 relative_error=0.010 smoothing="+str(smoothing_para)+
		" {}".format("limit_smoothing " if limit_smoothing else "")+
		"min_smoothing_points=10 max_smoothing_points=50 smoothing_iterations=1 "+
		"update_method=[New dataset] save_drift "+
		"drift_file=["+drift_file+"]")
	except:
		IJ.run("Close All", "")
		return False

	if not op.isfile(drift_file):
		IJ.run("Close All", "")
		return False
		
	wm.getWindow("Fit Results").close()
	IJ.redirectErrorMessages()
	IJ.run("Results Manager", "input=[Image (LVM LSE) (Corrected)] table=ImageJ "+
	"image=[Localisations (width=precision)] results_format=None results_file=[] save_to_memory "+
	"table_distance_unit=[pixel (px)] table_intensity_unit=photon table_angle_unit=[unknown (na)] "+
	"table_show_fitting_data table_show_noise_data table_show_precision table_precision=4 "+
	"equalised image_scale="+str(sr_scale))

	IJ.selectWindow("Fit Results")
	IJ.saveAs("text", op.join(directory, fit_name))
	IJ.run("Close All", "")
	return True
		
		
def Correct_fidicials_with_fid(img_dir, directory, sr_scale, w, h, fid_size, smoothing_para, limit_smoothing):
	xylist = []
	fid_file = op.join(img_dir, "Fiducials.txt")
	if not op.isfile(fid_file):
		return False
	with open(fid_file, 'r') as fid:
		for lines in fid:
			if 'X' not in lines:
				xylist.append((int(lines.split('\t')[0]), int(lines.split('\t')[1])))
	if len(xylist) < 1:
		return False
	sh.copyfile(fid_file, op.join(directory, "Fiducials.txt"))
	rm = RoiManager.getInstance()
	if not rm:
		rm = RoiManager()
	rm.reset()
	imp = IJ.createImage("Fiducials", "8-bit black", w, h, 1)
	for x, y in xylist:
		imp.setRoi(x, y, 1, 1)
		IJ.run(imp, "Enlarge...", "enlarge={}".format(fid_size)) 
		rm.addRoi(imp.getRoi())
		imp.killRoi()
		
	fit_name = "FitResults_Corrected.txt"
	IJ.run("Clear Memory Results", "All")
	fitresultfile = op.join(directory, "Image.results.xls")
	IJ.redirectErrorMessages()
	IJ.run("Results Manager", "input=File input_file=["+fitresultfile+"] table=ImageJ "+
	"image=[Localisations (width=precision)] results_format=None results_file=[] save_to_memory "+
	"table_distance_unit=[pixel (px)] table_intensity_unit=photon table_angle_unit=[unknown (na)] "+
	"table_show_fitting_data table_show_noise_data table_show_precision table_precision=4 "+
	"equalised image_scale="+str(sr_scale))
	
	drift_file = op.join(directory, "Drift.txt")
	try:
		IJ.redirectErrorMessages()
		IJ.run("Drift Calculator", "input=[Image (LVM LSE)] method=[Marked ROIs] "+
		"max_iterations=50 relative_error=0.010 smoothing="+str(smoothing_para)+
		" {}".format("limit_smoothing " if limit_smoothing else "")+
		"min_smoothing_points=10 max_smoothing_points=50 smoothing_iterations=1 "+
		"update_method=[New dataset] save_drift "+
		"drift_file=["+drift_file+"]")
	except:
		IJ.run("Close All", "")
		return False

	if not op.isfile(drift_file):
		IJ.run("Close All", "")
		return False
		
	wm.getWindow("Fit Results").close()
	IJ.redirectErrorMessages()
	IJ.run("Results Manager", "input=[Image (LVM LSE) (Corrected)] table=ImageJ "+
	"image=[Localisations (width=precision)] results_format=None results_file=[] save_to_memory "+
	"table_distance_unit=[pixel (px)] table_intensity_unit=photon table_angle_unit=[unknown (na)] "+
	"table_show_fitting_data table_show_noise_data table_show_precision table_precision=4 "+
	"equalised image_scale="+str(sr_scale))

	IJ.selectWindow("Fit Results")
	IJ.saveAs("text", op.join(directory, fit_name))
	IJ.run("Close All", "")
	return True
		
def Correct_fidicials(directory, sr_scale, smoothing_para, limit_smoothing):
	fit_name = "FitResults_Corrected.txt"
	IJ.run("Clear Memory Results", "All")
	fitresultfile = op.join(directory, "Image.results.xls")
	IJ.redirectErrorMessages()
	IJ.run("Results Manager", "input=File input_file=["+fitresultfile+"] table=ImageJ "+
	"image=[Localisations (width=precision)] results_format=None results_file=[] save_to_memory "+
	"table_distance_unit=[pixel (px)] table_intensity_unit=photon table_angle_unit=[unknown (na)] "+
	"table_show_fitting_data table_show_noise_data table_show_precision table_precision=4 "+
	"equalised image_scale="+str(sr_scale))
	
	drift_file = op.join(directory, "Drift.txt")
	try:
		IJ.redirectErrorMessages()
		IJ.run("Drift Calculator", "input=[Image (LVM LSE)] method=[Marked ROIs] "+
		"max_iterations=50 relative_error=0.010 smoothing="+str(smoothing_para)+
		" {}".format("limit_smoothing " if limit_smoothing else "")+
		"min_smoothing_points=10 max_smoothing_points=50 smoothing_iterations=1 "+
		"update_method=[New dataset] save_drift drift_file=["+drift_file+"]")
	except:
		IJ.run("Close All", "")
		return False

	if not op.isfile(drift_file):
		IJ.run("Close All", "")
		return False
		
	wm.getWindow("Fit Results").close()
	IJ.redirectErrorMessages()
	IJ.run("Results Manager", "input=[Image (LVM LSE) (Corrected)] table=ImageJ "+
	"image=[Localisations (width=precision)] results_format=None results_file=[] save_to_memory "+
	"table_distance_unit=[pixel (px)] table_intensity_unit=photon table_angle_unit=[unknown (na)] "+
	"table_show_fitting_data table_show_noise_data table_show_precision table_precision=4 "+
	"equalised image_scale="+str(sr_scale))
	
	IJ.selectWindow("Fit Results")
	IJ.saveAs("text", op.join(directory, fit_name))
	IJ.run("Close All", "")
	return True

# def create_image_pixel_median(img, savedir): 
# 	imp = IJ.openVirtual(img)
# 	IJ.run(imp, "Hyperstack to Stack", "")
# 	ave_imp = ZProjector.run(imp, "avg")
# 	stats = ave_imp.getStatistics(0x10000)
# 	IJ.run("Close All", "")
# 	with open(op.join(savedir, "median_intensity.txt"), 'w') as f:
# 		f.write(str(stats.median))
# 	return float(stats.median)
	
# def getMedian(numericValues):
# 	theValues = sorted(numericValues)
	
# 	if len(theValues) % 2 == 1:
# 		return float(theValues[(len(theValues)+1)/2-1])
# 	else:
# 		lower = theValues[len(theValues)/2-1]
# 		upper = theValues[len(theValues)/2]
# 		return (float(lower + upper)) / 2 
	
if __name__ in ['__builtin__', '__main__']:
	with open(jsonpath, 'r') as f:
		d = json.load(f)	# Configuration dictionary
	L = d['filelist']
	trim_track = d['trim_track']
	bg_measurement = d['BG_measurement']
	gdsc_smlm_xml = op.join(op.dirname(jsonpath), 'gdsc.smlm.settings.xml')
	results_dir = d['results_dir']
	pixel_size = d['pixel_size']
	camera_bias = d['camera_bias']
	camera_gain = d['camera_gain']
	frame_length = d['frame_length']
	signal_strength = d['signal_strength'] 
	precision = d['precision']
	min_photons = d['min_photons']
	sr_scale = d['sr_scale']
	run_fidicial_correction = d['fiducial_correction']['run']
	correction_method = d['fiducial_correction']['correction_method']
	
	fid_brightness = (d['fiducial_correction']['fid_brightness']*1.0)/camera_gain
	fid_size = d['fiducial_correction']['fid_size']
	fid_last_time = d['fiducial_correction']['fid_last_time']
	smoothing_para = d['fiducial_correction']['smoothing_para']
	limit_smoothing = d['fiducial_correction']['limit_smoothing']
	spot_filter = d['spot_filter']
	
#	with open(jsonpath, 'w') as f:
#		json.dump(d, f, indent=2) 	# Report back the results to python

	log_dict = {}
	i = 0
	List_of_mydir = []

	for img in L:
		r = op.dirname(img)
		f = op.basename(img)
		root_dir = op.dirname(results_dir)
		mydir = op.join(r.replace(root_dir, results_dir), str(i))
		if not op.isdir(mydir):
			os.makedirs(mydir)
		List_of_mydir.append((img, mydir))
		i += 1

	for img, mydir in List_of_mydir:
		r = op.dirname(img)
		f = op.basename(img)
		if bg_measurement:
			corrected_precision = pre_dict[img]
		else:
			corrected_precision = precision
		try:
			w, h = GDSC_SMLM(r, f, mydir, trim_track, bg_measurement, spot_filter, pixel_size, camera_gain, camera_bias, frame_length, signal_strength, corrected_precision, min_photons, sr_scale)
		except Exception as e:
			log_dict[img] = {}
			log_dict[img]['GDSC_SMLM'] = 'Error: {}'.format(e) # GDSC_SMLM peakfit failed
		else:
			if run_fidicial_correction:
				if correction_method == 'fid_file':
					try:
						flag = Correct_fidicials_with_fid(r, mydir, sr_scale, w, h, fid_size, smoothing_para, limit_smoothing)
						if not flag:
							log_dict[img] = {}
							log_dict[img]['Fiducials'] = 'none' # No fiducials detected
					except Exception as e:
						log_dict[img] = {}
						log_dict[img]['Fiducials'] = 'Error: {}'.format(e) # Correct fiducials failed
				
				elif correction_method == 'drift':
					try:
						flag = Correct_fidicials_with_drift(r, mydir, sr_scale, w, h, fid_size, smoothing_para, limit_smoothing)
						if not flag:
							log_dict[img] = {}
							log_dict[img]['Fiducials'] = 'none' # No fiducials detected
					except Exception as e:
						log_dict[img] = {}
						log_dict[img]['Fiducials'] = 'Error: {}'.format(e) # Correct fiducials failed
				
				elif correction_method == 'auto_fid':
					try:
						flag = Find_fidicials(mydir, fid_brightness, fid_size, fid_last_time)
						if flag:
							Correct_fidicials(mydir, sr_scale, smoothing_para, limit_smoothing)
						else:
							log_dict[img] = {}
							log_dict[img]['Fiducials'] = 'none' # No fiducials detected
					except Exception as e:
						log_dict[img] = {}
						log_dict[img]['Fiducials'] = 'Error: {}'.format(e) # Correct fiducials failed

	d['fit_name'] = "FitResults_Corrected.txt" if run_fidicial_correction else "FitResults.txt"
	d['ImageJ_log'] = log_dict
	with open(jsonpath, 'w') as f:
		json.dump(d, f, indent=2) 	# Report back the results to python

	IJ.run("Quit")
	
	
		
