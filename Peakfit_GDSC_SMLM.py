from __future__ import with_statement
import os
import shutil as sh
from ij import IJ
from ij import WindowManager as wm
import ij.plugin.frame.RoiManager as RoiManager
import ij.plugin.ZProjector as ZProjector
from ij.gui import YesNoCancelDialog as dia

gdsc_smlm_xml_location = "C://Users//yz520//gdsc.smlm.settings.xml"
Pixel_size = 237 # nm  # Artemisia pixel size
Signal_strength = 40 
Precision = 80 # nm
Min_photons = 0

def GDSC_SMLM(directory, file_name, dirnum, gdsc_smlm_xml_location, Pixel_size, Signal_strength, Precision, Min_photons):
	try:	
		wm.getWindow("Fit Results").close()	# clean the previous fit results
	except:
		pass

	mydir = os.path.join(directory, str(dirnum))
	try:
		os.mkdir(mydir)
	except:
		#sh.rmtree(mydir)
		raise OSError("Directory already exists, try clearing the previous analysis")
	pathfile = os.path.join(directory, file_name)
	imp = IJ.openVirtual(pathfile)
	imp.setTitle("Image")
	
	IJ.run(imp, "Peak Fit", "config_file=["+gdsc_smlm_xml_location+"] calibration="+str(Pixel_size)+
	" gain=55.50 exposure_time=50 initial_stddev0=2.000 initial_stddev1=2.000 initial_angle=0.000 "+
	"smoothing=0.50 smoothing2=3 search_width=3 fit_solver=[Least Squares Estimator (LSE)] "+
	"fit_function=Circular local_background camera_bias=0 fit_criteria=[Least-squared error] significant_digits=5 coord_delta=0.0001 "+
	"lambda=10.0000 max_iterations=20 fail_limit=10 include_neighbours neighbour_height=0.30 +"
	"residuals_threshold=1 duplicate_distance=0.50 shift_factor=2 signal_strength="+
	str(Signal_strength)+" width_factor=2 precision="+str(Precision)+" min_photons="+str(Min_photons)+" results_table=Uncalibrated "+
	"image=[Localisations (width=precision)] weighted equalised image_precision=5 image_scale=8 results_dir=["+mydir+"]"+
	"local_background camera_bias=0 fit_criteria=[Least-squared error] significant_digits=5 coord_delta=0.0001 lambda=10.0000 "+
	"max_iterations=20 stack")
	
	sr = wm.getWindow("Image (LSE) SuperRes")
	IJ.selectWindow("Image (LSE) SuperRes")
	IJ.saveAs("Tiff", mydir+"//"+"SR_"+f.split('.')[0]+"_"+str(Signal_strength)+"_"+str(Precision)+"nm"+str(Min_photons)+"photons.tif")
#		IJ.run("Scale...", "x=0.125 y=0.125 width=512 height=512 interpolation=Bilinear average create")
#		IJ.saveAs("Tiff", mydir+"//"+"SR_"+f.split('.')[0]+"_"+str(Signal_strength)+"_"+str(Precision)+"nm"+"_unscaled.tif")
#		wm.getCurrentWindow().close()
	sr.close()
	IJ.selectWindow("Fit Results")
	IJ.saveAs("text", mydir+"//"+"FitResults.txt")
	wm.getWindow("Fit Results").close()
	with open(os.path.join(mydir, "raw_file_path.txt"), 'w') as log:
		log.write(pathfile)

def cleanFeu(f):
	f_set = set(f)
	L = []
	for coords in f_set:
		if f.count(coords) <= 10:
			L.append(coords)
	for _ in L: f_set.remove(_)
	return list(f_set)

def Find_Feudicials(directory):
	to_load = os.path.join(directory, "FitResults.txt")
	dimensions = False
	Feu = []
	with open(to_load, 'r') as fit_file:
		for lines in fit_file:
			if lines.startswith("Image:"):
				origX, origY, origValue = lines.split('\t')[2:5]
				if int(float(origValue)) > 3000:
					Feu.append((int(float(origX)), int(float(origY))))
				if not dimensions:
					x0 = lines.split('\t')[0].split('x0 y0 w')[1].split(' ')[0]
					y0 = lines.split('\t')[0].split('x0 y0 w')[1].split(' h')[1]
					dimensions = (int(x0), int(y0))
	if len(Feu) < 100:
		return False
	else:
		Feu = cleanFeu(Feu)
		if len(Feu) == 0: return False
		imp = IJ.createImage("Feducials", "8-bit black", dimensions[0], dimensions[1], 1)
		for x, y in Feu:
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
			IJ.run(imp, "Enlarge...", "enlarge=6") 
			rm.addRoi(imp.getRoi())
			imp.killRoi()
		return True

def Correct_Feudicials(directory):
	IJ.run("Clear Memory Results", "All")
	fitresultfile = os.path.join(directory, "Image.results.xls")
	IJ.run("Results Manager", "coordinate=["+
	fitresultfile+"] input=File input_file=["+fitresultfile+"] results_table=Uncalibrated "+
	"image=[Localisations (width=precision)] weighted equalised image_precision=1.50 "+
	"image_scale=8 image_window=0 results_in_memory")
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
	"image_scale=8 image_window=0 results_file=[] results_in_memory")
	IJ.selectWindow("Fit Results")
	IJ.saveAs("text", os.path.join(directory, "FitResults_Corrected.txt"))
	wm.getWindow("Fit Results").close()
	wm.getWindow("Drift X").close()
	wm.getWindow("Drift Y").close()
	wm.getWindow("Image (LSE) SuperRes").close()
	wm.getWindow("Image (LSE) (Corrected) SuperRes").close()
	
if __name__ == '__builtin__':
	L = []
	i = 0
	pa = r"E:\Dropbox (Cambridge University)\Dropbox for Aloe vera\Jason\20181025\Image sequence combined_PrP_mE10 with beads"
	
	for roots, dirs, files in os.walk(pa):
		for f in files:
			if f.endswith('.tif') and '488' in f:
				L.append((roots, f))
	
	for r, f in L:
		GDSC_SMLM(r, f, i, gdsc_smlm_xml_location, Pixel_size, Signal_strength, Precision, Min_photons)
		flag = Find_Feudicials(os.path.join(r, str(i)))
		if flag:
			Correct_Feudicials(os.path.join(r, str(i)))
		i += 1
		
