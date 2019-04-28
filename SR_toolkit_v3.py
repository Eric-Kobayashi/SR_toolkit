from SR_toolkit_lib import Analysis
import time


'''
SR_toolkit_v3_3
This is the user runfile of SR_toolkit
For technical details please refer to SR_toolkit_lib.py

'''

def mainfunc():
    '''
    Input variables
    
    '''
       
    directory = r"C:\Users\Eric\OneDrive - University Of Cambridge\igorplay\feudicial"
    image_condition = lambda img: img.endswith('_561.tif')  # Change conditions to search for imagestacks
    input_dict = {
    # ==== Parameters for all analysis ====
    'pixel_size': 98.6 , # nm
    'sr_scale': 8 , # The scale used in length analysis and generation of super-res images 
    'camera_bias': 500.0 ,
    'camera_gain': 55.50 ,
    'frame_length': 50 , # ms per frame, i.e. exposure time
    'create_symlink_for_images': False, # Needs admin mode. Do not run if images are in Dropbox.
    'GDSC_SMLM_version': 1, # 1 or 2
    
    # ==== Parameters for GDSC SMLM fitting ====
    'GDSC_SMLM_peak_fit': True, # If False, the parameters for GDSC SMLM fitting will be ignored, only cluster_analysis will be run
    'trim_track': {'run': False, 'frame_number': 4000, 'from_end': True}, # trim the stack to the required frame number, from_end controls trim from end or beginning 
    'BG_measurement': {'run': True, 'correct_precision': False, 'kappa': 0.5}, # measure the background of the image_stack. (median pixel value) 
    # if correct_precision set to True, the fitting precision will be adjusted according to the background level: higher background -> higher precision threshold. kappa controls the extent of the adjustment.
    'signal_strength': 40, # Caution: Different measure in GDSCSMLM1 and GDSCSMLM 2
    'precision': 20.0, # nm
    'min_photons': 0,
    'fiducial_correction': {'run':True, 'fid_file':False, 'fid_brightness':50000,
    'fid_size': 6, 'fid_last_time':500, 'smoothing_para':0.25, 'limit_smoothing':True},
    # fid_file enables manually defined fiducials saved in 'Fiducials.txt', 
    # fid_brightness, fid_size, fid_last_time (frames) defines the criteria for automatically finding fiducials
    # smoothing_para and limit_smoothing controls the correction smoothing. If smoothing fails try setting limit_smotthing to False
    'spot_filter': {'filter_type': 'Difference', 'smoothing':0.20, 'smoothing2': 1.50}, # Only applicable to GDSC_SMLM_2, filter_type either Single or Difference,
    # smoothing, smoothing2 control the filter
    
    # ==== Parameters for cluster analysis and measurements ====
    'fitresults_file_name': 'default', # default: 'FitResults.txt' or 'FitResults_FeuRemoved.txt' if fiducial corrected
    'temporal_grouping': {'run': True, 'dThresh':20, 'min_loc':2, 'max_mol_area':float('inf'),
     'tThresh':2500, 'min_frame':2, 'min_burst':2, 'min_on_prop':0.0},
     # dThresh - nm, max_mol_area - nm2, tThresh - ms
    'cluster_analysis': {'run': True, 'cluster_subject':'burst', # burst, mol, loc
    'DBSCAN_eps_nm':200, 'DBSCAN_min_samples':2},
    # burst for STORM/PALM, mol for PAINT
    'length_measure': {'run':True, 'algorithm':'close', 'sigma':2}, # algorithm: blur or close; if blur, sigma is the gaussian sigma; if close, sigma is the closing square size
    'eccentricity_measure': True ,
    'convexhull_measure': True , # measure the area of the cluster
    'Rendering_SR': True , # rendering clustered SR images in imageJ
    'save_histogram':True,
    }

    Analysis(directory).run_fit(verbose=True, image_condition=image_condition, **input_dict)
if __name__ == '__main__':
    mainfunc()
#   time.sleep(1)