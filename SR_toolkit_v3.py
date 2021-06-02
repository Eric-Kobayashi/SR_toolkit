from SR_toolkit_lib import Analysis
import time


'''
SR_toolkit_v3_6
Please refer to README.md for parameter settings

'''

def mainfunc():
    '''
    Input variables
    
    '''
       
    directory = r"C:\Users\Eric\Documents\SR_toolkit_reboot\561"
    image_condition = lambda img: img.endswith('.tif')
    input_dict = {
    # ==== Parameters for all analysis ====
    'pixel_size': 107 , # nm
    'sr_scale': 8 , 
    'camera_bias': 400.0 ,
    'camera_gain': 84.40 ,
    'frame_length': 50 , # ms per frame, i.e. exposure time
    'GDSC_SMLM_version': 1, # 1 or 2
    
    # ==== Parameters for GDSC SMLM fitting ====
    'GDSC_SMLM_peak_fit': True, 
    'trim_track': {'run': False, 'frame_number': 6000, 'from_end': True}, 
    'BG_measurement': True, # measure the background of the image_stack. (median pixel value) 
    'signal_strength': 40, # Caution: Different measure in GDSCSMLM1 and GDSCSMLM 2
    'precision': 20.0, # nm
    'min_photons': 0,
    'spot_filter': {'filter_type': 'Difference', 'smoothing':0.20, 'smoothing2': 1.50}, # Only applicable to GDSC_SMLM_2, filter_type either Single or Difference,
    # smoothing, smoothing2 control the filter
    'fiducial_correction': {'run':False, 'correction_method':'auto_fid', # correction_method: 'auto_fid', 'fid_file', 'drift', 'corr' 
    'bin_size': 20, 'segpara':500, # parameters for 'corr' - correlation correction, bin_size in nm, segpara in frames
    'fid_brightness':10000, 'fid_size': 6, 'fid_last_time':500, # parameters for 'auto_fid' - automatically finding fiducials, fid_last_time in frames
    'smoothing_para':0.25, 'limit_smoothing':True},  

    
    # ==== Parameters for cluster analysis and measurements ====
    'fitresults_file_name': 'default', 
    'temporal_grouping': {'run': True, 'dThresh':15, 'min_loc':2, 'max_mol_area':float('inf'),
     'tThresh':101, 'min_frame':2, 'min_burst':2, 'min_on_prop':0.0},
     # dThresh - nm, max_mol_area - nm2, tThresh - ms
    'cluster_analysis': {'run': True, 'cluster_subject':'burst', # burst, loc
    'DBSCAN_eps_nm':75, 'DBSCAN_min_samples':2},
    'length_measure': {'run':True, 'algorithm':'close', 'sigma':2}, 
    'eccentricity_measure': True ,
    'convexhull_measure': True , 
    'Rendering_SR': 'loc_pre' , # loc, loc_pre, or False
    'save_histogram':True,
    }

    Analysis(directory).run_fit(verbose=True, image_condition=image_condition, **input_dict)
if __name__ == '__main__':
    mainfunc()
