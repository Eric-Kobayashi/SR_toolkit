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
    'camera_bias': 500 ,
    'camera_gain': 55.50 ,
    'frame_length': 50 , # ms per frame, i.e. exposure time
    'create_symlink_for_images': False, # Needs admin mode. Do not run if images are in Dropbox.
    
    # ==== Parameters for GDSC SMLM fitting ====
    'GDSC_SMLM_peak_fit': True, # If False, the parameters for GDSC SMLM fitting will be ignored, only cluster_analysis will be run
    'trim_track': {'run': False, 'frame_number': 4000, 'from_end': True}, # trim the stack to the required frame number, from_end controls trim from end or beginning 
    'BG_measurement': {'run': True, 'correct_precision': False, 'kappa': 0.5}, # measure the background of the image_stack. (median pixel value) 
    # if correct_precision set to True, the fitting precision will be adjusted according to the background level: higher background -> higher precision threshold. kappa controls the extent of the adjustment.
    'signal_strength': 100, 
    'precision': 20, # nm
    'min_photons': 0,
    'fiducial_correction': {'run':True, 'fid_file':True, 'fid_brightness':20000,
    'fid_size': 6, 'fid_last_time':500, 'smoothing_para':0.25, 'limit_smoothing':True},
    # fid_file enables manually defined fiducials saved in 'Fiducials.txt', 
    # fid_brightness, fid_size, fid_last_time defines the criteria for automatically finding fiducials
    # smoothing_para and limit_smoothing controls the correction smoothing. If smoothing fails try setting limit_smotthing to False
    
    # ==== Parameters for cluster analysis and measurements ====
    'cluster_analysis_measurement': False, # If False, only GDSC fitting will be run
    'fitresults_file_name': 'default', # default: 'FitResults.txt' or 'FitResults_FeuRemoved.txt' if fiducial corrected
    'DBSCAN_eps': 100 , # nm, not pixels!
    'DBSCAN_min_samples': 5 ,
    'burst_filter': {'run':False, 'fill_gap':50, 'min_burst':3} ,  # filter out the aggregates with less than min_bursts
    'burst_analysis': {'run':False, 'fill_gap':5, 'remove_single':True} ,  # analyse the burst number, burst length, dark length of the aggregates
    'length_measure': {'run':False, 'algorithm':'close', 'sigma':2}, # algorithm: blur or close; if blur, sigma is the gaussian sigma; if close, sigma is the closing square size
    'eccentricity_measure': False ,
    'convexhull_measure': False , # measure the area of the cluster
    'Rendering_SR': True , # rendering clustered SR images in imageJ
    'save_histogram':True,
    
    # ==== Setup the analysis, only change it when first run ====
    'ImageJ_path': r"C:\Users\Eric\fiji-win64\Fiji.app\ImageJ-win64.exe" ,
    'GDSCSMLM_script_path': r"C:\Users\Eric\Documents\data analysis package\Jython" ,
    'json_file': "C://Users//Eric//to_ImageJ.json" # Needs to be // otherwise it will fail !
    }

    Analysis(directory).run_fit(verbose=True, image_condition=image_condition, **input_dict)

if __name__ == '__main__':
    mainfunc()
#   time.sleep(1)