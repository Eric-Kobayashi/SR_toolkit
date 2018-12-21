from SR_toolkit_lib import Analysis

'''
SR_toolkit_v3_1
This is the user runfile of SR_toolkit
For technical details please refer to SR_toolkit_lib.py

'''

def mainfunc():
    '''
    Input variables
    
    '''
    directory = r"E:\Dropbox (Cambridge University)\Artemisia\Eric\SHSY5Y_Desab_imaging\Trimmed_images"
    image_condition = lambda img: img.endswith('_trimmed.tif')  # Change conditions to search for imagestacks
    input_dict = {
    # ==== Parameters for all analysis ====
    'pixel_size': 107 , # nm
    'sr_scale': 8 , # The scale used in length analysis and generation of super-res images 
    'frame_length': 50 , # ms per frame, i.e. exposure time
    'create_symlink_for_images': True, # Needs admin mode
    
    # ==== Parameters for GDSC SMLM fitting ====
    'trim_track': {'run': True, 'frame_number': 4000, 'from_end': True}, # trim the stack to the required frame number, from_end controls trim from end or beginning 
    'signal_strength': 150, 
    'precision': 30, # nm
    'min_photons': 0,
    'fiducial_correction': {'run':False, 'fid_brightness':20000},
    
    # ==== Parameters for cluster analysis and measurements ====
    'cluster_analysis_measurement': False, # If False, only GDSC fitting will be run
    'fitresults_file_name': 'default', # default: 'FitResults.txt' or 'FitResults_Corrected.txt' if fiducial corrected
    'DBSCAN_eps': 100 , # nm, not pixels!
    'DBSCAN_min_samples': 1 ,
    'burst_filter': {'run':True, 'fill_gap':50, 'min_burst':3} ,  # filter out the aggregates with less than min_bursts
    'burst_analysis': {'run':True, 'fill_gap':5, 'remove_single':True} ,  # analyse the burst number, burst length, dark length of the aggregates
    'length_measure': True ,
    'eccentricity_measure': True ,
    'convexhull_measure': True , # measure the area of the cluster
    'save_GDSC_header_file': True , # save GSDC header file for rendering, will consider remove since rendering will be incorporated in this code
    'save_histogram':True,
    
    # ==== Setup the analysis, only change it when first run ====
    'ImageJ_path': r"C:\Users\yz520\Downloads\fiji-win64-20170530\Fiji.app\ImageJ-win64.exe" ,
    'GDSCSMLM_script_path': r"C:\Users\yz520\Desktop\OneDrive - University Of Cambridge\Data analysis package\My Library\Jython\Peakfit_GDSC_SMLM.py" ,
    'json_file': "C://Users//yz520//to_ImageJ.json" # Needs to be // otherwise it will fail !
    }

    Analysis(directory).run_fit(verbose=True, image_condition=image_condition, **input_dict)

if __name__ == '__main__':
    mainfunc()