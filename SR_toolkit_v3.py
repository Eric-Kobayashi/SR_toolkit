import SR_toolkit_lib

'''
This is the user runfile of SR_toolkit
For technical details please refer to SR_toolkit_lib.py

'''

def mainfunc():
    '''
    Input variables
    
    '''
    directory = r"C:\Users\yz520\Desktop\OneDrive - University Of Cambridge\igorplay\training_set2"
    fit_filename = 'FitResults_feu_removed.txt'
    
    input_dict = {
    'pixel_size': 97.7 , # nm
    'DBSCAN_eps': 100 , # nm, not pixels!
    'DBSCAN_min_samples': 10 ,
    'sr_scale': 8 , # The scale used in length analysis and generation of super-res images 
    'frame_length': 50 , # ms, frame 
    'burst_filter': {'run':True, 'fill_gap':50, 'min_burst':2} ,
        # filter out the aggregates with less than min_bursts
    'burst_analysis': {'run':True, 'fill_gap':10, 'remove_single':True} ,
        # Analyse the burst number, burst length, dark length of the aggregates
    'length_measure': True ,
    'eccentricity_measure': True ,
    'convexhull_measure': True , # Measure the area of the cluster
    'copy_analysis': True , # copy skeleton images, dbscan_results into the analysis folder
    'save_GDSC_header_file': True , # save 
    'save_histogram':True}

    SR_toolkit_lib.Analysis(directory).run_fit(fit_name=fit_filename, verbose=True, **input_dict)


if __name__ == '__main__':
    mainfunc()