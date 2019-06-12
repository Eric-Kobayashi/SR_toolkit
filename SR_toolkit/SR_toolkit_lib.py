'''
SR_toolkit_v3_5

Issued: 28-04-2019
@author: Eric Kobayashi

Changes over v3_4:
    
    Major changes:
        
    1. Add temporal grouping: group localisations into molecules, then group
    them into individual bursts (using 2D and 1D DBSCAN, respectively.) Filters
    can be applied on burst number of molecule, molecule size and on time
    proportion of burst (blinking molecules can be removed). The X, Y positions
    and precisions in each group are averaged to enable further clustering.
    
    2. Add DBSCAN clustering with molecules or bursts. DBSCAN with bursts will 
    be suitable with STORM/PALM, while molecule suits with PAINT.
    
    3. Remove the burst filter and burst analysis because temporal grouping 
    does the same job.
    
    Minor changes:
    
    1. Move to_ImageJ.json to individual analysis folders, allowing multiple 
    analysis conducted at the same time. 
    
    2. Simplify the compatiablity code with GDSC SMLM 2.
    
    3. Simplify the data saving code.
    
'''

import os
import os.path as op
from shutil import copyfile, SameFileError
import sys
from datetime import datetime
from collections import OrderedDict
import numpy as np
import pandas as pd
from pandas import DataFrame as DF
import json
from copy import deepcopy
from subprocess import call
import warnings 
import ctypes
from SR_fit_lib import SR_fit

class Analysis(object):
    '''
    This class is used to package the analysis together.
    
    '''
    def __init__(self, path):
        self.path = path
        timestring = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        self.timestring = timestring
        self.results_dir = op.join(path, 'Analysis_'+timestring)
        os.mkdir(self.results_dir)
        sys.stdout.write("Time: {} Analysis starts!\n\n".format(timestring))
        self.json_file = op.join(self.results_dir, 'to_imageJ.json')
        self.json_log = OrderedDict()
        self.json_log['Time'] = timestring
        self.json_log['Path'] = path
        current_path = op.dirname(op.abspath(__file__))
        with open(op.join(current_path, 'config.txt'), 'r') as f:
            self.config_d = json.load(f)
        
    def _peak_fit_info(self, **kwargs):
        '''
        Return all the parameters related to the ImageJ peak fit and fiducial
        correction.
        '''
        
        peak_fit_para = {k: kwargs[k] for k in ['pixel_size', 'sr_scale', 
        'frame_length', 'camera_bias', 'camera_gain', 'BG_measurement', 'trim_track', 
        'signal_strength', 'precision', 'min_photons', 'fiducial_correction',
        'GDSC_SMLM_version', 'spot_filter']}
        peak_fit_para['filelist'] = self.imglist
        return peak_fit_para
        
    def _existing_peak_fit(self):
        '''
        Search all the existing peak fit in the folder.
        '''
        List_of_peak_fit_info_file = []
        Dict_of_peak_fit_info = {}
        for roots, dirs, files in os.walk(self.path):
            for f in files:
                if f == 'peak_fit_log.json':
                    List_of_peak_fit_info_file.append((roots, f))
        for json_root, json_log in List_of_peak_fit_info_file:
            with open(op.join(json_root, json_log), 'r') as f:
                Dict_of_peak_fit_info[json_root] = json.load(f)
        return Dict_of_peak_fit_info
            
    def _log_fit_info(self, fit_name, **kwargs):
        '''
        Write all the peak fit parameters into a config JSON file.
        To record the peak fit has been done and avoid repeating.
        Add fit_name to pass to cluster_fit function.
        '''
        peak_fit_para = self._peak_fit_info(**kwargs)
        peak_fit_para['fit_name'] = fit_name
        log_file = op.join(self.results_dir, 'peak_fit_log.json')
        with open(log_file, 'w') as f:
            json.dump(peak_fit_para, f, indent=2)

    def _search_images(self, condition):
        '''
        Search all the tif images that meet the condition.
        '''
        List_of_images = []
        for roots, dirs, files in os.walk(self.path):
            for f in files:
                if condition(f) and not f.endswith('.srf.tif'):  
                    #.srf.tif is the SR_fit image generated in this analysis
                    List_of_images.append(op.join(roots, f))
        List_of_images.sort()
        self.imglist = List_of_images
        self.n_img = len(List_of_images)
        
    def peak_fit(self, verbose=True, **kwargs):
        '''
        Run GDSC_SMLM peak fit on super-res tiff images and generate FitResults.txt
         
        '''

        d = deepcopy(kwargs)
        d['filelist'] = self.imglist
        d['results_dir'] = self.results_dir
        with open(self.json_file, 'w') as f:
            json.dump(d, f, indent=2)
        if d['GDSC_SMLM_version'] == 1:
            ImageJ_path = self.config_d['imageJ_GDSCSMLM1_path']
            GDSCSMLM_script_path = self.config_d['GDSCSMLM1_script_path']
        else:
            ImageJ_path = self.config_d['imageJ_GDSCSMLM2_path']
            GDSCSMLM_script_path = self.config_d['GDSCSMLM2_script_path']
        
        call([ImageJ_path, "--ij2", "--run", GDSCSMLM_script_path,
            'jsonpath="%s"' %self.json_file.replace('\\','//')])
        
        # Read the ImageJ log output
        with open(self.json_file, 'r') as f:
            d = json.load(f)
        assert 'fit_name' in d.keys(), 'Error in ImageJ script execution!\n'
        fit_name = d['fit_name']
        self._log_fit_info(fit_name, **kwargs)
        self.json_log['ImageJ_log'] = d['ImageJ_log']
         
        error_log = [(img, log) for img, log in d['ImageJ_log'].items() 
                    if 'Error' in ''.join(log.values())]
        if len(error_log) > 0:
            for i, l in error_log:
                for procedure, err in l.items():
                    warnings.warn("{}:\n\tError in ImageJ {} analysis, message: {}.".format(
                    i, procedure, err), UserWarning)
        
        return fit_name
        
    def _remove_fid(self, **kwargs):
        d = deepcopy(kwargs)
        s = d['fiducial_correction']['fid_size']
        b = d['fiducial_correction']['fid_brightness']
        gain = d['camera_gain']
        self._search_fitresults()
        List_of_corrected = []
        for fr in self.fitpathlist:
            file_dir = op.dirname(fr)
            if 'Fiducials.txt' in os.listdir(file_dir) and \
            'FitResults_FeuRemoved.txt' not in os.listdir(file_dir):
                fid_file = op.join(file_dir, 'Fiducials.txt')
                fit_res = op.join(file_dir, 'FitResults_Corrected.txt')
                List_of_corrected.append((fit_res, fid_file))
        
        for fit_res, fid_file in List_of_corrected:
            df = pd.read_table(fit_res)
            feus = pd.read_table(fid_file).values
            for feu in feus:
                feu = np.array(feu)
                feu = np.append(feu[0:2] - [s,s], feu[0:2] + [s,s])
                df = df[(df.origX<feu[0]) | (df.origY<feu[1]) | 
                 (df.origX>feu[2]) | (df.origY>feu[3])]
                if d['GDSC_SMLM_version'] == 1:
                    df = df[df.origValue<2*b]
                else:
                    df = df[df.origValue<(2*b/gain)]
            df = df.reset_index(drop=True)
            df.to_csv(op.join(op.dirname(fit_res), 
            'FitResults_FeuRemoved.txt'), index = False, sep = '\t')
        
        self.fit_name = 'FitResults_FeuRemoved.txt'
        
    def _create_symlinks(self):
        assert ctypes.windll.shell32.IsUserAnAdmin(), 'You need Admin'\
            ' permission to create symlink on Windows, try running Python as an'\
            ' administrator.\n\n'
        List_of_image_path_file = []
        for roots, dirs, files in os.walk(self.results_dir):
            for f in files:
                if f.endswith('.path.txt'):
                    List_of_image_path_file.append((roots, f))
        
        for r, f in List_of_image_path_file:
            with open(op.join(r, f)) as path_f:
                src = path_f.read()
            link = op.join(r, 'Raw_image.lnk')
            os.symlink(src, link)
            
    def _rendering(self, **kwargs):
        d = deepcopy(kwargs)
        if d['GDSC_SMLM_version'] == 1:
            ImageJ_path = self.config_d['imageJ_GDSCSMLM1_path']
            Rendering_script_path = self.config_d['Rendering_script_path1']
        elif d['GDSC_SMLM_version'] == 2:
            ImageJ_path = self.config_d['imageJ_GDSCSMLM2_path']
            Rendering_script_path = self.config_d['Rendering_script_path2']

        d2 = {'sr_scale': d['sr_scale'], 'results_dir': self.results_dir}
        with open(self.json_file, 'w') as f:
            json.dump(d2, f, indent=2)

        call([ImageJ_path, "--ij2", "--run", Rendering_script_path,
            'jsonpath="%s"' %self.json_file.replace('\\','//')])
                      
    def run_fit(self, image_condition=None, verbose=True, **kwargs):
        '''
        Run fit analysis on tiff image stacks acquired by super-resolution imaging. 
        Workflow all the analysis.
        Verbose set to False will supress all stdout logs but the start and finish
        times.
        
        '''

        # Check if the matched fitresults exists
        if kwargs['GDSC_SMLM_version'] not in [1,2]:
            raise NameError("GDSC_SMLM_version must be 1 or 2")

        if kwargs['GDSC_SMLM_peak_fit']:
            if verbose:
                sys.stdout.write('Looking for tiff image files meeting '\
                'defined image condition in \npath: {}...\n'.format(self.path))  
            self._search_images(image_condition)
            assert self.n_img > 0, 'No image found!\n\n'
            if verbose:
                sys.stdout.write("{} images found.\n\n".format(self.n_img))
    
            this_peak_fit_info = self._peak_fit_info(**kwargs)
            existing_peak_fit = self._existing_peak_fit()
            self.fitresults_folder = None
            for peak_fit_location, peak_fit_dict in existing_peak_fit.items():
                peak_fit_info = deepcopy(peak_fit_dict)
                del peak_fit_info['fit_name']
                if this_peak_fit_info == peak_fit_info:
                    self.fit_name = peak_fit_dict['fit_name']
                    self.fitresults_folder = peak_fit_location
                    break
            if not self.fitresults_folder:
                if verbose:
                    sys.stdout.write('No existing GDSC SMLM fit results found, '
                    'start running peak fit in ImageJ...\n\n') 
                self.fit_name = self.peak_fit(**kwargs)
                self.fitresults_folder = self.results_dir
            else:
                if verbose:
                    sys.stdout.write('Existing GDSC SMLM fit results found in path:\n'
                    '{}\n\n'.format(self.fitresults_folder)) 
        else:
            if verbose:
                sys.stdout.write('GDSC SMLM peak fit turned off.\n\n')
            self.fit_name = kwargs['fitresults_file_name']
            self.fitresults_folder = self.path
            
        self.cluster_fit(verbose=verbose, **kwargs)               
        self._logging(**kwargs)    
    
        if kwargs['GDSC_SMLM_peak_fit'] and kwargs['create_symlink_for_images']:
            self._create_symlinks()
        
        try: os.remove(self.json_file)
        except: pass
        
        sys.stdout.write("\nTime: {} Analysis finished!\n\n".format(
        datetime.now().strftime('%Y-%m-%d_%H-%M-%S')))

    def cluster_fit(self, verbose=True, **kwargs):           
        '''
        Run fit analysis on fitresults file generated by imageJ GDSC SMLM plugin.

        Accepted keywords in kwargs:
            pixel_size, DBSCAN_eps, DBSCAN_min_samples, sr_scale=8, frame_length=50,
            burst_filter={'run':True, 'min_burst':2, 'fill_gap':50},
            burst_analysis={'run':True, 'fill_gap':10, 'remove_single':True},
            length_measure=True, eccentricity_measure=True, convexhull_measure=True,
            save_histogram=True, save_GDSC_header_file = True
        '''
                    
        assert kwargs['GDSC_SMLM_peak_fit'] or kwargs['fitresults_file_name'] != 'default', \
        'Please enter fitresults_file_name'
        if kwargs['fitresults_file_name'] != 'default':
            self.fit_name = kwargs['fitresults_file_name']
        else:
            # Remove fiducials
            if kwargs['fiducial_correction']['run'] and kwargs['fiducial_correction']['fid_file'] != 'drift':
                self._remove_fid(**kwargs)
            elif kwargs['fiducial_correction']['run'] and kwargs['fiducial_correction']['fid_file'] == 'drift':
                self.fit_name = 'FitResults_Corrected.txt'
                
        self._search_fitresults()
        assert self.n_fit > 0, 'No fit results files found!\n'

        if verbose:
            sys.stdout.write("Looking for super-res fit result files '{}' in \npath: {}\n".format(
                self.fit_name, self.fitresults_folder))  
            
        self.resultslist = [SR_fit(fitpath) for fitpath in self.fitpathlist]
        if verbose:
            sys.stdout.write("{} fit results found.\n\n".format(self.n_fit))
        
        for fit in self.resultslist:
            new_dir = fit.root.replace(self.fitresults_folder, self.results_dir)
            if not op.isdir(new_dir):
                os.makedirs(new_dir)
                
            if kwargs['GDSC_SMLM_peak_fit']:
                for roots, dirs, files in os.walk(fit.root):
                    for f in files:
                        if f.endswith('.path.txt'):
                            f_root = roots
                            break
                    try:
                        copyfile(op.join(f_root, f), op.join(new_dir, f))
                    except SameFileError:
                        pass
    
            fit.input_parameters(results_dir=new_dir, **kwargs)
            fit_ID = 'T'+ self.timestring + '_' + fit.path # This serves as the
                        # primary key for the summary database
            fit.update_ID(fit_ID)
            fit.update_peakfitinfo()
            
        # Run temporal grouping
        paras = kwargs['temporal_grouping'].copy()
        if paras.pop('run', False):
            self._shinchoku_run(SR_fit.temporal_grouping, namae='temporal grouping',
            verbose=verbose, **paras)

        # Run cluster analysis
        paras = kwargs['cluster_analysis'].copy()
        if paras.pop('run', False):
            if paras['cluster_subject'] in ['mol', 'burst']:
                assert kwargs['temporal_grouping']['run'], '\nNo temporal grouping run, cannot set cluster analysis on mol or burst.'
            self._shinchoku_run(SR_fit.cluster_info, namae='cluster analysis',
            verbose=verbose, **paras)                

            paras = kwargs['length_measure'].copy()
            if paras.pop('run', False):
                self._shinchoku_run(SR_fit.length_measure, namae='length measurement',
                verbose=verbose, **paras)  
  
            if kwargs['eccentricity_measure']:
                self._shinchoku_run(SR_fit.eccentricity_measure, namae='eccentricity measurement',
                 verbose=verbose)                 
            if kwargs['convexhull_measure']:            
                self._shinchoku_run(SR_fit.convexhull_meausure, namae='convex hull area measurement',
                 verbose=verbose) 
               
        if kwargs['Rendering_SR']:
            if verbose:
                sys.stdout.write("\nRendering clustered SR images results...")
            for fit in self.resultslist:
                if kwargs['temporal_grouping']['run']:
                    fit.save_with_header(subj='burst', **kwargs)
                    fit.save_with_header(subj='mol', **kwargs)
                if kwargs['cluster_analysis']['run']: 
                    fit.save_with_header(subj='cluster', **kwargs)
                    fit.labelled_cluster()                    

            self._rendering(**kwargs)     
                       
        if verbose:
            sys.stdout.write("\nSaving results...")       
                 
        if kwargs['save_histogram']:
                self._save_histogram()        
        self._save_summary()
    
    def _shinchoku_run(self, func, namae='analysis', verbose=True, **paras):
        '''
        Loop over all images. Show progress if verbose.
        '''
        if verbose:
            sys.stdout.write("\nRun {}:\n".format(namae))             
            i = 1  
            for fit in self.resultslist:
                func(fit, **paras)
                sys.stdout.write('\r')
                sys.stdout.write("%d%%" % (i/self.n_fit*100))
                sys.stdout.flush()
                i += 1
        else:
            for fit in self.resultslist:
                func(fit, **paras)
        
    def _search_fitresults(self):
        '''
        This function finds the all the fit results file in fitresults_folder.

        '''
        List_of_fitresults = []
        for roots, dirs, files in os.walk(self.fitresults_folder):
            for f in files:
                if f == self.fit_name:
                    List_of_fitresults.append(op.join(roots, f))
        self.fitpathlist = List_of_fitresults
        self.n_fit = len(List_of_fitresults)
        
    def _logging(self, **kwargs):
        '''
        Write the analysis parameters into the log.json file
        
        '''
        d = deepcopy(kwargs)
        try: d['num_of_images'] = self.n_img
        except: d['num_of_images'] = np.nan
            
        try: d['num_of_fitresults'] = self.n_fit
        except: d['num_of_fitresults'] = np.nan
            
        try: d['fitresults_source'] = self.fitresults_folder
        except: d['fitresults_source'] = np.nan
            
        try: d['fitresults_name'] = self.fit_name
        except: d['fitresults_name'] = np.nan
            
        self.json_log.update(d)
        if not kwargs['GDSC_SMLM_peak_fit']:
            try: 
                with open(op.join(self.fitresults_folder, 'log.json'), 'r') as info:
                    peak_fit_log = json.load(info)
                self.json_log.update({k: peak_fit_log[k] for k in ['trim_track',
                 'signal_strength', 'precision', 'min_photons', 'fiducial_correction']})
            except:
                self.json_log.update({k: None for k in ['trim_track', 'signal_strength',
                'precision', 'min_photons', 'fiducial_correction']})
            
        with open(op.join(self.results_dir, "log.json"), 'w') as log:
            json.dump(self.json_log, log, indent=2)
        
    def _save_summary(self):
        '''
        Save the aggregated summary file
        '''
        to_save = []
        for fit in self.resultslist:
            to_save.append(fit.summarise())
            
        DF(to_save).to_csv(op.join(self.results_dir, 'Summary.csv'), index=False)
            
    def _save_histogram(self):
        for subj in ['burst', 'mol', 'cluster']:
            to_hist = []
            for fit in self.resultslist:
                try:
                    to_hist.append(fit.output_histogram(subj))
                    pd.concat(to_hist, ignore_index=True).to_csv(op.join(self.results_dir,
                     'Histogram_{}.csv'.format(subj)), index=False)
                except: pass
                                               
if __name__ == '__main__':
    import importlib
    import SR_fit_lib
    importlib.reload(SR_fit_lib)
    from SR_fit_lib import SR_fit
    
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
    'trim_track': {'run': True, 'frame_number': 4000, 'from_end': True}, # trim the stack to the required frame number, from_end controls trim from end or beginning 
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
     'tThresh':2500, 'min_frame':2, 'min_burst':2, 'min_on_prop':0},
     
    'cluster_analysis': {'run': True, 'cluster_subject':'burst', # burst, mol, loc
    'DBSCAN_eps_nm':200, 'DBSCAN_min_samples':2},
    'length_measure': {'run':False, 'algorithm':'close', 'sigma':2}, # algorithm: blur or close; if blur, sigma is the gaussian sigma; if close, sigma is the closing square size
    'eccentricity_measure': True ,
    'convexhull_measure': True , # measure the area of the cluster
    'Rendering_SR': True , # rendering clustered SR images in imageJ
    'save_histogram':True,
    }

    Analysis(directory).run_fit(verbose=True, image_condition=image_condition, **input_dict)