'''
SR_toolkit_v3_3

Issued: 01-02-2019
@author: Eric Kobayashi

Changes over v3_2:
    
    1. Give the option of running the fiducial correction based on given fiducial coordinates.
    
    2. Add fiducial correction parameters (fiducial last time, smoothing parameter, limit_smoothing).
    
    3. Add camera bias, gain in the settings.
    
    4. Add labelled images for clusters.
    
    5. Add background calculation and variable background correction.
    
    6. Rendering SR images after cluster analysis.
    
    7. Add closing algorithm to calculate length. Instead of blurring the skeletion,
    the holes between the localisations are filled so that no loops will appear in the skeleton. 
    
    8. Fix some bugs that cause error_log to malfunction.
    
Notice:
    
    Please put Rendering_SR.py in the same folder as GDSC_SMLM_peak_fit.py
    
    
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
import SR_fit_lib



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
        self.json_log = OrderedDict()
        self.json_log['Time'] = timestring
        self.json_log['Path'] = path
        
    def _peak_fit_info(self, **kwargs):
        '''
        Return all the parameters related to the ImageJ peak fit and fiducial
        correction.
        '''
        
        peak_fit_para = {k: kwargs[k] for k in ['pixel_size', 'sr_scale', 
        'frame_length', 'camera_bias', 'camera_gain', 'BG_measurement', 'trim_track', 
        'signal_strength', 'precision', 'min_photons', 'fiducial_correction']}
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
        for json_root, json_file in List_of_peak_fit_info_file:
            with open(op.join(json_root, json_file), 'r') as f:
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
        json_file = d['json_file']
        with open(json_file, 'w') as f:
            json.dump(d, f, indent=2)

        ImageJ_path = d['ImageJ_path']
        GDSCSMLM_script_path = d['GDSCSMLM_script_path']
        if not GDSCSMLM_script_path.endswith('.py'):
            GDSCSMLM_script_path = op.join(GDSCSMLM_script_path, 'Peakfit_GDSC_SMLM.py')
        call([ImageJ_path, "--ij2", "--run", GDSCSMLM_script_path,
            'jsonpath="%s"' %json_file])
        
        # Read the ImageJ log output
        with open(json_file, 'r') as f:
            d = json.load(f)
        assert 'fit_name' in d.keys(), 'Error in ImageJ script execution!\n'
        fit_name = d['fit_name']
        self._log_fit_info(fit_name, **kwargs)
        self.json_log['ImageJ_log'] =  d['ImageJ_log']
         
        error_log = [(img, log) for img, log in d['ImageJ_log'].items() 
                    if 'Error' in ''.join(log.values())]
        if len(error_log) > 0:
            for i, l in error_log:
                for procedure, err in l.items():
                    warnings.warn("{}:\n\tError in ImageJ {} analysis, message: {}.".format(
                    i, procedure, err), UserWarning)
        
        return fit_name
        
    def _remove_fid(self, b, s):
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
                df = df[(df.X<feu[0]) | (df.Y<feu[1]) | (df.X>feu[2]) | (df.Y>feu[3])]
                df = df[df.origValue<2*b]
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
        ImageJ_path = d['ImageJ_path']
        Rendering_script_path = d['GDSCSMLM_script_path']
        json_file = d['json_file']
        with open(json_file, 'r') as f:
            d2 = json.load(f)
        d2['results_dir'] = self.results_dir
        with open(json_file, 'w') as f:
            json.dump(d2, f, indent=2)

        if not Rendering_script_path.endswith('.py'):
            Rendering_script_path = op.join(Rendering_script_path, 'Rendering_SR.py')
        else:
            Rendering_script_path = op.join(op.dirname(Rendering_script_path), 'Rendering_SR.py')
            
        call([ImageJ_path, "--ij2", "--run", Rendering_script_path,
            'jsonpath="%s"' %json_file])
                      
    def run_fit(self, image_condition=None, verbose=True, **kwargs):
        '''
        Run fit analysis on tiff image stacks acquired by super-resolution imaging. 
        Workflow all the analysis.
        Verbose set to False will supress all stdout logs but the start and finish
        times.
        
        '''

        # Check if the matched fitresults exists
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
            self.fit_name = 'NotApplicable'
            self.fitresults_folder = self.path
            
        self.cluster_fit(verbose=verbose, **kwargs)               
        self._logging(**kwargs)    
    
        if kwargs['GDSC_SMLM_peak_fit'] and kwargs['create_symlink_for_images']:
            self._create_symlinks()
        
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
            if kwargs['fiducial_correction']['run']:
                self._remove_fid(kwargs['fiducial_correction']['fid_brightness'],
                kwargs['fiducial_correction']['fid_size'])
                
        self._search_fitresults()
        assert self.n_fit > 0, 'No fit results files found!\n'

        if verbose:
            sys.stdout.write("Looking for super-res fit result files '{}' in \npath: {}\n".format(
                self.fit_name, self.fitresults_folder))  
            
        self.resultslist = [SR_fit_lib.SR_fit(fitpath) for fitpath in self.fitpathlist]
        # OrderedDefaultDict is used here. (defined below)
        self.to_summary = OrderedDefaultDict(list) 
        self.to_hist_cluster = OrderedDefaultDict(list) 
        self.to_hist_len = OrderedDefaultDict(list) 
        self.to_hist_darklen = OrderedDefaultDict(list) 
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
            self.to_summary['Analysis_ID'].append(fit.fit_ID)
            self.to_summary['Raw_loc_number'].append(fit.loc_num())
            self.to_summary['Frame_number'].append(fit.framenum)
            if kwargs['BG_measurement']['run']:
                self.to_summary['BG_level'].append(
                 fit.output_attr_fromfile('median_intensity.txt'))
                if kwargs['BG_measurement']['correct_precision']:
                    self.to_summary['Corrected_precision'].append(
                     fit.output_attr_fromfile('corrected_precision.txt'))
            if 'Fiducials.txt' in os.listdir(fit.root):
                fid_file = op.join(fit.root, 'Fiducials.txt')
                feus = pd.read_table(fid_file).values
                self.to_summary['Fiducials_number'].append(len(feus))
            else:
                self.to_summary['Fiducials_number'].append(0)
            
        if kwargs['cluster_analysis_measurement']:
            # Run cluster analysis
            if verbose:
                sys.stdout.write("Run cluster analysis:\n")        
                i = 1
            for fit in self.resultslist:
                fit.cluster_info()
                
                if verbose:
                    sys.stdout.write('\r')
                    sys.stdout.write("%d%%" % (i/self.n_fit*100))
                    sys.stdout.flush()
                    i += 1
                    
            if verbose:
                sys.stdout.write("\nRun burst filter:\n")             
                i = 1
            # Run burst filter if applied
            paras = kwargs['burst_filter'].copy()
            if paras.pop('run', False):
                for fit in self.resultslist:
                    fit.burst_filter(**paras)
                    self.to_summary['Filtered_cluster'].append(fit.filtered_cluster)
                    
                    if verbose:
                        sys.stdout.write('\r')
                        sys.stdout.write("%d%%" % (i/self.n_fit*100))
                        sys.stdout.flush()
                        i += 1
                    
            # Output cluster info after burst filter        
            for fit in self.resultslist:
                self.to_summary['Clustered_loc_number'].append(fit.loc_num())
                self.to_summary['Number_of_clusters'].append(fit.cluster_num())
                self.to_summary['Average_loc_per_cluster'].append(fit.ave_cluster_locnum())
                
                if kwargs['save_histogram']:
                    self.to_hist_cluster['Analysis_ID'] += [fit.fit_ID]*fit.cluster_num()
                    self.to_hist_cluster['Cluster_ID'] += fit.output_histogram('num')
                    self.to_hist_cluster['Loc_number'] += fit.output_histogram('loc_num')
                    
            # Run burst analysis if applied    
            if verbose:
                sys.stdout.write("\nRun burst analysis:\n")             
                i = 1  
            
            paras = kwargs['burst_analysis'].copy()
            if paras.pop('run', False):
                for fit in self.resultslist:
                    fit.burst_info(**paras)
    
                    self.to_summary['Average_burst_number'].append(fit.ave_burstnum)
                    self.to_summary['Average_burst_length_(s)'].append(fit.ave_burstlen)
                    self.to_summary['Average_dark_length_(s)'].append(fit.ave_darklen)
                    self.to_summary['Average_lighttodark_ratio'].append(fit.ave_lighttodark)
                    
                    if kwargs['save_histogram']:
                        self.to_hist_cluster['Burst_number'] += fit.output_histogram('burstnum')
                        self.to_hist_cluster['Average_burst_length_(s)'] += fit.output_histogram('ave_burstlen')
                        self.to_hist_cluster['Average_dark_length_(s)'] += fit.output_histogram('ave_darklen')
                        
                        alllen = fit.output_noncluster_histogram('alllen')
                        self.to_hist_len['Analysis_ID'] += [fit.fit_ID]*len(alllen[0])
                        self.to_hist_len['Cluster_ID'] += alllen[0]
                        self.to_hist_len['All_burst_length_(s)'] += alllen[1]
                        
                        alldarklen = fit.output_noncluster_histogram('alldarklen')
                        self.to_hist_darklen['Analysis_ID'] += [fit_ID]*len(alldarklen[0])
                        self.to_hist_darklen['Cluster_ID'] += alldarklen[0]
                        self.to_hist_darklen['All_dark_length_(s)'] += alldarklen[1]
                    
                    if verbose:
                        sys.stdout.write('\r')
                        sys.stdout.write("%d%%" % (i/self.n_fit*100))
                        sys.stdout.flush()
                        i += 1
    
            if verbose:
                sys.stdout.write("\nRun length measurement:\n")             
                i = 1     
            paras = kwargs['length_measure'].copy()       
            if paras.pop('run', False):
                for fit in self.resultslist:
                    fit.length_measure(**paras)
                    self.to_summary['Average_length_(nm)'].append(fit.ave_length)
                    self.to_summary['Average_density_1D'].append(fit.ave_density_1D)
                    
                    if kwargs['save_histogram']:
                        self.to_hist_cluster['Length_(nm)'] += fit.output_histogram('nm_length')
                        
                    if verbose:
                        sys.stdout.write('\r')
                        sys.stdout.write("%d%%" % (i/self.n_fit*100))
                        sys.stdout.flush()
                        i += 1
                        
            if verbose:
                sys.stdout.write("\nRun eccentricity measurement:\n")             
                i = 1      
            if kwargs['eccentricity_measure']:
                for fit in self.resultslist:
                    fit.eccentricity_measure()
                    self.to_summary['Average_eccentricity'].append(fit.ave_ecc)
                    self.to_summary['Average_flattening'].append(fit.ave_flattening)
                    
                    if kwargs['save_histogram']:
                        self.to_hist_cluster['Eccentricity'] += fit.output_histogram('ecc')
                        self.to_hist_cluster['Flattening'] += fit.output_histogram('flattening')
                        
                    if verbose:
                        sys.stdout.write('\r')
                        sys.stdout.write("%d%%" % (i/self.n_fit*100))
                        sys.stdout.flush()
                        i += 1
    
            if verbose:
                sys.stdout.write("\nRun convex hull area measurement:\n")             
                i = 1           
            if kwargs['convexhull_measure']:
                for fit in self.resultslist:
                    fit.convexhull_meausure()
                    self.to_summary['Average_convexhull_area_(nm2)'].append(fit.ave_area)
                    self.to_summary['Average_density_2D'].append(fit.ave_density_2D)
                    
                    if kwargs['save_histogram']:
                        self.to_hist_cluster['Convex_hull_area_(nm2)'] += fit.output_histogram('nm2_area')
                        
                    if verbose:
                        sys.stdout.write('\r')
                        sys.stdout.write("%d%%" % (i/self.n_fit*100))
                        sys.stdout.flush()
                        i += 1

            for fit in self.resultslist:
                fit.save_with_header()
                
            if kwargs['Rendering_SR']:
                if verbose:
                    sys.stdout.write("\nRendering clustered SR images results...")    
                for fit in self.resultslist:
                    fit.labelled_cluster()
                self._rendering(**kwargs)
                 
            if verbose:
                sys.stdout.write("\nSaving results...")       
            if kwargs['save_histogram']:
                df = self._save_histogram()
                if kwargs['Rendering_SR']:
                    for fit in self.resultslist:
                        all_clusters_info = op.join(fit.results_root, 'all_clusters_info.csv')
                        try:
                            fit_hist = df[df['Analysis_ID'] == fit.fit_ID].copy()
                        except TypeError:
                            continue
                        fit_hist.drop(columns=['Analysis_ID']).to_csv(all_clusters_info, index=False)
            
        else:
            if verbose:
                sys.stdout.write("\nNo cluster analysis run. Saving results...")         
        self._save_summary()
                

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
        try:
            d['num_of_images'] = self.n_img
        except:
            d['num_of_images'] = np.nan
            
        try:
            d['num_of_fitresults'] = self.n_fit
        except:
            d['num_of_fitresults'] = np.nan
            
        try:
            d['fitresults_source'] = self.fitresults_folder
        except:
            d['fitresults_source'] = np.nan
            
        try:
            d['fitresults_name'] = self.fit_name
        except:
            d['fitresults_name'] = np.nan
            
        self.json_log.update(d)
        if not kwargs['GDSC_SMLM_peak_fit']:
            self.json_log['trim_track'] = 'NotApplicable'
            self.json_log['signal_strength'] = 'NotApplicable'
            self.json_log['precision'] = 'NotApplicable'
            self.json_log['min_photons'] = 'NotApplicable'
            self.json_log['fiducial_correction'] = 'NotApplicable'
            
        with open(op.join(self.results_dir, "log.json"), 'w') as log:
            json.dump(self.json_log, log, indent=2)
        
    def _save_summary(self):
        DF(self.to_summary).to_csv(op.join(self.results_dir, 'Summary.csv'),
            columns=self.to_summary.keys(), index=False)
            
    def _save_histogram(self):
        df_cluster = DF(self.to_hist_cluster)
        df_cluster.to_csv(op.join(self.results_dir, 'Histogram_clusters.csv'),
            columns=self.to_hist_cluster.keys(), index=False)
            
        DF(self.to_hist_len).to_csv(op.join(self.results_dir, 'Histogram_all_burstlen.csv'),
            columns=self.to_hist_len.keys(), index=False)
            
        DF(self.to_hist_darklen).to_csv(op.join(self.results_dir, 'Histogram_all_darklen.csv'),
            columns=self.to_hist_darklen.keys(), index=False)
            
        return df_cluster

class OrderedDefaultDict(OrderedDict):
    '''
    Source: http://stackoverflow.com/a/6190500/562769
    Order is preserved so the output dataframe columns follow the set order. 
    Default method is used for dynamically adding new fields.
    
    '''
    def __init__(self, default_factory=None, *a, **kw):
        if (default_factory is not None and
           not hasattr(default_factory, '__call__')):
            raise TypeError('first argument must be callable')
        OrderedDict.__init__(self, *a, **kw)
        self.default_factory = default_factory

    def __getitem__(self, key):
        try:
            return OrderedDict.__getitem__(self, key)
        except KeyError:
            return self.__missing__(key)

    def __missing__(self, key):
        if self.default_factory is None:
            raise KeyError(key)
        self[key] = value = self.default_factory()
        return value

    def __reduce__(self):
        if self.default_factory is None:
            args = tuple()
        else:
            args = self.default_factory,
        return type(self), args, None, None, self.items()

    def copy(self):
        return self.__copy__()

    def __copy__(self):
        return type(self)(self.default_factory, self)

    def __deepcopy__(self, memo):
        import copy
        return type(self)(self.default_factory,
                          copy.deepcopy(self.items()))

    def __repr__(self):
        return 'OrderedDefaultDict(%s, %s)' % (self.default_factory,
                                               OrderedDict.__repr__(self))
                                               
if __name__ == '__main__':
    import importlib
    importlib.reload(SR_fit_lib)
    
    directory = r"C:\Users\Eric\OneDrive - University Of Cambridge\igorplay\feudicial"
    image_condition = lambda img: img.endswith('_561.tif')  # Change conditions to search for imagestacks
    input_dict = {
    # ==== Parameters for all analysis ====
    'pixel_size': 98.6 , # nm
    'sr_scale': 8 , # The scale used in length analysis and generation of super-res images 
    'camera_bias': 500 ,
    'camera_gain': 55.50 ,
    'frame_length': 50 , # ms per frame, i.e. exposure time
    'create_symlink_for_images': False, # Needs admin mode
    
    # ==== Parameters for GDSC SMLM fitting ====
    'GDSC_SMLM_peak_fit': True, # If False, the parameters for GDSC SMLM fitting will be ignored, only cluster_analysis will be run
    'trim_track': {'run': False, 'frame_number': 4000, 'from_end': True}, # trim the stack to the required frame number, from_end controls trim from end or beginning 
    'BG_measurement': {'run': True, 'correct_precision': True, 'kappa': 2}, # measure the background of the image_stack. (median pixel value) 
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