'''
The library includes the classes and methods related to DBSCAN cluster analysis
and subsequent measurements.

'''

import os.path as op
import pandas as pd
import numpy as np
import warnings
import os
import scipy
import json
from math import ceil
from pandas import Series
from pandas import DataFrame as DF
from imageio import imwrite
from skimage.morphology import skeletonize_3d, closing, square
from scipy.ndimage.filters import gaussian_filter
from scipy.spatial import ConvexHull
from skimage.measure import regionprops
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN

class SR_fit(object):
    '''
    This class is to store all the information and methods on one fitresults 
    file. 
    '''
        
    def __init__(self, filepath):
        '''
        Initialise the fit file into a dataframe.
        
        '''    
            
        self.path = filepath
        self.root = op.dirname(filepath)
        
        # For compatibility with GDSC SMLM 2
        with open(filepath, 'r') as f:
            contents = f.read()
        with open(filepath, 'w') as f:
            f.write(contents.replace('*',''))

        self.to_summary = {}         # initialise output

        try:
            self.df = pd.read_csv(filepath, sep='\t')       # default delimitor: space
            if len(self.df.columns) <= 2:    # if it is a comma delimitor file
                self.df = pd.read_csv(filepath)
        except:
            self._isfit = False          # Wrong filetype
            return

        self.df = self.df.astype('float64', errors='ignore')
        self._isfit = True
        if len(self.df) == 0:
            self._empty = True
            return
        self._empty = False
        
        # For compatibility with GDSC SMLM 2
        self.df = self.df.rename(columns={'T':'Frame', 'X (px)':'X', 'Y (px)':'Y'})
        
        # Shift negative values
        shift = (min(self.df['X'].min(), 0), min(self.df['Y'].min(), 0))
        self.df['X'] = self.df['X'] - shift[0]
        self.df['Y'] = self.df['Y'] - shift[1]
        self.to_summary['xy_shift'] = shift

        self.overhead = 200 # For clusters that are more than 200 subpixels out of edge, remove them
            
    def input_parameters(self, results_dir=None, pixel_size=None, sr_scale=8, 
     frame_length=50, GDSC_SMLM_version=1, **kwargs):
        '''
        Initialise all the analysis parameters.
        This method is essential.
        
        '''
        
        self.pixel_size = pixel_size
        self.sr_scale = sr_scale
        self.frame_length = frame_length/1000 # convert to sec
        self.results_root = results_dir
        self.version = GDSC_SMLM_version
        
        # Extract image information from image.results.xls
        image_details_file = "Image.results.xls"
        image_details = op.join(self.root, image_details_file)
        assert op.isfile(image_details), "Image metadata file not found!"
        with open(image_details, 'r') as metadata:
            metadata_content = metadata.read()
            self.framenum = int(SR_fit._parse_xml(metadata_content, 'frames'))
            self.width = int(SR_fit._parse_xml(metadata_content, 'width'))
            self.height = int(SR_fit._parse_xml(metadata_content, 'height'))
            self.header = metadata_content.split('#Frame')[0]
                                
    def cluster_info(self, cluster_subject='loc', DBSCAN_eps_nm=100, 
     DBSCAN_min_samples=5):
        '''
        Extract information from the cluster analysis results.        
        '''
        
        if self._isnan():
            # No localisations
            return
        
        cluster_analysis_file = "DBSCAN_eps_{:.2f}nm_min_{}.csv".format(
                                DBSCAN_eps_nm, DBSCAN_min_samples)
        cluster_analysis_file_path = op.join(self.results_root, 
                                    cluster_analysis_file)

        self._cluster_analysis(cluster_subject=cluster_subject, 
         DBSCAN_eps_nm=DBSCAN_eps_nm, DBSCAN_min_samples=DBSCAN_min_samples)
         
        # Save cluster analysis result
        self.df.to_csv(cluster_analysis_file_path, index=False)
        
        clustered_locs = self.df[self.df.Cluster > 0].copy()

        # Remove overflown clusters due to possible fiducial correction error
        dim_limitx = self.width + int(self.overhead/self.sr_scale)
        dim_limity = self.height + int(self.overhead/self.sr_scale)
        if len(clustered_locs[(clustered_locs['X'] >=dim_limitx) | (clustered_locs['Y'] >= dim_limity)]) > 0:
            print('Index overflow, fiducial corrections may be incoorect, please check:\n{}'.format(self.path))
            clustered_locs = clustered_locs[(clustered_locs['X'] < dim_limitx) & (clustered_locs['Y'] < dim_limity)].copy()

        self.unclustered_df = self.df[~self.df.isin(clustered_locs)].copy()
        self.df = clustered_locs
        if len(clustered_locs) == 0:
            # No clustered localisations
            self._clustered = False
            return
        else:
            self._clustered = True
            
        self.clusterlist = []
        self.all_cluster_info = {}
        cluster_labels = sorted(clustered_locs.Cluster.unique())
        clu_num = 1
        for clu_labels in cluster_labels:
            self.clusterlist.append(cluster_track(clu_num, cluster_subject,
            self.frame_length, clustered_locs[clustered_locs.Cluster == clu_labels]))
            clu_num += 1  # Cluster analysis could miss some cluster number label
            
        # Save for summary
        self.mean_cluattr('num', 'Cluster_ID')
        try: self.to_summary['Average_cluster_locnum'] = self.mean_cluattr(
             'loc_num', 'Num_localisation')
        except: pass
        
        try: self.to_summary['Average_cluster_burstnum'] = self.mean_cluattr(
             'burst_num', 'Num_burst')
        except: pass
            
        try: self.to_summary['Average_cluster_molnum'] = self.mean_cluattr(
             'mol_num', 'Num_molecule')
        except: pass
        
        self.to_summary['Number_of_cluster'] = clu_num - 1
      
    def length_measure(self, algorithm='blur', sigma=2):
        if self._isnan():
            # No localisations
            self.ave_length, self.ave_density_1D = (np.nan,)*2
            return
            
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=FutureWarning)
            self._skeletonise(algorithm, sigma)
            for clu in self.clusterlist:
                clu.length_measure(self.pixel_size/self.sr_scale)
        
        self.to_summary['Average_cluster_length'] = self.mean_cluattr(
         'nm_length', 'Length')    
                           
    def eccentricity_measure(self):
        if self._isnan():
            # No localisations
            self.ave_ecc, self.ave_flattening = (np.nan,)*2
            return

        self._regionprops()
        for clu in self.clusterlist:
            setattr(clu, 'ecc', self.props[clu.num-1]['eccentricity'])
            setattr(clu, 'flattening', 1 - np.sqrt(1 - clu.ecc**2))
        self.to_summary['Average_cluster_eccentricity'] = self.mean_cluattr(
         'ecc', 'Eccentricity')
        self.to_summary['Average_cluster_flattening'] = self.mean_cluattr(
         'flattening', 'Flattening')
   
    def convexhull_meausure(self):    
        '''
        Calculate the area of all clusters' convex hull
        This analysis replaces the nearest neighbour analysis in ver.2.
        This method is optional.
        
        '''
        
        if self._isnan():
            # No localisations
            self.ave_area, self.ave_density_2D = (np.nan,)*2
            return
            
        self._regionprops()
        scale = self.pixel_size/self.sr_scale
        for clu in self.clusterlist:
            setattr(clu, 'convex_area', self.props[clu.num-1]['convex_area']*(scale**2))
        self.to_summary['Average_cluster_area'] = self.mean_cluattr(
         'convex_area', 'Convex_area')  
            
    def update_ID(self, keystring):
        self.fit_ID = keystring
        self.to_summary['Analysis_ID'] = keystring
    
    def update_peakfitinfo(self):
        '''
        Save the information of the images gained from peak fit.
        '''
        if self._isnan():
            self.to_summary['Raw_loc_number'] = 0
        else:
            self.to_summary['Raw_loc_number'] = len(self.df)
        self.to_summary['Frame_number'] = self.framenum
        try:
            self.to_summary['BG_level'].append(
             self.output_attr_fromfile('median_intensity.txt'))
        except:
            pass
        
        try:
            self.to_summary['Corrected_precision'].append(
             self.output_attr_fromfile('corrected_precision.txt'))
        except:
            pass
            
        if 'Fiducials.txt' in os.listdir(self.root):
            fid_file = op.join(self.root, 'Fiducials.txt')
            feus = pd.read_table(fid_file).values
            self.to_summary['Fiducials_number'] = len(feus)
        else:
            self.to_summary['Fiducials_number'] = 0
               
    def labelled_cluster(self):
        if self._isnan():
            # No localisations
            return

        roi = DF([np.array(list(map(ceil, (clu.bbox()*self.sr_scale)))) 
            for clu in self.clusterlist])
        roi_file = op.join(self.results_root, 'clusters_roi.txt')
        roi.to_csv(roi_file, index=False)
        
    def temporal_grouping(self, dThresh=20, min_loc=2, max_mol_area=float('inf'),
     tThresh=2500, min_frame=2, min_burst=1, min_on_prop=0):
        '''
        New in version 3.5. Conduct a spatial and temporal grouping to identify
        individual binding events.
        '''
        if self._isnan():
            # No localisations
            return

        df_ = self.df.copy()
        fit_xy = df_[['X', 'Y']].values.astype('float64')
        
        SUBPIXEL = self.pixel_size/self.sr_scale
        tThresh_frame = tThresh/self.frame_length/1000
        dThresh_pixel = dThresh/self.pixel_size
        
        # spatial grouping at molecular scale (dThresh ~ precision)
        db = DBSCAN(eps=dThresh_pixel, min_samples=min_loc).fit(fit_xy) 
        
        df_['Molecule_ID'] = db.labels_ + 1
        clu_num = max(db.labels_) + 1
        burst_df = []
        track_save = []
        
        for i in range(1, clu_num+1):
            track = df_[df_['Molecule_ID'] == i].copy()
            try:
                mol_area = ConvexHull(track[['X','Y']].values).area*(self.pixel_size**2)
            except scipy.spatial.qhull.QhullError:
                mol_area = SUBPIXEL**2 # subpixel area
            if mol_area > max_mol_area: # max_mol_area filters out overlapping molecules (inspired by Dave!)
                continue
                
            fit_t = track['Frame'].values.reshape(-1, 1)
            dbt = DBSCAN(eps=tThresh_frame, min_samples=min_frame).fit(fit_t)
            burst_num = max(dbt.labels_) + 1
            if burst_num < min_burst: # min_burst filters out orphan events
                continue
            track['Burst_ID'] = dbt.labels_ + 1
            start = 1
            
            for b, charas in track.groupby('Burst_ID'):
                if b == 0: continue # Non-clustered events
                burst = charas.mean()
                burst['ON_time'] = len(charas)*self.frame_length # The actual ON time (excluding blinking time)
                burst['ON_span'] = (max(charas['Frame'])-min(charas['Frame'])+1)*self.frame_length # The apparent ON time (including blinking time)
                burst['ON_prop'] = burst['ON_time']/burst['ON_span']
                if burst['ON_prop'] < min_on_prop: # min_on_prop filters out 'fake' bursts that simply join orphans
                    continue
                burst['OFF_time'] = (min(charas['Frame'])-start)*self.frame_length
                if burst['OFF_time'] == 0: burst['OFF_time'] = np.nan # The dark time between last and this blink
                burst['Area'] = mol_area
                start = max(charas['Frame']) + 1
                burst_df.append(burst)
                
            track_save.append(track)
        try:
            pd.concat(track_save, ignore_index=True).to_csv(op.join(
            self.results_root, 'temporal_grouping.csv'), index=False)
        except ValueError:
            self._empty = True
            return
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=np.ComplexWarning)
            self.burst_df = DF(burst_df).astype({'Frame':int, 'origX':int, 'origY':int, 
                'origValue':float, 'Error':float, 'Noise':float, 'SNR':float,  
                'Background':float,  'Signal':float,  'Angle':float, 'X':float,   
                'Y':float,   'X SD':float, 'Y SD':float,'Precision (nm)':float,  
                'Molecule_ID':int, 'Burst_ID':int,'ON_time':float, 
                'ON_span':float, 'ON_prop':float, 'OFF_time':float,'Area':float}, errors='ignore')
        
        mol_df = []
        for m, charas in self.burst_df.groupby('Molecule_ID'):
            mol = charas.mean(skipna=True)
            mol['ON_time'] = charas['ON_time'].sum()
            mol['ON_span'] = charas['ON_span'].sum()
            mol['OFF_time'] = charas['OFF_time'].sum(skipna=True)
            mol['Burst_number'] = len(charas)
            mol_df.append(mol)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=np.ComplexWarning)
            self.mol_df = DF(mol_df).astype({'Frame':int, 'origX':int, 'origY':int, 
                'origValue':float, 'Error':float, 'Noise':float, 'SNR':float,  
                'Background':float,  'Signal':float,  'Angle':float, 'X':float,   
                'Y':float,   'X SD':float, 'Y SD':float,'Precision (nm)':float,  
                'Molecule_ID':int, 'Burst_ID':int,'ON_time':float, 
                'ON_span':float, 'ON_prop':float, 'OFF_time':float,'Area':float,
                'Burst_number':int}, errors='ignore')
        
        # Save results in self.to_summary
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self.to_summary['Total_burstnum'] = len(self.burst_df)
            self.to_summary['Total_molnum'] = len(self.mol_df)
            self.to_summary['Average_burst_time'] = np.nanmean(self.burst_df['ON_time'])
            self.to_summary['Average_burst_span'] = np.nanmean(self.burst_df['ON_span'])
            self.to_summary['Average_burst_prop'] = np.nanmean(self.burst_df['ON_prop'])
            self.to_summary['Average_burst_darktime'] = np.nanmean(self.burst_df['OFF_time'])
            self.to_summary['Average_burst_precision'] = np.nanmean(self.burst_df['Precision (nm)'])
            self.to_summary['Average_mol_burstnum'] = np.nanmean(self.mol_df['Burst_number'])
            self.to_summary['Average_mol_area'] = np.nanmean(self.mol_df['Area'])
                   
    # Output methods
    def save_with_header(self, subj=None, **kwargs):
        if self._isnan():
            # No localisations
            return
            
        if subj == 'burst':
            try: df = self.burst_df.copy()
            except: return
        elif subj == 'mol':
            try: df = self.mol_df.copy()
            except: return
        elif subj == 'cluster':
            try: df = self.df.copy()
            except: return
        else:
            return 

        hf = op.join(self.results_root, 'All_{}_header.txt'.format(subj))
        to_output = df.drop(columns=['Source', 'Molecule_ID', 'Burst_ID', 
         'ON_time', 'ON_span', 'ON_prop', 'OFF_time', 'Burst_number', 'Area',
         'Cluster', 'SNR', 'Precision (nm)'], errors='ignore')
        if self.version == 2:
            to_output = to_output.rename(columns={'X': 'X (px)', 'Y':'Y (px)'})
            to_output['Z (px)'] = 0

        to_output.to_csv(hf, sep = '\t',  index = False)
        with open(hf, 'r+') as log:
            content = log.read()
            log.seek(0)
            log.write(self.header + '#' + content)
        
    def cluster_num(self):
        if hasattr(self, 'clusterlist'):
            return len(self.clusterlist)
        else:
            return 0
                
    def output_histogram(self, subj):
        '''
        Output the histogram for subj, either bursts, molecules or clusters
        add Analysis_ID and save in the result folder
        '''

        if subj == 'burst':
            try: # Cluster_ID exists
                to_save = self.burst_df[['Molecule_ID', 'Burst_ID', 'Cluster_ID', 
                'ON_time', 'ON_span', 'ON_prop', 'OFF_time', 'Area', 'Frame', 
                'origValue', 'Precision (nm)']].copy()
            except:
                to_save = self.burst_df[['Molecule_ID', 'Burst_ID', 
                'ON_time', 'ON_span', 'ON_prop', 'OFF_time', 'Area', 'Frame', 
                'origValue', 'Precision (nm)']].copy()
        elif subj == 'mol':
            try: # Cluster_ID exists
                to_save = self.mol_df[['Molecule_ID', 'Burst_ID', 'Cluster_ID', 
                'ON_time', 'ON_span', 'ON_prop', 'OFF_time', 'Area', 'Frame', 
                'origValue', 'Precision (nm)']].copy()
            except:
                to_save = self.mol_df[['Molecule_ID', 'Burst_ID', 
                'ON_time', 'ON_span', 'ON_prop', 'OFF_time', 'Area', 'Frame', 
                'origValue', 'Precision (nm)']].copy()
        elif subj == 'cluster':
            to_save = DF(self.all_cluster_info)
        else:
            return None
        to_save['Analysis_ID'] = self.fit_ID
        to_save.to_csv(op.join(self.results_root, 'All_{}_info.csv'.format(
         subj)), index=False) 
            
        return to_save
        
    def mean_cluattr(self, attr, title):
        '''
        Calculate the mean value of the attribute and save it as histograms
        '''
        all_attr = [getattr(clu, attr) for clu in self.clusterlist]
        self.all_cluster_info[title] = Series(all_attr)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            return np.nanmean(all_attr)  
        
    # Private methods
    def _isnan(self):
        if not hasattr(self, '_clustered'):
            if not self._isfit or self._empty:
                return True
            else:
                return False
        else:
            if not self._isfit or self._empty or not hasattr(self, 'clusterlist'):
                return True
            elif not self._clustered or len(self.clusterlist) == 0:
                return True
            else:
                return False
            
    def _cluster_analysis(self, cluster_subject='loc', DBSCAN_eps_nm=100, 
         DBSCAN_min_samples=5):
        '''
        Do DBSCAN cluster analysis on the dataframe.
        This method is only called inside the class because it is essential 
        to all the following analysis.
        
        '''
        if cluster_subject == 'mol': # for PAINT
            self.df = self.mol_df
        elif cluster_subject == 'burst': # for STORM/PALM
            self.df = self.burst_df 
        else: # default: localisation
            pass # self.df contains all localisations

        fit_xy = self.df[['X', 'Y']].values

        db = DBSCAN(eps=(DBSCAN_eps_nm/self.pixel_size), min_samples=DBSCAN_min_samples).fit(fit_xy)
        self.df['Cluster'] = db.labels_ + 1  # Start labelling with 1, makes analysis easier

        # Feed the cluster_ID info back into the cluster_subject df
        if cluster_subject == 'mol': 
            self.mol_df['Cluster_ID'] = db.labels_ + 1
        elif cluster_subject == 'burst': 
            self.burst_df['Cluster_ID'] = db.labels_ + 1

    def _cluster_dict(self):
        '''
        Create a cluster number -> cluster_track object dictionary.
        
        '''
        if self._isnan():
            # No localisations
            return {}
            
        return {clu.num: clu for clu in self.clusterlist}
                
    def _regionprops(self):
        '''
        run the skimage.measure.regionprops function to measure the eccentricity
        and convex hull area of the aggregates in the labelled image.
        
        '''
        
        if not hasattr(self, 'props'):
            self.props = regionprops(self._labelled_image(), coordinates='xy')
            
    def _skeletonise(self, algorithm, sigma):
        '''
        Skeletonise the fibrils to one-dimensional shape
        Save the results in cluster_track.skele
        
        '''
        cluster_labels = {} # Store the (x, y) -> cluster number info
        binary_image = np.zeros((self.width*self.sr_scale + self.overhead, 
        self.height*self.sr_scale + self.overhead)).astype('uint8')
        
        # record the labels of the cluster and generate the binoary image
        for clu in self.clusterlist:
            cluster_xy = (clu.xy_coordinates()*self.sr_scale).astype('int')
            for xy in cluster_xy:
                pos = tuple(xy)
                binary_image[pos] = 1
                cluster_labels[pos] = clu.num
        if algorithm == 'close':
            closed = closing(binary_image, square(int(sigma*self.sr_scale)))
            second_skele = skeletonize_3d(closed)
        else:
            first_skele = skeletonize_3d(binary_image)*100 # Used for gaussian smoothing
            skele_blurred = gaussian_filter(first_skele, sigma).astype(bool)
            second_skele = skeletonize_3d(skele_blurred)
        
        # Assign the cluster number back to the skeletonised image
        cluster_dict = self._cluster_dict()
        skele_locs = list(zip(*np.where(second_skele > 0)))
        for xy in skele_locs:
            for xy_neighbours in SR_fit._nearby(xy):
                if xy_neighbours in cluster_labels.keys():
                    clu_num = cluster_labels[xy_neighbours]
                    cluster_dict[clu_num].update_skele_list(xy)
                    break
                    
        # save skeletonised image
        skele_image_path = op.join(self.results_root, 'skele.png')
        closed_image_path = op.join(self.results_root, 'closed.png')
        if algorithm == 'close':
            imwrite(skele_image_path, (second_skele*255).transpose())
            imwrite(closed_image_path, (closed*255).transpose())
        else:
            imwrite(skele_image_path, second_skele.transpose())

    def _labelled_image(self):
        '''
        Generate scaled image with aggregates labelled with their aggregate num
        This is for analysis of eccentricity and convex hull
        '''
        
        assert hasattr(self, 'clusterlist')
        labelled_image = np.zeros((self.width*self.sr_scale + self.overhead, 
        self.height*self.sr_scale + self.overhead))
        for clu in self.clusterlist:
            cluster_xy = (clu.xy_coordinates()*self.sr_scale).astype('int')
            for xy in cluster_xy:
                pos = tuple(xy)
                labelled_image[pos] = clu.num
        return labelled_image.astype('int')        
    
    def summarise(self):
        '''
        Output method for summarising the clusters, return the results and also
        save it within the result folder.
        '''
        with open(op.join(self.results_root, 'to_summary.txt'), 'w') as s:
            json.dump(self.to_summary, s, indent=2)
        return self.to_summary
        
    @classmethod
    def _parse_xml(cls, source_string, kw):
        '''
        Static method to extract information of the raw image from the 
        imageJ GDSC SMLM metadata xml file.
        
        '''
        flag = True
        try:
            to_extract = source_string.split('</{}>'.format(kw))[0].split(
                                '<{}>'.format(kw))[1]
        except Exception: flag = False
        if not to_extract: flag = False
        assert flag, "Extract image metadata failed!"
            
        return to_extract
         
    @classmethod
    def _nearby(cls, xy, neighbours=3):
        x, y = xy
        L = []
        for i in range(x-neighbours+1, x+neighbours):
            for j in range(y-neighbours+1, y+neighbours):
                L.append((i, j))

        centerList = [] # To generate the list that searches in the order of 
                        # distance to the center point
        for i in range(len(L)):
            centerList.append(L[len(L)-1-i])
            centerList.append(L[i])  
        return centerList[len(centerList)//2:]

        
class cluster_track(object):
    '''
    This class stores information associated with a DBSCAN cluster
    
    '''
    
    def __init__(self, num, subj, framelen, df):
        self.num = num
        self.df = df
        self.loc_num = len(df)
                
        if subj == 'mol': # for PAINT
            self.mol_num = len(df)
            self.burst_num = sum(df['Burst_number'])
            self.loc_num = sum(df['ON_time'])/framelen
        elif subj == 'burst': # for STORM/PALM
            self.mol_num = len(df['Molecule_ID'].drop_duplicates())
            self.burst_num = len(df)
            self.loc_num = sum(df['ON_time'])/framelen
        else: # legacy: localisation
            self.loc_num = len(df)
    
    def bbox(self):
        a = self.df[['X','Y']]
        box = np.min(a['X']), np.min(a['Y']), \
            np.max(a['X'])-np.min(a['X']), \
            np.max(a['Y'])-np.min(a['Y'])
        return np.array(box)
        
    def update_num(self, update_dict):
        self.num = update_dict[self.num]
        return self
        
    def update_skele_list(self, skele_xy):
        '''
        Store the dataframe of the skeleton of the cluster.
        
        '''
        if not hasattr(self, 'skele'):
            self.skele = DF(columns=['X','Y'])

        self.skele = self.skele.append({'X':skele_xy[0], 'Y':skele_xy[1]}, ignore_index=True)
        
    def xy_coordinates(self):
        return self.df[['X','Y']].values 
        
    def length_measure(self, scale):
        '''
        Measure the length of the skeletonised aggregate.
        Return the length in scaled pixels
        
        '''
        length = 0             
   
        try: # skip arrays with no sample
            xy = self.skele.values
            nbrs = NearestNeighbors(radius = 1.5, algorithm='auto').fit(xy)
            rng = nbrs.radius_neighbors(xy)
            for i in rng[0]:
                length += sum(i)
        except (ValueError, AttributeError):
            pass
        finally:
            length = length/2 + 1
        self.nm_length = length*scale

