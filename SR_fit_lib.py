'''
The library includes the classes and methods related to DBSCAN cluster analysis
and subsequent measurements.

'''

import os
import os.path as op
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings
from math import ceil
from pandas import Series
from pandas import DataFrame as DF
from imageio import imwrite
from skimage.morphology import skeletonize_3d, closing, square
from scipy.ndimage.filters import gaussian_filter
from skimage.measure import regionprops
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN
from collections import OrderedDict


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

        try:
            self.df = pd.read_table(filepath)       # default delimitor: space
            if len(self.df.columns) <= 2:    # if it is a comma delimitor file
                self.df = pd.read_csv(filepath)
        except:
            self._isfit = False          # Wrong filetype
            return
        
        self._isfit = True
        
        if len(self.df) == 0:
            self._empty = True
            return
            
        self._empty = False
        
        # Define logfile to store the analysis detail
        # self.logfile
            
    def input_parameters(self, results_dir=None, pixel_size=None, camera_gain=None,
     camera_bias=None, DBSCAN_eps=None, DBSCAN_min_samples=None, sr_scale=8, 
     frame_length=50, **kwargs):
        '''
        Initialise all the analysis parameters.
        This method is essential.
        
        '''
        
        self.pixel_size = pixel_size
        self.camera_gain = camera_gain
        self.camera_bias = camera_bias
        self.DBSCAN_eps_nm = DBSCAN_eps
        self.DBSCAN_eps = DBSCAN_eps/pixel_size
        self.DBSCAN_min_samples = DBSCAN_min_samples
        self.sr_scale = sr_scale
        self.frame_length = frame_length/1000
        self.results_root = results_dir

        # self.timestring = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        
        # Extract image information from image.results.xls
        image_details_file = "Image.results.xls"
        image_details = op.join(self.root, image_details_file)
        assert op.isfile(image_details), "Image metadata file not found!"
        with open(image_details, 'r') as metadata:
            metadata_content = metadata.read()
            self.framenum = int(SR_fit._parse_xml(
                                metadata_content, 'frames'))
            self.width = int(SR_fit._parse_xml(
                                metadata_content, 'width'))
            self.height = int(SR_fit._parse_xml(
                                metadata_content, 'height'))
                                
    def cluster_info(self):
        '''
        Extract information from the cluster analysis results.
        This method is essential. 
        Replace 'easy_numbers' in v2.
        
        '''
        
        if not self._isfit or self._empty:
            # No localisations
            return
        
        cluster_analysis_file = "DBSCAN_eps_{:.2f}nm_minlocs_{}.txt".format(
                                self.DBSCAN_eps_nm, self.DBSCAN_min_samples)
        cluster_analysis_file_path = op.join(self.results_root, 
                                    cluster_analysis_file)
        if op.isfile(cluster_analysis_file_path):
            # Cluster analysis is already done
            self.df = pd.read_table(cluster_analysis_file_path)
        else:
            self._cluster_analysis(self.DBSCAN_eps, self.DBSCAN_min_samples)
            self.df.to_csv(cluster_analysis_file_path,
                index=False, sep='\t')
            
        clustered_locs = self.df[self.df.Cluster > 0]
        self.unclustered_df = self.df[self.df.Cluster == 0]
        self.df = clustered_locs
        if len(clustered_locs) == 0:
            # No clustered localisations
            self._clustered = False
            return
        else:
            self._clustered = True
            
        self.clusterlist = []
        cluster_labels = sorted(clustered_locs.Cluster.unique())
        clu_num = 1
        for clu_labels in cluster_labels:
            self.clusterlist.append(cluster_track(clu_num,
            clustered_locs[clustered_locs.Cluster == clu_labels]))
            clu_num += 1  # Cluster analysis could miss some cluster number label
    
    def burst_filter(self, min_burst, fill_gap=50):
        '''
        Filter the bursts that have less than min_burst bursts
        Recommend fill_gap to be set to 50, for removing blinking non-specific bindings
        Single bursts are by default removed
        '''
        with pd.option_context('mode.chained_assignment', None):
                # Suppress the SettingwithCopyWarning 
            
            if self._isnan():
                # No localisations
                self.filtered_cluster = 0
                return
    
            #print('\t\tNumber of spot before filtering: {}'.format(n_spot))
            remain_cluster = []
            for clu in self.clusterlist:
                clu.burst_profile(fill_gap=fill_gap, frame_length=self.frame_length)
                if clu.burstnum >= min_burst:
                    remain_cluster.append(clu.num)
            self.filtered_cluster = len(self.clusterlist)-len(remain_cluster)
            
            #print('\t\tNumber of spot after filtering: {}'.format(n_spot))
            output_df = self.df[self.df.Cluster.isin(remain_cluster)].copy()
            # rename all the cluster number to 1, 2, 3,...
            cluster_dict = {v: k+1 for k, v in enumerate(remain_cluster)}
            output_df['Cluster'] = output_df['Cluster'].map(cluster_dict)
            self.clusterlist = [clu.update_num(cluster_dict) for 
                clu in self.clusterlist if clu.num in remain_cluster]
            self.df = output_df
            # 
            # output_df.to_csv(op.join(self.root,
            # 'DBSCAN_Results_filtered_minburst_{}_fillgap_{}.txt'.format(min_burst, fill_gap)), 
            # index = False, sep = '\t')
            output_df.to_csv(op.join(self.results_root, 'DBSCAN_filtered.txt'),
                index = False, sep = '\t')
    
    def burst_info(self, fill_gap, remove_single):
        if self._isnan():
            # No localisations
            self.ave_burstlen, self.ave_burstnum, \
                self.ave_darklen, self.ave_lighttodark = (np.nan,)*4
            return
        
        for clu in self.clusterlist:
            clu.burst_analysis(fill_gap=fill_gap, max_frame=self.framenum, 
                frame_length=self.frame_length, remove_single=remove_single)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self.ave_burstlen = np.nanmean([clu.ave_burstlen for clu in self.clusterlist])
            self.ave_burstnum = np.nanmean([clu.burstnum for clu in self.clusterlist])
            self.ave_darklen = np.nanmean([clu.ave_darklen for clu in self.clusterlist])
            self.ave_lighttodark = np.nanmean([clu.lighttodark for clu in self.clusterlist])
        
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
                
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self.ave_length = np.nanmean([clu.nm_length for clu in self.clusterlist])
            self.ave_density_1D = np.nanmean([clu.loc_num / clu.nm_length for 
                                            clu in self.clusterlist])
                                 
    def eccentricity_measure(self):
        if self._isnan():
            # No localisations
            self.ave_ecc, self.ave_flattening = (np.nan,)*2
            return

        self._regionprops()
        for clu in self.clusterlist:
            setattr(clu, 'ecc', self.props[clu.num-1]['eccentricity'])
            setattr(clu, 'flattening', 1 - np.sqrt(1 - clu.ecc**2))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self.ave_ecc = np.nanmean([clu.ecc for clu in self.clusterlist])
            self.ave_flattening = np.nanmean([clu.flattening for clu in self.clusterlist])

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
            setattr(clu, 'nm2_area', self.props[clu.num-1]['convex_area']*(scale**2))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self.ave_area = np.nanmean([clu.nm2_area for clu in self.clusterlist])
            self.ave_density_2D = np.nanmean([clu.loc_num / clu.nm2_area for 
                                            clu in self.clusterlist])
            
    def update_ID(self, keystring):
        self.fit_ID = keystring
        
    def labelled_cluster(self):
        if self._isnan():
            # No localisations
            return

        roi = DF([np.array(list(map(ceil, (clu.bbox()*self.sr_scale)))) 
            for clu in self.clusterlist])
        roi_file = op.join(self.results_root, 'clusters_roi.txt')
        roi.to_csv(roi_file, index=False)
             
    # Output methods
    def save_with_header(self):
        if self._isnan():
            # No localisations
            return

        header = '''#Localisation Results File
#FileVersion Text.D0.E0.V2

#Name Image (LSE)
#Source <gdsc.smlm.ij.IJImageSource><singleFrame>0</singleFrame><extraFrames>0</extraFrames><path>{0}</path></gdsc.smlm.ij.IJImageSource>
#Bounds x0 y0 w{1} h{2}
#Calibration <gdsc.smlm.results.Calibration><nmPerPixel>{3}</nmPerPixel><gain>{4}</gain><exposureTime>{5}</exposureTime><readNoise>0.0</readNoise><bias>{6}</bias><emCCD>false</emCCD></gdsc.smlm.results.Calibration>
#Configuration <gdsc.smlm.engine.FitEngineConfiguration><fitConfiguration><fitCriteria>LEAST_SQUARED_ERROR</fitCriteria><delta>1.0E-4</delta><initialAngle>0.0</initialAngle><initialSD0>2.0</initialSD0><initialSD1>2.0</initialSD1><computeDeviations>false</computeDeviations><fitSolver>LVM</fitSolver><minIterations>0</minIterations><maxIterations>20</maxIterations><significantDigits>5</significantDigits><fitFunction>CIRCULAR</fitFunction><flags>20</flags><backgroundFitting>true</backgroundFitting><notSignalFitting>false</notSignalFitting><coordinateShift>4.0</coordinateShift><signalThreshold>1665.0</signalThreshold><signalStrength>30.0</signalStrength><minPhotons>30.0</minPhotons><precisionThreshold>625.0</precisionThreshold><precisionUsingBackground>false</precisionUsingBackground><nmPerPixel>{3}</nmPerPixel><gain>{4}</gain><emCCD>false</emCCD><modelCamera>false</modelCamera><noise>0.0</noise><widthFactor>2.0</widthFactor><fitValidation>true</fitValidation><lambda>10.0</lambda><computeResiduals>false</computeResiduals><duplicateDistance>0.5</duplicateDistance><bias>{6}</bias><readNoise>0.0</readNoise><maxFunctionEvaluations>1000</maxFunctionEvaluations><searchMethod>POWELL</searchMethod><gradientLineMinimisation>false</gradientLineMinimisation><relativeThreshold>1.0E-6</relativeThreshold><absoluteThreshold>1.0E-16</absoluteThreshold></fitConfiguration><search>3.0</search><border>1.0</border><fitting>3.0</fitting><failuresLimit>10</failuresLimit><includeNeighbours>true</includeNeighbours><neighbourHeightThreshold>0.3</neighbourHeightThreshold><residualsThreshold>1.0</residualsThreshold><noiseMethod>QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES</noiseMethod><dataFilterType>SINGLE</dataFilterType><smooth><double>0.5</double></smooth><dataFilter><gdsc.smlm.engine.DataFilter>MEAN</gdsc.smlm.engine.DataFilter></dataFilter></gdsc.smlm.engine.FitEngineConfiguration>
#'''.format(self.path, self.width, self.height, self.pixel_size, self.camera_gain, self.frame_length, self.camera_bias)
        
        hf = op.join(self.results_root, 'All_fits_with_header_py.txt')
        to_output = self.df.drop(columns=['Source', 'Cluster', 'SNR'], errors='ignore')
        to_output.to_csv(hf, sep = '\t',  index = False)
        with open(hf, 'r+') as log:
            content = log.read()
            log.seek(0)
            log.write(header+content)
                
    def loc_num(self):
        return len(self.df)
        
    def cluster_num(self):
        if hasattr(self, 'clusterlist'):
            return len(self.clusterlist)
        else:
            return 0
        
    def ave_cluster_locnum(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            if hasattr(self, 'clusterlist'):
                return np.nanmean(self.output_histogram('loc_num'))
            else:
                return np.nan
            
    def output_histogram(self, attr):
        List_of_attr = []
        if hasattr(self, 'clusterlist'):
            for clu in self.clusterlist:
                assert hasattr(clu, attr)
                List_of_attr.append(getattr(clu, attr))
        return List_of_attr
        
    def output_noncluster_histogram(self, attr):
        List_of_attr = []
        List_of_clusternum = []
        if hasattr(self, 'clusterlist'):
            for clu in self.clusterlist:
                assert hasattr(clu, attr)
                attrlist = getattr(clu, attr)
                List_of_attr += attrlist
                List_of_clusternum += [clu.num]*len(attrlist)
        return List_of_clusternum, List_of_attr
        
    def output_attr_fromfile(self, attrfile):
        attrpath = op.join(self.root, attrfile)
        if not op.isfile(attrpath):
            return np.nan
        else:
            with open(attrpath, 'r') as f:
                attr = float(f.read())
            return attr
        
    # Private methods
    def _isnan(self):
        if not self._isfit or self._empty or not hasattr(self, 'clusterlist'):
            return True
        elif len(self.clusterlist) == 0:
            return True
        else:
            return False

            
    def _cluster_analysis(self, epsilon, min_samples):
        '''
        Do DBSCAN cluster analysis on the dataframe.
        This method is only called inside the class because it is essential 
        to all the following analysis.
        
        '''       
        fit_xy = self.df[['X', 'Y']].values
        db = DBSCAN(eps=epsilon, min_samples=min_samples).fit(fit_xy)
        self.df['Cluster'] = db.labels_ + 1  # Start labelling with 1, makes analysis easier

            
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
        overhead = 50 # Make extra room for fiducial corrections
        binary_image = np.zeros((self.width*self.sr_scale + overhead, 
        self.height*self.sr_scale + overhead)).astype('uint8')
        
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
        overhead = 50 # Make extra room for fiducial corrections
        labelled_image = np.zeros((self.width*self.sr_scale + overhead, 
        self.height*self.sr_scale + overhead))
        for clu in self.clusterlist:
            cluster_xy = (clu.xy_coordinates()*self.sr_scale).astype('int')
            for xy in cluster_xy:
                pos = tuple(xy)
                labelled_image[pos] = clu.num
        return labelled_image.astype('int')        
      
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
    This class stores information associated with a localisation cluster
    
    '''
    
    def __init__(self, num, df):
        self.num = num
        self.df = df
        self.loc_num = len(df)
        self.frames = np.array(df.Frame)
        self.frames.sort()        # sort the frame numbers
    
    def bbox(self):
        a = self.df[['X','Y']]
        box = np.min(a['X']), np.min(a['Y']), np.max(a['X'])-np.min(a['X']), np.max(a['Y'])-np.min(a['Y'])
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
  
    def burst_profile(self, fill_gap, frame_length, remove_single = True):
        '''
        Generate time profile of the burst track at a certain fill_gap
        fill_gap: the minimum number of frames between each burst
        remove_single: if Ture, remove the burst that has only one frame
        
        ''' 
        F = self.frames
        worktrack = [0]  # label the burst each element in F belongs to
        burst = 0
        for i in range(1, len(F)):
            if (F[i] - F[i-1]) < fill_gap:  # Fill in the gap between the frames 
                worktrack.append(burst)
            else:
                burst += 1
                worktrack.append(burst)
        self.burstnum = burst + 1
        worktrack = Series(worktrack)
        bursts = {}     # store the #burst -> burstlength details
        burstframes = OrderedDict()
        for i in range(burst + 1):
            rangeburst = worktrack[worktrack == i].index
            lenburst = F[rangeburst.max()] - F[rangeburst.min()] + 1
            bursts[i] = lenburst
            burstframes[i] = (F[rangeburst.min()], F[rangeburst.max()])
        self.allburstframe = burstframes

        if remove_single:  # bursts that last only one frame are most likely to 
                            # be camera noises
            temp_dict = bursts.copy()
            for numburst, lenburst in bursts.items():
                if lenburst <= 1:                                          
                    del temp_dict[numburst]
                    del self.allburstframe[numburst]
            bursts = temp_dict.copy()
        self.alllen = [frame_length*l for l in bursts.values()]
        if bursts == {}: 
            self.ave_burstlen = np.nan
            self.burstnum = 0
            return
        # elif sum(self.alllen) > 2000:
        #     self.ave_burstlen = 0
        #     self.burstnum = 0
        #     self.avgdark = 0
        #     self.lighttodark = 0
        #     print("Constant burst detected, possibly a fiducial marker. The position is X = %d, Y = %d" %self.positions)
        #     return None
        self.burstnum = len(self.alllen)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self.ave_burstlen = np.nanmean(self.alllen)
        
    def burst_analysis(self, fill_gap, max_frame, frame_length, remove_single = True):
        '''
        Measure the burst number, burst length, dark number, dark length
        of each track.
        
        '''
        
        self.burst_profile(fill_gap=fill_gap, frame_length=frame_length, 
            remove_single=remove_single)
        darkwork = []
        if self.burstnum == 0: 
            self.lighttodark = np.nan
            self.ave_darklen = np.nan
            self.alldarklen = []
            return
        for k, v in self.allburstframe.items():
            darkwork.append(v[0])
            darkwork.append(v[1])
        darkwork = [0] + darkwork + [max_frame + 1]
        # generate the dark begin and end
        alldark = [darkwork[i:i+2] for i in range(0, len(darkwork),2)] 
        # the number of dark frame is end - begin - 1
        self.alldarklen = [frame_length*(i[1] - i[0] - 1)
                            for i in alldark if i[1] != (i[0] + 1)]          
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            self.ave_darklen = np.nanmean(self.alldarklen)
        # self.darknum = len(self.alldarklen)
        if self.ave_darklen != 0:
            self.lighttodark = float(self.ave_burstlen)/self.ave_darklen
        else:
            self.lighttodark = float(max_frame)     
        
    def length_measure(self, scale):
        '''
        Measure the length of the skeletonised aggregate.
        Return the length in scaled pixels
        
        '''
        
        #assert hasattr(self, 'skele'), "Run skeletonisation before measuring length!"
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

#     def intensity_profile(self, max_frame):
#         multidots = {}
#         for f, g in [(_, g) for _, g in self.intensityprofile.groupby('Frame') if len(g) > 1]:
#             multidots[f] = g['origValue'].mean()
#         to_graph = self.intensityprofile[self.intensityprofile.Frame.map(lambda x: x not in multidots.keys())].copy()
#         for keys, values in multidots.items():
#             to_graph.append({'Frame':keys, 'origValue':values}, ignore_index=True)
#         for f in range(max_frame + 1):
#             if f not in to_graph.Frame.values:
#                 to_graph = to_graph.append({'Frame':f, 'origValue':0}, ignore_index=True)
#         to_graph.sort_values(by='Frame', inplace=True)
#         
#         temp = 0
#         filldict = {}
#         for index, row in to_graph.iterrows():
#             if row['origValue'] != -1:
#                 temp = row['origValue']
#             else:
#                 filldict[row['Frame']] = temp
#         to_graph_fill = to_graph[to_graph.Frame.map(lambda x: x not in filldict.keys())].copy()
#         for keys, values in filldict.items():
#             to_graph_fill = to_graph_fill.append({'Frame':keys, 'origValue':values}, ignore_index=True)
#         to_graph_fill.sort_values(by='Frame', inplace=True)
#         to_graph_fill.plot(x='Frame', y='origValue', kind='line', linewidth=1, color='r', figsize=(30,10))
#         plt.title('Trace at x={}, y={}'.format(int(self.X), int(self.Y)))
#         # to_graph.to_csv(os.path.join(savepath, 'fill_gap_raw.csv'), index=False)
#         # plt.savefig(os.path.join(savepath, 'fill_gap.pdf'), format='pdf')
#         plt.close()
# 
#         return to_graph, to_graph_fill
        
