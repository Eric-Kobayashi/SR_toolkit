# SR_toolkit
Integrated ImageJ and Python workflow for super-resolution fluoresence image analysis.

Author: Eric Hidari

Version: 3.6

Issue Date: 02 June 2021

Summary
----
This workflow software uses [GDSC SMLM ImageJ plugin](https://gdsc-smlm.readthedocs.io/en/latest/) for localisation fitting/fiduical drift correction/rendering super-res images, [Redundant cross-correlation](https://github.com/yinawang28/RCC) for correlation based drift correction, and DBSCAN for cluster analysis.

The software supports both GDSC SMLM 1 and GDSC SMLM 2, however some functions are not available for GDSC SMLM 2.

Before you start
----
Make sure ImageJ and the GDSC SMLM plugin are installed. If you wish to use the redundant cross-correlation drift correction, make sure MATLAB (version not earlier than 2018b) is installed and included in system PATH.

Run the `installer.exe` file, fill in the requested paths and press `Install`. It will copy the source code into the destination folder, create configuration file, and set the `PYTHONPATH` to the destination folder. 

Parameters to be set
----
- `directory`: dirctory containing the images to be analysed.

- `image_condition`: search criteria for images. If set to `'separate_images'`, images in each folder would be combined into a TIF stack before the analysis.

- `pixel_size`: pixel size in nanometer.

- `sr_scale`: default set to 8, scale for generating super-resolution images.

- `camera_bias`: in ADU, raw values recorded when no photons detected.

- `camera_gain`: ADUs/photon, see [here](https://gdsc-smlm.readthedocs.io/en/latest/fitting_plugins.html#imaging-calibration-parameters) for details.

- `frame_length`: ms per frame for exposure time.

- `GDSC_SMLM_version`: the version of GDSC SMLM used for the analysis.

- `GDSC_SMLM_peak_fit`: controls whether to run the GDSC SMLM peak fit in ImageJ or not, if False, only cluster analysis will be run.

- `BG_measurement`: controls whether to measure the background of the image stack (median pixel value). Value stored in `median_intensity.txt`.

- `trim_track`: `run` - whether to run the track trimmer. `frame_number` - how many frames remains for the analysis (the rest are trimmed). `from_end` - whether to keep the frames from the end or the beginning.

- `signal_strength`, `precision`, `min_photons`: parameters set for GDSC SMLM peak fit filters, see [here](https://gdsc-smlm.readthedocs.io/en/latest/fitting_plugins.html#filtering-parameters) for details. `precision` is in nm.

- `spot_filter`: only used with GDSC SMLM 2. `filter_type` can be either `'Single'` or `'Difference'`, see [here](https://gdsc-smlm.readthedocs.io/en/latest/fitting_plugins.html#spot-filter-type) for details. `smoothing` is the size of the first smoothing filter, and `smoothing2` is the size of the second filter for the `'Difference'` spot filter.

- `fiducial_correction`: `run` - whether to run the drift correcter. `correction_method` - the method used to correct the drift. Four modes are supported for GDSC SMLM 1: `'auto_fid'` (automatic fiducial correction), `'drift'` (drift file), `'fid_file'` (predefined fiducial coordinates), and `'corr'` (redundant cross-correlation based drift correction); while only `'auto_fid'`, `'drift'` and `'fid_file'` are supported for GDSC SMLM 2. Furthermore, the drift plot and corrected image rendering are not supported for GDSC SMLM 2. `smoothing_para` and `limit_smoothing` control the smoothing of curve for the correction, see [here](https://gdsc-smlm.readthedocs.io/en/latest/analysis_plugins.html#drift-calculator) for details.

- Here the 4 correction modes for `fiducial_correction` are briefly explained.
    - `'auto_fid'`: locate the fiducial markers in the FOV by thresholding out the pixels with intensity higher than `fid_brightness`, draw a square with `fid_size` around the center, and select fiducials lasting longer than `fid_last_time` in frames. Then it uses the GDSC SMLM builtin [Drift calculator](https://gdsc-smlm.readthedocs.io/en/latest/analysis_plugins.html#drift-calculator) to correct the drift with the method setting `Fiducial markers within an image`.

    - `'fid_file'`: directly uses the predefined fiducial location in XY coordinates to correct the drift using the above method. Each image must be contained as the only image file in a folder, and the XY coordinates must be stored in the file `Fiducials.txt` within the folder, with the format as (delimited by tab):
    ```
    X    Y
    100    100
    200    200
    ...
    ```

    - `'drift'`: directly applies the predefined drift profile, and uses [Drift calculator](https://gdsc-smlm.readthedocs.io/en/latest/analysis_plugins.html#drift-calculator) to correct the drift with the method setting `Drift file`. Each image must be contained as the only image file in a folder, and the drift file must be stored as `Drift.txt` within the folder, with the format as (delimited by tab):
    ```
    Time    X    Y
    1    0    0
    2    0.1    0.1
    3    0.2    0.2
    ...
    ```

    - `'corr'`: calculate the drift using redundant cross-correlation in MATLAB, save it in the drift file and correct the drift using the above method. `bin_size` sets the spatial bin pixel size in nm, and `segpara` sets the segmentation parameter in frame. See [here](https://github.com/yinawang28/RCC) for details and the original paper.


> The rest of the parameters are for the python code of cluster analysis and characterisation. The detailed explanation of the analysis can be found in my thesis section 2.3.2 and 2.3.3, which can be downloaded [here](https://www.repository.cam.ac.uk/handle/1810/315556).

- `fitresults_file_name`: search criteria for fit results files. If set to `'default'`, it will be automatically set based on whether fiducial correction is enabled.

- `temporal_grouping`: group neighbouring localisations (in time and space) into the same burst. `run` - whether to run the temporal grouping. `dThresh` and `min_loc` control the DBSCAN parameters for grouping localisations in space as molecules; while `tThresh` and `min_frame` are the ones for grouping localisations in the time dimension as bursts. `max_mol_area`, `min_on_prop` and `min_burst` are extra filters can be applied: `max_mol_area` (in nm$^2$) filters out localisations incorrectly grouped in space and thus occupying too large of an area; `min_on_prop` filters out localisations that are incorrectly grouped in time and thus leaving too long between bursts. `'min_burst'` filters out bursts that blink less than the threshold during the track.

- `cluster_analysis`: group neighbouring subjects (`cluster_subject`, either `'burst'` for bursts generated in temporal grouping, or `'loc'` for the raw localisations) into clusters that represent protein aggregates. `run` - whether to run the cluster analysis. `DBSCAN_eps_nm` and `DBSCAN_min_samples` - DBSCAN parameters.

- `length_measure`: `run` - whether to run the length measurement. `algorithm` - algorithm to treat holes in skeletonised aggregate structures. `'blur'` or `'close'`; if `'blur'`, `sigma` is the Gaussian sigma; if `'close'`, `sigma` is the closing square size.

- `eccentricity_measure` and `convexhull_measure`: whether to run the eccentricity and convex hull area measurements.

- `Rendering_SR`: method to render the super-resolution images after cluster analysis. Set to `False` for no rendering, `'loc'` for rendering with `localisations` only, `'loc_pre'` for rendering with `localisations (width = precision)`. See [here](https://gdsc-smlm.readthedocs.io/en/latest/fitting_plugins.html#results-parameters) for details.

- `save_histogram`: whether to save the information of individual molecules and bursts (generated in temporal grouping) into summary tables. Set to `False` for large imaging experiments can significantly save time and storage space.



-----




Past versions 
----------------
### SR_toolkit_v3_6

Issued: 02-06-2021

Changes over v3_5:

1. Add a GUI installer, including automated installation of Python dependencies.

2. Add support for cross-correlation fiducial correction.

3. Add support for reading images of different folder.

4. Save drift profile plots.

5. Fix the conflict with newer versions of ImageJ.

6. Remove symlink and variable background functions.

7. Document the parameter settings.


### SR_toolkit_v3_5

Issued: 28-04-2019

Changes over v3_4:
    
#### Major changes:
        
1. Add temporal grouping: group localisations into molecules, then group
them into individual bursts (using 2D and 1D DBSCAN, respectively.) Filters
can be applied on burst number of molecule, molecule size and on time
proportion of burst (blinking molecules can be removed). The X, Y positions
and precisions in each group are averaged to enable further clustering.

2. Add DBSCAN clustering with molecules or bursts. DBSCAN with bursts will 
be suitable with STORM/PALM, while molecule suits with PAINT.

3. Remove the burst filter and burst analysis because temporal grouping 
does the same job.
    
#### Minor changes:
    
1. Move to_ImageJ.json to individual analysis folders, allowing multiple 
analysis conducted at the same time. 

2. Simplify the compatiablity code with GDSC SMLM 2.

3. Simplify the data saving code.



### SR_toolkit_v3_4

Issued: 10-03-2019

> SR_toolkit_v3_5 is under development with a more accurate clustering strategy and therefore 
more punctuate reconstruction of super-resolution imge thanks to @StevetheChemist's advice.
If you would like to use this toolkit for publication, please wait for v3_5 (before end of 
April 2019)

Changes over v3_3:

1. Add support of GDSC SMLM 2, although the package is still an unstable release,
it gives the difference Gaussian filter which could be useful. The other parameters
have not yet been included in SR_toolkit however if any of them is significant it 
will be added in a future version.

2. Add setup file for easy installation.

3. Simplify the code of generating fit header file.    

### SR_toolkit_v3_3

Issued: 01-Feb-2019

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
    
> Note: Please put Rendering_SR.py in the same folder as GDSC_SMLM_peak_fit.py


=============================================================================

SR_toolkit_v3_2

Issued: 17-Jan-2019 

Changes over v3_1:

    1. Add compatibility with v3_0, i.e. run cluster fit only.
    
    2. Add fiducial markers removal.
    
    3. Save fiducial marker number in the summary.
    
    4. Fix some bugs in the library.
    
=============================================================================
    
SR_toolkit_v3_1

Issued: 20-Dec-2018

Major changes over v3:

    1. Incorporate ImageJ GDSC SMLM peak fit and the brightness-based fiducial
    correction code into the analysis. Re-organise the main run_fit function
    so it workflows the ImageJ part and the python part. The ImageJ part is
    'lazy'--- it only runs if a previous analysis using the same parameters
    is not found in the same root path. Notice that if only fitresults file are 
    available (images are not available), this script will not run properly. 
    Please use SR_toolkit_v3_0 for that. 
    
    2. Trimming the frame to a certain frame_number is supported. This allows a
    certain number of frames to be processed. 
    
    3. Cluster analysis measurement can be turned off. This allows only GDSC SMLM
    to run and therefore a quick check of the localisation number.
    
    4. A symbolic link to the original image stack can be generated in the 
    analysis folder. This feature needs python to be run in administrator mode.
     
Some minor changes:    

    1. All the analysis now are by default copied. This means no analysis files
    will be appended to the raw image folder. 
    2. Use JSON to output the logging information rather than xml. More human
    readable.
    3. Fix some bugs in the library.
    
Potential features update:

    1. Integrate more ImageJ analysis including GDSC SMLM, background correction, 
    find maxima (ThT counting), rendering SR images into the main analysis.
    2. Doing DBSCAN on bursts (events lasting multiple frames count as one burst)
    rather than localisations can be more precise.
    3. Use burst analysis for fiducial correction. (need to check papers)
    4. Length measurement is slow. Could be improved.
    5. Integrate sqlite to pack the outputs into database automatically. (ongoing)

===============================================================================

SR_toolkit_v3_0

First published: 03-Dec-2018

Major changes over v2:

    1. Pack the functions in the v2 scripts into 3 classes and their methods:
    a. Analysis class - logs the analysis parameters, collects and saves the
    outputs, with the flexibility of packing more analysis in the class.
    A log file is saved to keep track of the parameters used in the analysis.
    
    b. SR_fit class - represents a fit result file aka. an image.
    Records the parameters passed by the analysis class, runs the 
    DBSCAN cluster and the subsequent characterisations of the clusters.
    
    c. cluster_track class - for each cluster identified in the cluster 
    analysis. Stores the characteristics of the cluster.
    
    2. Add functionalities: burst filter and burst analysis
    These methods analyse the time dimension of each cluster with a certain 
    fill_gap parameter to account for dye blinking.
       
    Burst filter (set with a higher fill_gap) tries to remove clusters that 
    have less than N bursts. This is assuming the real events are recurring 
    and the non-specific bindings are random, non-recurring. One could play 
    around the balance between min_burst and DBSCAN_min_samples options to 
    see how the results changes.
    
    Burst analysis (usually set with a lower fill_gap) measures the on and 
    off time of each cluster. This is related with qPAINT analysis, which 
    determines how many labelled molecules are at one spot. See Jungmann paper
    for details. The remove_single option removes bursts that only lasts for
    1 frame. They are likely to be camera noises.
        
    3. Remove Nearest neighbour analysis: The NN analysis made assumption such 
    as defining 5 times nearest neighbour distances as the 'neighbourhood';
    the nearest neighbour distance is not really informative since the nearest 
    neighbours usually come from the same binding spanning in mulitple frames.
    
    Since the anlysis was intended to measure the density of localisations,
    convex hull area is added here so that the 2D density of the aggregates
    can be estimated. 
    
    More elegant methods, such as a 'Concave hull' (with a smoothing constant)
    or Ripley's k function could be added in the features.
    
Some minor changes:

    1. Remove the generate SR_mage, ThT maxima, cropping functionality. 
    The former two are more convenient done in imageJ, the latter will be added
    in a future version, cropping either space or time.
    
    2. The new analysis makes some uses of the Image.results.xls file which 
    records the frame number, width and height. However it will throw an error
    if the file does not exist. 
    
    3. Timestamp is used to log the analysis, giving each image a unique ID 
    which will be handy in SQL operations. A relationship can be easily built 
    between the summary and histograms file. 
    
    4. A folder will be generated to hold the analysis outputs everytime an 
    Analysis object is generated. There is also an option of copying all analysis 
    generated files to the results folder. This avoids overwriting the previous 
    results.
    
    5. More to cover...


Potential features update:
    1. Create symlink for raw images.
    2. Integrate ImageJ analysis including GDSC SMLM, background correction, 
    find maxima (ThT counting), fiducial correction into the main analysis.
    3. Doing DBSCAN on bursts (events lasting multiple frames count as one burst)
    rather than localisations can be more precise.
    4. Use burst analysis for fiducial correction.
    5. Length measurement is slow. Could be improved.
    6. Integrate sqlite to pack the outputs into database automatically.
    7. See minor changes 1, major changes 3.
    8. Output JSON log file
