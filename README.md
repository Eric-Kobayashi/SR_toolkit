SR_toolkit_v3_2

Issued: 17-Jan-2019 
@author: Eric Kobayashi

Changes over v3_1:

    1. Add compatibility with v3_0, i.e. run cluster fit only.
    
    2. Add fiducial markers removal.
    
    3. Save fiducial marker number in the summary.
    
    4. Fix some bugs in the library.
    
=============================================================================
    
SR_toolkit_v3_1

Issued: 20-Dec-2018
@author: Eric Kobayashi

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

Created: 12-Nov-2018 17:40
First published: 03-Dec-2018

@author: Eric Kobayashi (Based on DRW's SR_toolkit_v2)

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
