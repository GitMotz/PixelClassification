---
output:
  pdf_document: default
  html_document: default
---
This code relates to

##### A systems-level study reveals regulators of membrane-less organelles in human cells
by Doris Berchtold, Nico Battich, and Lucas Pelkmans.  

*** 
  
#### Dependencies
*ClassifyPixels.m* and supporting functions have been tested on Windows and MATLAB 2014a.  

#### Example data set
The example data set comprises 10 images of unperturbed cells of the splicing speckle screen that were stained with antibodies against SRRM2 (in *TIFF*), a measurement file to correct these images for uneven illumination, and the segmentations of cells and nuclei (in *SEGMENTATION*). For comparison, the folder *segmented_splicing_speckles* contains the segmentation of splicing speckles obtained with the model that was trained on the screen.  

#### How to train a segmentation model with the ClassifyPixels GUI and apply it on the example data set  

1.	Launch MATLAB 2014a on a Windows machine  
2.	Add the folder *PixelClassification* to your MATLAB path  
3.	Open *CreateSettingsFile.m* and run  
4.	Launch the GUI by typing *ClassifyPixels* in the MATLAB command line  
5.	Read through the detailed GUI description below. In brief:  
6.	Load the settings file, click next image  
7.	Train the classifier by clicking on single pixels of multiple images and mark them as signal or background  
8.	Save the classification model (default: *PCProject_01.mat*)  
9.	Open *run_pixel_classification_example.m*, change the name of the classification model and run the script. This will apply the trained model on all images of the example data set. The segmentation files will be saved into the folder *PixelClassification\ExampleDataSet\TEST_seg*.  

#### Detailed description of the ClassifyPixels Graphical User Interface  
##### Top section  
* In */path/to/settings* insert the full path to the settings file  
* In *PCProject01.mat* insert the name of the output file, which will contain the classification model
* Click *Load* to load the settings file
* *Save Classification* saves the classification model after training  

##### Bottom section
* Click *Next image* to display the first or next image
* Set the *Low Quantile* and *High Quantile* and click *refresh* to rescale the intensities of the displayed image
* *Go to section* allows to move to another field of view of the same image  

##### Drawing tools
* Select *Signal* (green) before clicking on pixels of the object to segment; select *Background* (red) before clicking on background pixels; select *eraser* to remove and correct wrongly classified pixels 
* The *Brush* scale bar allows to change the number of pixels selected upon clicking on the image  

##### Modeling
All segmentations of MLOs in Berchtold et al. were achieved using the *Gaussian* option of SVMs and the *Smooth Probability* checked. Other implementations have not been tested extensively.  

* Select the type of SVM model that will be trained on the selected pixels
* Click *Train* after marking sufficient pixels as signal or background
* Click *Classify* to classify all pixels of the current image; to visualize the result click *Show Prob*, see below
* The *Smooth Prob* and *More* options allow smoothing of the estimated probabilities of being an object 
* Click *Optimize (SVM only)* to run a feature selection to optimize the misclassification error; this option only works for *Linear*, *Gaussian*, or *Polynomial* modeling  

##### Visualization
* Define the probability threshold for segmentation in the *Seg* box; default = 0.5
* Click *Show Image* to see the current image
* Click *Show Prob* to see the current probability of pixels of being an object; click *Classify* before 
* Click *Show Seg* to see the current segmentation of objects according to the *Seg* threshold
* Click *Original Image* to generate a new figure of the full original image
* Check *Show cell seg* to also plot the cell segmentations
* Click *Plot Features* to generate figures of all pixel features used for classification; looking through these features helps to select the most informative pixel features for object segmentation  

*** 
  
#### Pixel features used for classification
**Feature name** : **Description**  
**log_5_2**: Laplacian of Gaussian filter, size = 5 pixels, sigma = 2  
**log_6_3**: Laplacian of Gaussian filter, size = 6 pixels, sigma = 3  
**log_8_4**: Laplacian of Gaussian filter, size = 8 pixels, sigma = 4  
**log_10_5**: Laplacian of Gaussian filter, size = 10 pixels, sigma = 5  
**log_20_10**: Laplacian of Gaussian filter, size = 20 pixels, sigma = 10  
**gau_20_10**: Gaussian filter, size = 20 pixels, sigma = 10  
**gau_9_4**: Gaussian filter, size =9 pixels, sigma = 4  
**gau_5_2**: Gaussian filter, size = 5 pixels, sigma = 2  
**diffgau_1**: difference of Gaussians, 1st size = 15 pixels, sigma = 8, 2nd size = 15 pixels, sigma = 3   
**diffgau_2**: difference of Gaussians, 1st size = 15 pixels, sigma = 9, 2nd size = 15 pixels, sigma = 5  
**average_3**: averaging filter, size = 3 pixels  
**average_5**: averaging filter, size = 5 pixels  
**average_10**: averaging filter, size = 10 pixels  
**average_20**: averaging filter, size = 20 pixels  
**edges**:  Sobel and Prewitt filters, and edges on Laplacian transformations  
**simple_edges**: Sobel and Prewitt filters  
**edge_log_s8_4**: Prewitt filter on Laplacian-of-Gaussian transformation  
**hole_log_8_4**: A hollow filter of 5 pixel size applied on the Laplacian-of-Gaussian transformation  
**entropy_log_8_4**: Entropy filter applied on the Laplacian-of-Gaussian transformation size = 8 pixels, sigma = 4  
**entropy_log_5_2**: Entropy filter applied on the Laplacian-of-Gaussian transformation size = 5 pixels, sigma = 2  
**morph_opening**: Morphological opening with filter sized of 5, 8, and 10 pixels  
**morph_close**: Morphological opening with filter sized of 3, 5, 7, 10,15, and 20 pixels  
**morph_close_small**: Morphological opening with filter sized of 3, and 5 pixels  
**loc_back_small** and **loc_back_large**: local background estimation  
**deconv_10_7**: Lucy-Richardson deblur with Gaussian filter size = 10 pixels, sigma = 7  
**deconv_10_3**: Lucy-Richardson deblur with Gaussian filter size = 10 pixels, sigma = 3  
**deconv_5_3**: Lucy-Richardson deblur with Gaussian filter size = 5 pixels, sigma = 3  
**std_dif_1**:  standard deviation differences of deblur image  
**zscore_1**: z-scoring filters  
**watershed_lines**, **watershed_lines_smooth**, **watershed_lines_smooth_low**, and **watershed_lines_smooth_very_low**: Filters based on watershed lines  
**tophat_12**: Tophat filter of size 12 applied on Gaussian filtered image  
**orig_tophat_30**: Tophat filter of size 12 applied on original image  

*** 
  
##### Pixel features used for classification of MLOs in Berchtold et al.
**Cajal bodies**:  
log_5_2, log_6_3, log_8_4, log_10_5, gau_9_4, gau_5_2, diffgau_1,diffgau_2, average_3, average_5, average_10, morph_opening, morph_close, edges, deconv_10_7, deconv_10_3, deconv_5_3, std_dif_1, zscore_1, watershed_lines_smoothed  
**Nucleoli**:  
gau_20_10, gau_9_4, gau_5_2, average_3, average_5, average_10, morph_opening, morph_close_small, orig_tophat_30  
**P-bodies**:  
log_5_2, log_5_3, log_6_3, log_8_4, gau_9_4, gau_5_2, diffgau_1, diffgau_2, average_3, average_5, average_10, edges, edge_log_8_4, hole_log_8_4, deconv_5_3, std_dif_1, zscore_1, watershed_lines  
**PML nuclear bodies**:  
log_5_2, log_5_3, log_6_3, log_8_4, gau_9_4, gau_5_2, diffgau_1, diffgau_2, average_3, average_5, average_10, edges, edge_log_8_4, hole_log_8_4, deconv_5_3, std_dif_1, zscore_1, watershed_lines  
**Stress granules**:  
log_8_4, log_10_5, log_20_10, average_20, simple_edges, entropy_log_8_4, entropy_log_5_2, morph_close, tophat_12, orig_tophat_30  
**Slicing speckles**:  
log_5_2, log_8_4, log_20_10, gau_20_10, gau_9_4, gau_5_2, diffgau_1, diffgau_2, average_5, average_10, average_20, edges, hole_log_8_4, morph_opening, morph_close, deconv_10_7, deconv_10_3, std_dif_1, zscore_1  
