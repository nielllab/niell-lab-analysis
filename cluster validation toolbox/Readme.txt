Help manual of Cluster Validity Analysis Platform (CVAP) (Version 3.7)


You are welcome to give some comments, which are helpful to the improvements of this tool.
Citation information: Kaijun Wang, Baijie Wang, and Liuqing Peng. CVAP: Validation for Cluster Analyses. Data Science Journal, Vol. 8, pp.88-93, 2009.
E-mail: wangkjun@yahoo.com

Note 1: To start CVAP, type "mainCVAP" in Matlab command window (Please do not click "mainCVAP.fig", which is a figure file). The CVAP was tested under Matlab 7.2.

Note 2: Statistics Toolbox of Matlab needs to be installed, since it contains routines such as K-means and Silhouette index.

Note 3: The authors of the contributed files/codes hold the copyrights (see "Notice_contributedFiles.txt"):
(a) the PAM codes are from LIBRA (http://wis.kuleuven.be/stat/robust/LIBRA.html)
(b) the hierarchical clustering codes (the original file is altered by adding some codes) are 
     from Clustering Toolbox by David Corney (http://www.dcorney.com/)
(c) the SOM codes are from Netlab neural network software (www.ncrg.aston.ac.uk/netlab/), 
(d) the codes of Affinity Propagation clustering are from "Brendan J. Frey and Delbert Dueck. Clustering by Passing Messages Between Data Points. Science 315, 972 (2007)"
(e) and the xlim2.m (the original file is altered by adding some codes) is from Matlab

------------------------------------------------------------------------------------------
(1) Contents of CVAP
Purpose
    The evaluation of clustering results is important in cluster analysis. It is the validity indices that are usually used to evaluate clustering results. To supply a visual analysis environment to facilitate cluster validation, we developed cluster validity analysis platform (CVAP) and implemented most of the validity indices listed below. The CVAP based on a graphical user interface (GUI) is a convenient working environment for cluster validation, integrating 18 validity indices and several clustering algorithms (PAM, K-means, hierarchical, SOM and etc.).

External and internal validity indices included
i)  External validity indices are the measures of the agreement between two partitions, one of which is usually a known/golden partition, e.g. true class labels, and another is from the clustering procedure. 
    Rand index                 (codes from Clustering Toolbox by David Corney)
    Adjusted Rand index   (codes from Clustering Toolbox by David Corney)
    Jaccard index
    Fowlkes-Mallows (FM) index
ii) Internal validity indices evaluate clustering results by using only features and information inherent in a dataset. They are usually used in the case that true solutions are unknown.
    Silhouette index          (an inner function of Matlab)
    Davies-Bouldin
    Calinski-Harabasz
    Dunn index
    R-squared index
    Hubert-Levin (C-index)
    Krzanowski-Lai index
    Hartigan index
    Root-mean-square standard deviation (RMSSTD) index
    Semi-partial R-squared (SPR) index
    Distance between two clusters (CD) index
    weighted inter-intra index
    Homogeneity index
    Separation index

The function parts of CVAP:
Part 1: Data/solution loading and result saving. 
Part 2: Clustering process. 
Part 3: Cluster validation process.
Part 4: Result plotting and visualization.


(2) Part 1: Data/solution loading and result saving
    Data/solution loading: In CVAP GUI, a data file like "yourdata.txt" or demo data sets may be loaded from File menu: click "Load Data File" to find and load a data file, lick "Load Demo Data" for demo using of CVAP (see (6) Demo using of CVAP), or click "Load Solution File" to load a solution file with class labels from other clustering algorithms.
    
    Data format: input data file is tab delimited text file with numeric tabular data (e.g. rows denote data points and columns are dimensions), and all the data should be numeric values and without characters and missing values. Such tab-delimited text files can be created and exported by such as Microsoft Excel.
    For example, Let there be four 4-dimensional samples (data points):
sample1:	(5.1, 3.5, 1.4, 0.2)
sample2:	(4.9, 3.0, 1.4, 0.2)
sample3:	(4.7, 3.2, 1.3, 0.2)
sample4:	(1.7, 1.2, 1.3, 0.2)
    Then, the data file "mydata.txt" looks like:
5.1	3.5	1.4	0.2
4.9	3	1.4	0.2
4.7	3.2	1.3	0.2
1.7	1.2	1.3	0.2
and  you should choose "All columns/rows are data" in the pop-up menu "Position".  

    When true class labels (class 1 and class 2) are available, put them in the first column of the data file:
2	5.1	3.5	1.4	0.2
2	4.9	3	1.4	0.2
1	4.7	3.2	1.3	0.2
1	1.7	1.2	1.3	0.2
and  you should choose "Class labels in first column/row" (default) in the pop-up menu "Position", and not change "Clustering Rows" (default) in Option menu.
    
    Solution-file format: it is tab delimited text file with numeric tabular data; its first row consists of numbers of clusters, and each column stores class labels at a number of clusters from a clustering procedure (refer to sample solution file included); if true class labels are available, they must be in the first column and the first element of this column should be 0 to indicate true lables. How to use solution files, refer to Note 1 of following (4) Part 3: Cluster validation.
    For example, for the above example the solution file "mysolution.txt" looks like:
1	2	3	4
1	2	2	1
1	2	2	2
1	2	3	3
1	1	1	4
The first row shows the number of clusters, e.g., 2 means that the data set is divided into 2 clusters. Each column (except the first row) gives class labels of every data point.
    When true class labels are available, they are in the first column and the first element of this column is 0 (instead of 1) to indicate true lables. The file looks like:
0	2	3	4
2	2	2	1
2	2	2	2
1	2	3	3
1	1	1	4


    Result saving: After running a clustering algorithm, one may save the output class labels in a Matlab ".mat" data file by clicking "Save Results" in File menu. Similarly, the values of validity indices and the output class labels may be saved after running validation.


(3) Part 2: Clustering process
The related items for a clustering process:
  # The PAM (partitioning around medoids), K-means, hierarchical and SOM (Self-Organizing Map) clustering algorithms are included, and one may select one of them in the pop-up menu "Algorithm".
  # "Clustering Rows" (default) or "Clustering Columns" of tabular data is optional in Option menu, and you can change this choice only before data loading.
  # Data normalization (namely, each row/column <matching clustering-rows/columns option> is scaled to zero mean and unit variance.) or standardization (values are scaled within [0 1]) is optional in Option menu, click or unclick it if applicable or not. The default is none.
  # Specifying position of known class labels from the pop-up menu "Position" 
  "Class labels in first column/row" (default), "Class labels in end column/row" or "All columns/rows are data" is optional.
  # Selecting a similarity metric from the pop-up menu "Measure"
  The Euclidean distance and Pearson correlation coefficient are supported. 
  # Specifying single or multiple (default) clustering from the pop-up menu "Times"
  the single clustering runs only under a given number of clusters (NC), the multiple clustering runs in a range of NC. 
  # Giving (a range of) number of clusters (NC)
  set a NC in the first NC box for single clustering or set NCs in two NC boxes for multiple clustering, since pre-assigned NCs are needed for PAM, K-means and SOM. To calculate validity indices under different NCs, the NCs need to be given for hierarchical algorithm so that the hierarchical tree at different levels can be cut to obtain corresponding clusters.

  Once all these items are specified, one can perform the clustering by pressing button "Run Clustering".

Note 1: Initialization of K-means (from Statistics Toolbox of Matlab) is to select K centroids from data at random, for other choices, refer to valid_clusteringAlgs.m (row 31) of CVAP and kmeans.m of Matlab.

Note 2: Initialization and training methods of SOM are similar to the ones in "demsom1.m" of Netlab except adaptive iterations and initial neighborhood size, see "som_netlab.m" of CVAP for detail.

Note 3: A new clustering algorithm may be added in last 4 rows of "valid_clusteringAlgs.m" if you are familiar to code making, and then change "ialg >= 8" to "ialg > 8" in row 2 of "valid_runclustering.m". It is to select "New algorithm" in pop-up menu "Algorithm" when using this added algorithm.

(4) Part 3: Cluster validation
    Once a clustering algorithm divides a data set into k (k=2,3,4,¡­) clusters, the validity indices may be selected from the pop-up menu "Validity indices" and used to evaluate the clustering solutions. The following two tasks are common in the literature.

i)  Comparison/evaluation of clustering algorithms
    After clicking "Multi-Algorithm Validation" in the Validation menu, choose and run a clustering algorithm, and then specify a validity index and run the validation. When the validation processes under a same index are carried out for several clustering algorithms, the index values corresponding to different clustering algorithms are displayed in the plotting window for further analyses. These indices were used in the literature:
  # Rand, Adjusted Rand, Jaccard and Fowlkes-Mallows indices
     higher the score, better the solution [Dudoit et al. 2002; Halkidi 2001; Sharan et al. 2003].
  # Silhouette index (overall average silhouette)
     a larger Silhouette value indicates a better quality of a clustering result [Chen et al. 2002].
  # Davies-Bouldin index
     a low value indicates good cluster structures [Kasturi et al. 2003; Bolshakova et al. 2003].
  # Calinski-Harabasz index
     it is the pseudo F statistic, which evaluates the clustering solution by looking at how similar the objects are within each cluster and how well the objects of different clusters are separated [Zhao et al. 2005; Shu et al. 2003].
  # Dunn index
     large values indicate the presence of compact and well-separated clusters [Bolshakova et al. 2003; Halkidi et al. 2001].
  # R-squared index
     large R-squared statistic indicates large difference between clusters [Halkidi et al. 2001; Shu et al. 2003].
  # Homogeneity & Separation indices
     improving Homogeneity & Separation suggests an improvement in clustering results [Sharan et al. 2003; Chen et al. 2002].

Note: For some examples of using these indices, refer to the following references:
> G. Chen, S. A. Jaradat, N. Banerjee, T. S. Tanaka, M. S. H. Ko and M. Q. Zhang. Evaluation and Comparison of Clustering Algorithms in Anglyzing ES Cell Gene Expression Data. Statistica Sinica 12(2002), 241-262 
> M. Halkidi, Y. Batistakis, M. Vazirgiannis. On Clustering Validation Techniques. Intelligent Information Systems Journal, 17(2-3): 107-145, 2001
> J. Kasturi, R. Acharya and M. Ramanathan. An information theoretic approach for analyzing temporal patterns of gene expression. Bioinformatics, Vol. 19 no. 4 2003, pages 449¨C458
> R. Sharan, A. Maron-Katz and R. Shamir. CLICK and EXPANDER: A System for Clustering and Visualizing Gene Expression Data. Bioinformatics 19, pp. 1787-1799, 2003
> G. Shu, B. Zeng, Y. P. Chen, O. H. Smith. Performance assessment of kernel density clustering for gene expression profile data. Comparative and Functional Genomics, 1 June 2003, vol. 4, no. 3,  pp. 287-299(13)
> Y. Zhao and G. Karypis. Data Clustering in Life Sciences. Molecular Biotechnology, 31(1), pp. 55¡ª80, 2005. 

ii) Estimation of the optimal number of clusters (NC)
    After clicking "Estimate Number of Clusters" in the Validation menu, specify and run a clustering algorithm, and then select a validity index and run the validation. The index values across NCs will appear in the plotting window, and the optimal NC is indicated by a square symbol (except the RMSSTD group and R-squared indices). The following indices were used in the literature:
  # Silhouette index (overall average silhouette)
     the largest silhouette indicates the optimal NC [Dudoit et al. 2002; Bolshakova et al. 2003].
  # Davies-Bouldin
     minimum value determines the optimal number of clusters [Bolshakova et al. 2003; Dimitriadou et al. 2002].
  # Calinski-Harabasz
     maximum value indicates optimal NC [Dudoit et al. 2002].
  # Dunn index
     maximum value indicates optimal NC [Bolshakova et al. 2003; Halkidi et al. 2001].
  # C index (Hubert-Levin)
     minimal C-index indicates optimal NC [Bolshakova et al. 2003; Bolshakova et al. 2006; Dimitriadou et al. 2002].
  # Krzanowski-Lai index
     maximum value indicates optimal NC [Dudoit et al. 2002].
  # Hartigan (Ha) index
     the estimated number of clusters is the smallest k ¡Ý 1 such that Ha ¡Ü 10 [Dudoit et al. 2002]. There will be no estimation if this condition is not satisfied. However, an elbow might indicate an optimal NC.
  # RMSSTD group: Root-mean-square standard deviation (RMSSTD), R-squared, Semi-partial R-squared (SPR), Distance between two clusters (CD) indices.
     the steepest knee indicates optimal NC, i.e., the greater jump of these indices¡¯ values from larger to smaller NC [Halkidi et al. 2001; Frossyniotis et al. 2005].
  # weighted inter-intra index
     search forward (k=2,3,4,...) and stop at the first down-tick of the index, which indicates optimal NC [Strehl 2002].
   
Note 1: Cluster validation process on a solution file is as same as on a data file except that it need not run a clustering process. For cluster validation with internal indices, the corresponding data file needs to be loaded at the same time.

Note 2: RMSSTD index group is designed specially and recommended to estimate the optimal NC for hierarchical algorithms.

Note 3: The inter-cluster distance for calculating the Dunn¡¯s and Davies-Bouldin indices in CVAP is the centroid linkage distance, the distance between centers of two clusters; and the intra-cluster distance is the centroid diameter distance, the average distance between all of the samples and the cluster¡¯s center.

Note 4: The calculation of C-index, weighted inter-intra index and Homogeneity & Separation indices for large datasets could be very time-consuming.

Note 5: For the detail of these indices, refer to the following references:
> N. Bolshakova, F. Azuaje. 2003. Cluster validation techniques for genome expression data, Signal Processing. V.83. N4, P.825-833.
> N. Bolshakova, F. Azuaje. Estimating the number of clusters in DNA microarray data. Methods of Information in Medicine, Feb 2006.
> E. Dimitriadou, S. Dolnicar, A. Weingessel. An examination of indexes for determining the Number of Cluster in binary data sets. Psychometrika, 67(1):137-160, 2002.
> S. Dudoit, J. Fridlyand. (2002). A prediction-based resampling method for estimating the number of clusters in a dataset. Genome Biology. 3(7): 0036.1-21.
> D. S. Frossyniotis, C. Pateritsas, A. Stafylopatis. A multi-clustering fusion scheme for data partitioning. Int. J. Neural Syst. 15(5): 391-401 (2005)
> M. Halkidi, Y. Batistakis, M. Vazirgiannis, ¡°On Clustering Validation Techniques¡±, Intelligent Information Systems Journal, 17(2-3): 107-145, 2001
> A. Strehl. Relationship-based Clustering and Cluster Ensembles for High-dimensional Data Mining. Ph.D thesis, The University of Texas at Austin, May 2002.

iii) Other functions
    "Error rate" in Tool menu: it is the percentage of error class labels of a clustering solution compared with true labels (if any) under true NC. The error rate might be inaccurate if a clustering solution under true NC has an error rate > 20%, since it is hard to count matching labels properly in complex cases. Hence, this function is only for reference, and a better tool is external indices.

(5) Part 4: Result printing, plotting and visualization.
  Clustering solutions (class labels) may be printed in Matlab command window by clicking "Print Class Labels" in Tool menu after the clustering is carried out, while true class labels (if any) are printed when "0" is given as input.

  In a validation process, values of validity indices are printed in Matlab command window and shown in Plotting window, and you may redraw them for a better display by pressing button "Redrawing" or clear them by "Clear Plotting".

  When one needs copy, paste or save the result plotting, click "Plot in New Figure" in View menu (unclick to cancel), and then the plotting and redrawing will work in a new figure.

  When "All validity indices" in pop-up menu "Validity indices" is selected (not applicable to "Multi-Algorithm Validation") in a validation process, values of all validity indices are shown in several new figures or in the plotting window.

  When "Homogeneity-Separation" indices are plotted together, the two-row title (if any) above the plotting window indicates which index (Homogeneity or Separation) is in upper position.

  If you are interested in Principal component analysis (PCA) of data, click "Plot Data by PCA" or "Plot Class Labels" in the View menu to see data or class labels in space of first two PCs (for reference only). The data will be plotted directly when there are two dimensional data. When clicking "Plot Class Labels" before running clustering, true labels (if any) instead of clustering solutions will be plotted.

  To see a hierarchical dendrogram of data, click "Plot Dendrogram" in the View menu, and there appears a dendrogram corresponding to single linkage, complete linkage or centroid linkage when 'Hierarchical-single', 'Hierarchical-complete' or 'Hierarchical-centroid' in pop-up menu "Algorithm" is selected, otherwise the dendrogram corresponds to average linkage. (This function is from the Statistics Toolbox of Matlab.)

(6) Demo using of CVAP
Data set "leuk72_3k" has 72 samples in 3 classes, and true class labels are in the first column.

Type "mainCVAP" in Matlab command window to start CVAP. The clustering process should be performed before a cluster validation process, but no clustering process is needed for a solution file. 
<1> Clustering process
  # Loading a data file
  Click "Load Demo Data" in File menu bar, and then "Running state: Default data set ... is loaded, ..." appears in the running state window.
  We may select demo data set "leuk72_3k.txt" in pop-up menu "Data set", or click "Load Data File" to load a new data file.
  # Specifying "Clustering Rows" (default) or "Clustering Columns" in Option menu
  no action needed here, since the default is applicable.
  # Specifying data normalization (optional) in Option menu (click it if applicable and unclick it to cancel)
  no action needed here, since the data has been normalized. 
  # Specifying the position of known class labels from the pop-up menu "Position" 
  no action needed here, since default "Class labels in first column/row" is applicable.
  # Selecting a similarity metric from the pop-up menu "Measure"
  no action needed here, since default Pearson correlation is applicable.
  # Choosing a clustering algorithm from the pop-up menu "Algorithm"
  no action needed here, since default PAM is applicable.
  # Specifying single or multiple (default) clustering from the pop-up menu "Times"
  no action needed here, since default multiple clustering is applicable.
  # Giving the (range of) number of clusters
  "2" & "10" in Edit boxes is given as a range of number of clusters for running multiple clustering. The "10" will be neglected if "Single Clustering" is selected.
  # Running clustering
  press button "Run Clustering", the following appears in the Running state window:
  "Clustering state: PAM clustering on leuk72_3k.txt at k = 2 - 10 runs over, ..."

  # Printing and saving clustering solutions (optional)
  class labels at k = 3 are printed in Matlab command window when "Print Class Lables" in Tool menu is clicked and the 3 is given as input.
  class labels at k = 2 - 10 are saved in "test.mat" when "Save Results" in File menu is clicked and file name "test.mat" is given. You may check them by using "load test.mat" in Matlab command window.
  # Plotting data and clustering solutions (for reference only)
  Click "Plot Data by PCA" or "Plot Class Labels" in View menu to see data or solutions in space of first two PCs.

<2> Validation process
  After running the clustering process stated in <1>, turn to
i) Clustering quality evaluation for a clustering algorithm
  # Specifying Clustering Quality Evaluation (default) or Estimate Number of Clusters in Validation menu
  no action needed here, the default is applicable.
  # Choosing Single-algorithm (default) or Multi-algorithm validation in Validation menu
  no action needed here, the default is applicable.
  # Selecting a validity index from pop-up menu "Validity Indices"
  Silhouette index is selected here.
  # Running cluster validation
  press button "Run Validation", output index values with regard to the number of clusters are displayed in the plotting window. The largest Silhouette value at k = 3 indicates the best clustering solution at k = 3.
  # Redrawing index values or clear plotting (optional)
  press button "Redrawing" or "Clear Plotting". Redrawing works in a new figure if "Plot in New Figure" in View menu is clicked first (unclick to cancel).
  # Error rate of a clustering solution under true NC
  it is 4.1667% in Running state window when clicking "Error rate" in Tool menu.

ii) Estimation of the number of clusters
  # Be sure that a clustering process is multiple clustering in a proper range of NCs. here we skip it, since PAM result exists.
  # Click "Estimate Number of Clusters"  (applicable to single-algorithm validation only) in Validation menu
  # Select Silhouette index (or All validity indices) from pop-up menu "Validity Indices"
  # Press button "Run Validation"
  index values across the number of clusters are displayed in the plotting window, and the optimal NC (k = 3) is indicated by a square symbol.
  # Redrawing in a new figure
  click "Plot in New Figure" in View menu, and then press button "Redrawing".

iii) Comparison analysis between several algorithms
  # Click "Multi-Algorithm Validation" in Validation menu
  # Run the clustering process stated in <1>. Here we skip it, since PAM result exists.
  # Select Silhouette index from pop-up menu "Validity Indices"
  # Press button "Clear Validation" if starting a new work. Here we skip it.
  # Press button "Run Validation"
  index values across the number of clusters are displayed for PAM algorithm. 

  # Choose a new algorithm such as k-means from pop-up menu "Algorithm", but load a solution file which stores clustering solutions from an external clustering algorithm/procedure.
  # Press button "Run Clustering", but skip this step for a solution file.
  # Press button "Run Validation".
  Silhouette values across the number of clusters are displayed for PAM and k-means.
  # For redrawing in a new figure, click "Plot in New Figure" in View menu first.
  # Repeat last four steps for a new algorithm.

  Now the Silhouette values corresponding to PAM, k-means and etc. are ready for further analyses.

---------------------------------------------------------------------------------------------------
This free software is distributed under BSD license.  (see license.txt)
Copyright (C) 2006-2008.
Last modified: July 25, 2009
