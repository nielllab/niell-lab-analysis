% FuzME
% Fuzzy k-means matlab code
% Budiman (Dec. 2003)
% budiman@acss.usyd.edu.au
%
% Contents:
% run_fuzme.m       : main command to run fuzzy k-means 
% fuzme.m           : fuzzy k-means algorithm
% fuzall.m          : allocate data to existing class centroid
% run_fuzme.m       : main command to run fuzzy k-means with extragrades FkME
% fkme.m            : fuzzy k-means with extragrades
% fkme_all.m        : allocate data to existing class of FkME
% gustafson.m       : fuzzy k means with Gustafson-Kessel (GK) algorithm
% gk_all.m          : allocate data to existing class of GK
%
% flda.m            : fuzzy linear discriminant analysis
%
% fvalidity.m       : cluster performance measure
% confusion.m       : calculate confusion index
% 
% initmember.m      : calculate initial guess membership
% distmat0.m        : calculate Euclidean distance, code from Peter J. Acklam
% mahaldist.m       : calculate Mahalanobis distance, code from Peter J. Acklam
% GKdist.m          : calculate distance of GK algorithm, based on code from Peter J. Acklam
% fkme_obj          : objective function to minimise for FkME
%
% Example
% test_fuzme.m      : clustering data of irises
% irises.txt        : text data of irises
% test_gk_fkm.m     : test for GK algorithm
%
% For thery see:
%       http://www.usyd.edu.au/su/agric/acpa/fkme/FuzME_Theory.pdf
%
% Acknowldegement:
% distance calculation matlab files are from Peter J. Acklam