% sample command to run Fuzme
clear all;
load irises.txt               % Load data from a text file
data = irises;                % data to be clustered

nclass=3;                   % number of class
phi=2;                      % fuzzy exponent >1
maxiter=300;                % maximum iterations
toldif=0.000001;            % convergence criterion
distype=3;                  %   distance type:        1 = euclidean, 2 = diagonal, 3 = mahalanobis
scatter=0.2;                % scatter around initial membership
ntry=10;                    % number of trial to choose an optimal solution

% run fuzme
[U, centroid, dist, W, obj] = run_fuzme(nclass,data,phi,maxiter,distype,toldif,scatter,ntry);

% output:
%   U           = membership matrix
%   centroid    = centroid              centroid(nclass, ndim)
%   dist        = distance matrix       dist(ndata,nclass)
%   W           = distance norm matrix
%   obj         = objective function

% calculate validity 
[fpi mpe S djdphi]=fvalidity(U,W,centroid,dist,nclass,phi);

% calculate confusion index
ci = confusion(nclass,data,U);

% perform fuzzy linear discriminant analysis
scaling=2;
[F,WC,BC,wilks,V,e,vexp,z,zc]=flda(data,nclass,U,centroid,phi,scaling);

% To test the allocate function 
% to allocate say new data into existing centroid 
[U, dist, obj] = fuzall(data,phi,centroid,distype,W);