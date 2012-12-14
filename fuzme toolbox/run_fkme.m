function [U, Ue, centroid, dist, W, alfa, obj] = run_fkme(nclass,data,phi,maxiter,distype,toldif,scatter,ntry)
% Fuzzy k means with extragrades
% [U, centroid, dist, W, obj] = run_fkme(nclass,data,phi,maxiter,distype,toldif,scatter,ntry)
% input
%   nclass      = number of class
%   data        = data matrix           data(ndata,ndim)
%   phi         = fuzzy exponent        >1
%   maxiter     = maximum iterations
%   distype     = distance type:        1 = euclidean, 2 = diagonal, 3 = mahalanobis
%   toldif      = convergence tolerance
%   scatter     = scatter around initial membership ~0.1
%   ntry        = number of trial to choose optimal solution
%
% output:
%   U           = membership matrix
%   Ue          = extragrade membership
%   centroid    = centroid              centroid(nclass, ndim)
%   dist        = distance matrix       dist(ndata,nclass)
%   W           = distance norm matrix  
%   alfa        = extragrade parameter
%   obj         = objective function
%
% Budiman (2003)

ndata = size(data, 1);         % number of data 
ndim = size(data, 2);         % number of dimension
Uereq=1/(nclass+1);             % required mean extragrade membership
atry=1/(nclass+1);              % trial value for alfa
%   find best alfa
[alfa,falfa] = fzero(@fkme_obj,atry,[],Uereq,nclass,data,phi,maxiter,distype,toldif);
%   find appropriate centroid & membership function
Uinit= initmember(scatter,nclass,ndata);
[U, Ue, centroid, dist, W, obj] = fkme(nclass,data,Uinit,phi,alfa,maxiter,distype,toldif);    
    
