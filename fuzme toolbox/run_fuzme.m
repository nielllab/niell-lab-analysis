function [U, cbest, dist, W, obj] = run_fuzme(nclass,data,phi,maxiter,distype,toldif,scatter,ntry)
% Fuzzy k means clustering
% cluster data into nclass with exponent phi
% [U, cbest, dist, W, obj] = run_fuzme(nclass,data,phi,maxiter,distype,toldif,scatter,ntry)
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
%   U           = new membership matrix
%   cbest       = centroid              cbest(nclass, ndim)
%   dist        = distance matrix       dist(ndata,nclass)
%   W           = distance norm matrix  output from fuzme
%   obj         = objective function
%
% Budiman (2003)

ndata = size(data, 1);         % number of data 
ndim = size(data, 2);         % number of dimension
ftry = zeros(ntry, 1);
fmin=1.e30;
cbest=ones(nclass,ndim);

for i=1:ntry,

  	%fprintf('Try no = %d\n', i);

    Uinit= initmember(scatter,nclass,ndata);
    [U,centroid,dist,W,obj] = fuzme(nclass,data,Uinit,phi,maxiter,distype,toldif);
    if obj<fmin,
        fmin=obj;
        cbest=centroid;
    end;
    ftry(i)=obj;
    
end;

[U,dist,obj] = fuzall(data,phi,cbest,distype,W);
