function [U, dist] = gk_all(nclass,data,phi,centroid, W)
% allocation of fuzzy k means Gustafson-Kessel algorithm
% [U, dist] = gk_fkm(nclass,data,U,phi,centroid,W)
% input
%   nclass      = number of class
%   data        = data matrix                   data(ndata,ndim)
%   U           = initial membership matrix     U(ndata,nclass)
%   phi         = fuzzy exponent        >1
%   W           = distance norm matrix          W(nclass,ndim,ndim)
%   centroid    = centroid
%
% output:
%   U           = new membership matrix
%   dist        = distance matrix           dist(ndata,nclass)
    
ndata = size(data, 1);         % number of data 
ndim = size(data, 2);         % number of dimension

dist=zeros(ndata,nclass);

% calculate distance between data & centroid
dist = GKdist(data, centroid, W);
    
% calculate membership matrix
tmp = dist.^(-1/(phi-1));      
t1=sum(tmp')';
t2=t1(:,ones(nclass,1));
U = tmp./t2;

dist=sqrt(dist);