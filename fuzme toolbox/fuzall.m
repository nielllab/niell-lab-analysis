function [U,dist,obj] = fuzall(data,phi,centroid,distype,W)
% fuzzy k means allocation
% calculate membership U & distance matrix, given centroid & phi
% [U,dist,obj] = fuzall(data,phi,centroid,distype,W)
%
% input
%   data        = data matrix           data(ndata,ndim)
%   phi         = fuzzy exponent        >1
%   centroid    = centroid              centroid(nclass, ndim)
%   distype     = distance type:        1 = euclidean, 2 = diagonal, 3 = mahalanobis
%   W           = distance norm matrix  output from fuzme
%
% output:
%   U           = membership matrix
%   dist        = distance matrix       dist(ndata,nclass)
%   obj         = objective function
%
% Budiman (2003)

nclass = size(centroid,1);     % number of class
ndata = size(data, 1);         % number of data 
ndim = size(data, 2);         % number of dimension
U =zeros(ndata,nclass);


obj=0;

% calculate distance of data to centroid
if(distype==1),      % euclidean distance
    dist=distmat0(data, centroid);
else,
    dist = sqrt(mahaldist(data, centroid, W));
end;

    
% calculate new membership matrix
tmp = dist.^(-2/(phi-1));      
t1=sum(tmp')';
t2=t1(:,ones(nclass,1));
U = tmp./t2;

% calculate objective function
uphi = U.^phi;   
o1=(dist.^2).*uphi;
obj = sum(sum(o1')); 
    
