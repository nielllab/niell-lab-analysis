function [U, Ue, dist, obj] = fkme_all(nclass,data,centroid,W,phi,alfa,distype)
% allocate fuzzy k means with extragrades
% input
% [U, centroid, dist, W, obj] = fuzme(nclass,data,U,phi,maxiter,distype,toldif)
%   nclass      = number of class
%   data        = data matrix           data(ndata,ndim)
%   centroid    = centroid of clasees
%   phi         = fuzzy exponent        >1
%   alfa        = extragrade parameter
%   maxiter     = maximum iterations
%   distype     = distance type:        1 = euclidean, 2 = diagonal, 3 = mahalanobis
%
% output:
%   U           = membership matrix
%   Ue          = extragrade membership
%   dist        = distance matrix       dist(ndata,nclass)
%   obj         = objective function
%

ndata = size(data, 1);         % number of data 
ndim = size(data, 2);         % number of dimension
dist=zeros(ndata,nclass);

% calculate distance of data to centroid
    if(distype==1),      % euclidean distance
        dist=distmat0(data, centroid);
    else,
        dist = sqrt(mahaldist(data, centroid, W));
    end;

a1=(1-alfa)/alfa;
% calculate membership matrix
    tmp = dist.^(-2/(phi-1)); 
    tm2 = dist.^(-2);
    s2=(a1.*sum(tm2')').^(-1/(phi-1));
    
    t1=sum(tmp')';
    t2=repmat(t1,1,nclass)+repmat(s2,1,nclass);
    U = tmp./t2;
    Ue=ones(ndata,1)-sum(U')';
    uphi = U.^phi;   
    uephi= Ue.^phi;
    
% calculate objective function
    o1=(dist.^2).*uphi;
    d2=dist.^(-2);
    o2=uephi.*sum(d2')';
    obj = alfa*sum(sum(o1'))+(1-alfa)*sum(o2); 