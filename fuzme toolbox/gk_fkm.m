function [U, centroid, dist, W, F, obj] = gk_fkm(nclass,data,U,phi,maxiter,toldif)
% fuzzy k means with Gustafson-Kessel algorithm
% [U, centroid, dist, W, F, obj] = gk_fkm(nclass,data,U,phi,maxiter,toldif)
% input
%   nclass      = number of class
%   data        = data matrix                   data(ndata,ndim)
%   U           = initial membership matrix     U(ndata,nclass)
%   phi         = fuzzy exponent        >1
%   maxiter     = maximum iterations
%   toldif      = convergence tolerance
%
% output:
%   U           = new membership matrix
%   centroid    = centroid                  centroid(nclass, ndim)
%   dist        = distance matrix           dist(ndata,nclass)
%   W           = distance norm matrix      W(nclass,ndim,ndim)
%   F           = fuzzy covariance matrix   F(nclass,ndim,ndim)
%   obj_fcn     = objective function

printing=1;
if(phi<=1), phi=1.01, end; 
    
ndata = size(data, 1);         % number of data 
ndim = size(data, 2);         % number of dimension
centroid=zeros(nclass,ndim);
dist=zeros(ndata,nclass);

obj=0;
uphi = U.^phi;   

for i = 1:maxiter,

    % calculate centroid
    c1=uphi'*data;
    t1=sum(uphi)';
    t1=t1(:,ones(ndim,1));
    centroid=c1./t1;

    % calculate fuzzy covariance matrix
    for k=1:nclass,
        ufi=U(:,k).^phi;
        c1=data-repmat(centroid(k,:),ndata,1);
        c2=repmat(ufi,1,ndim).*c1;
        c3=c2'*c1;
        c3=c3./sum(ufi);
        S=(det(c3)).^(1/2)*inv(c3);
        F{k}=c3;
        W{k}=S;
    end
    
    % calculate distance between data & centroid
    dist = GKdist(data, centroid, W);
    
    %   save previous iterations result
    U_old=U;
    obj_old=obj;

    % calculate new membership matrix
    tmp = dist.^(-1/(phi-1));      
    t1=sum(tmp')';
    t2=t1(:,ones(nclass,1));
    U = tmp./t2;
    uphi = U.^phi;   

    % calculate objective function
    o1=(dist).*uphi;
    obj = sum(sum(o1')); 
    
	% check for convergence
    dif=(obj_old-obj);
    difU=sqrt((U - U_old).*(U - U_old));
    Udif=sum(sum(difU));
    if printing==1,
 	    fprintf('Iteration = %d, obj. fcn = %f.  diff = %f\n', i, obj, Udif);
    end
    if and(dif<toldif,Udif < toldif), break; end,
    
end

dist=sqrt(dist);