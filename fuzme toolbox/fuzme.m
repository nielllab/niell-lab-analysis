function [U, centroid, dist, W, obj] = fuzme(nclass,data,U,phi,maxiter,distype,toldif)
% fuzzy k means
% input
% [U, centroid, dist, W, obj] = fuzme(nclass,data,U,phi,maxiter,distype,toldif)
%   nclass      = number of class
%   data        = data matrix           data(ndata,ndim)
%   U           = initial membership matrix     U(ndata,nclass)
%   phi         = fuzzy exponent        >1
%   maxiter     = maximum iterations
%   distype     = distance type:        1 = euclidean, 2 = diagonal, 3 = mahalanobis
%   toldif      = convergence tolerance
%
% output:
%   U           = new membership matrix
%   centroid    = centroid              centroid(nclass, ndim)
%   dist        = distance matrix       dist(ndata,nclass)
%   W           = distance norm matrix  W(ndim,ndim)
%   obj_fcn     = objective function
%

printing=0;
if(phi<=1), phi=1.01, end; 
    
ndata = size(data, 1);         % number of data 
ndim = size(data, 2);         % number of dimension
centroid=zeros(nclass,ndim);
dist=zeros(ndata,nclass);

% check distance type
if(distype==1)      % euclidean distance
    W=eye(ndim);
elseif(distype==2)  % diagonal
    W=eye(ndim).*cov(data);
elseif(distype==3)  % mahalanobis
    W=inv(cov(data));
end

obj=0;
uphi = U.^phi;   

for i = 1:maxiter,

    % calculate centroid
    c1=uphi'*data;
    t1=sum(uphi)';
    t1=t1(:,ones(ndim,1));
    centroid=c1./t1;

    % calculate distance of data to centroid
    if(distype==1),      % euclidean distance
        dist=distmat0(data, centroid);
    else,
        dist = sqrt(mahaldist(data, centroid, W));
    end;
    
    %   save previous iterations
    U_old=U;
    obj_old=obj;

    % calculate new membership matrix
    tmp = dist.^(-2/(phi-1));      
    t1=sum(tmp')';
    t2=t1(:,ones(nclass,1));
    U = tmp./t2;
    uphi = U.^phi;   

    % calculate objective function
    o1=(dist.^2).*uphi;
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
