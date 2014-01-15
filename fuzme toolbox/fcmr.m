function [U, centroid, dist, b, obj] = fcmr(nclass,data,y,U,phi,maxiter,toldif)
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

if(phi<=1), phi=1.01, end; 
    
ndata = size(data, 1);         % number of data 
ndim = size(data, 2);         % number of dimension
centroid=zeros(nclass,ndim);
dist=zeros(ndata,nclass);

W=inv(cov(data));

obj=0;
%x=[ones(ndata,1), data];
%nd = ndim+1;       % number of dimension
x=[data];
nd=ndim;

b=zeros(nd,nclass);
yp=zeros(ndata,nclass);

for i = 1:maxiter,

    %   save previous iterations
    U_old=U;
    obj_old=obj;
        
    for k=1:nclass
    % weighted LS
        wei(:,1)=U(:,k); 
        A=wei(:,ones(nd,1)).*x;
        y1=wei.*y;
        b(:,k)=A\y1;
        yp(:,k)=x*b(:,k);
        dist(:,k)=(yp(:,k)-y).^2;
    end
    
    % calculate new membership matrix
    tmp = dist.^(-1/(phi-1));      
    t1=sum(tmp')';
    t2=t1(:,ones(nclass,1));
    U = tmp./t2;

    
    % calculate objective function
    uphi = U.^phi;   
    o1=(dist.^2).*uphi;
    obj = sum(sum(o1')); 
    
	% check for convergence
    dif=(obj_old-obj);
    difU=sqrt((U - U_old).*(U - U_old));
    Udif=sum(sum(difU));
 	
    fprintf('Iteration = %d, obj. fcn = %f.  diff = %f\n', i, obj, Udif);

    if and(dif<toldif,Udif < toldif), break; end,
    
end


    % calculate centroid
    uphi = U.^phi;   
    c1=uphi'*data;
    t1=sum(uphi)';
    t1=t1(:,ones(ndim,1));
    centroid=c1./t1;
%    npar=nclass*ndim;
%    x0=reshape(centroid,npar,1);
%    pars={U,W,data,phi,ndata,nclass,ndim};
%    opts=[1 1e-10 1e-12 500 1e-6]; 
    %[xx] = SMarquardt(@fciter,pars, x0, opts);
%    xx=lsqnonlin(@fciter,x0,[],[],[],U,W,data,phi,ndata,nclass,ndim)

%    centroid=reshape(xx,nclass,ndim);
    % calculate U according to FCM
%    dis2 = mahaldist(data, centroid, W);
%    tmp = dis2.^(-1/(phi-1));      
%    t1=sum(tmp')';
%    t2=t1(:,ones(nclass,1));
%    U = tmp./t2;
