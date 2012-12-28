function [fpi, nce, S, dJdphi]=fvalidity(U,W,centroid,dist,nclass,phi)
% calculate cluster validity value
% [fpi, mpe, S, dJdphi]=fvalidity(U,dist,nclass,phi)
% Input:
%   U           = membership matrix     U(ndata,nclass)
%   W           = distance norm matrix  W(ndim,ndim)
%   centroid    = cluster centre        centroid(nclass,ndim)
%   dist        = distance matrix       dist(ndata,nclass)
%   nclass      = no. class
%   phi         = fuzzy exponent        >1
% Output:
%   fpi = fuzzy performance index
%   nce = normalised classification entropy
%   S   = cluster separation 
%   dJdphi = derivative of objective function over phi dJ/dphi

ndata = size(U, 1);         % number of data 

% fuzzy performance index
F = 1/ndata*sum(sum(U.*U));
f = (nclass*F-1)/(nclass-1);
fpi = 1-f;

% mpe
H = -1/ndata*sum(sum(U.*log(U)));
nce = H/log(nclass);

% S
disc = mahaldist(centroid, centroid, W)+eye(nclass).*1/eps;
dmin=min(min(disc));
ud=(U.*U).*(dist.*dist); 
j2=sum(sum(ud));
S=j2./(ndata*dmin);

% dJ/dphi
uphi = U.^phi;   
d1=(dist.^2).*log(U).*uphi;
dJdphi = sum(sum(d1')); 