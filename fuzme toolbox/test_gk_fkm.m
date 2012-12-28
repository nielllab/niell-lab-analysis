% test gk_fkm on cross data
%%%%%%%%%%
clear all
% Gustafson's cross
data = [
   -9.7500   -0.1500
   -6.4400    0.3400
   -4.6900   -0.3000
   -2.0400    0.3700
   -1.2400    0.4500
    0.3300   -0.0800
    5.0400   -0.2100
    5.8600   -0.2500
    7.5400    0.1600
    7.6700    0.2400
   -0.3000   -8.0700
    0.1300   -7.1300
   -0.3700   -5.1800
    0.0300   -3.3300
    0.3500   -2.6300
    0.2300   -2.6800
   -0.0500   -2.0000
    0.4100    0.3700
    0.6900    4.7500
    0.7400    8.8700];

nclass=2;                       % number of class
phi=2;                          % fuzzy exponent

maxiter=200;
toldif=0.0000001;

ndata = size(data, 1);         % number of data 
% initialise random 
Uinit= initmember(0.1,nclass,ndata);
[U, centroid, dist, W, F, obj] = gk_fkm(nclass,data,Uinit,phi,maxiter,toldif);
centroid
F{1}
F{2}

% plot the result
maxU=max(U');
figure(1);
hold on;
k=1
    index = find(U(:, k) == maxU');
    cluster = data(index', :);
    plot(cluster(:, 1), cluster(:, 2),'ro',centroid(1,1),centroid(2,1),'r+');
k=2
    index = find(U(:, k) == maxU');
    cluster = data(index', :);
    plot(cluster(:, 1), cluster(:, 2),'bo',centroid(1,2),centroid(2,2),'b+');
    
% centroid should be like:
%    0.1875   -1.6712
%    0.2399    0.0681
%
% F{1}
%    0.1240    1.4643
%    1.4643   23.8535
%
% F{2}
%   37.8132   -0.0199
%   -0.0199    0.0782
