
clear all;
n1 = 200; n2 = 100;
n = n1+n2;

x = zeros(n,2);

% generate Gaussian samples of class 1
randn('seed',13);   % set seed
x(1:n1,:) = ones(n1,1)*[3 3]+(ones(n1,1)*[1.7 0.5]).*randn(n1,2);

% generate Gaussian samples of class 2
randn('seed',5431); % set seed
x(n1+1:n,:) = ones(n2,1)*[-4 -4]+(ones(n2,1)*[1 1.3]).*randn(n2,2);

data=x;

nclass=2;
maxiter=200;
phi=1.5;
toldif=0.0001;
ndata = size(data, 1);         % number of data 
Uinit= initmember(0.1,nclass,ndata);
[U, centroid, dist, W, F, obj] = gustafson(nclass,data,Uinit,phi,maxiter,toldif)

W{1}
W{2}
    
% xc should be like
%
%    3.2590   -3.9495
%    2.9463   -4.1688
%
% F1
%
%    2.7873    0.0111     ---> 1.67 and 0.48 as estimates of std's
%    0.0111    0.2316
%
% F2
%    1.0437   -0.0901     ---> 1.02 and 1.38 as estimates of std's
%   -0.0901    1.9201


t = 1:ndata;
plot(t,U(:,1),'y',t,U(:,1),'c');
title('membership functions');



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

nclass=2;
maxiter=500;
phi=2;
toldif=0.0000001;
ndata = size(data, 1);         % number of data 
Uinit= initmember(0.1,nclass,ndata);
[U, centroid, dist, W, F, obj] = gk_fkm(nclass,data,Uinit,phi,maxiter,toldif);
centroid
F{1}
F{2}

plot(data(:,1),data(:,2),'bo', centroid(1,:),centroid(2,:),'ro');

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
