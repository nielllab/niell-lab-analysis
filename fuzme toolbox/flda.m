function [F,WC,BC,wilks,V,e,vexp,z,zc]=flda(data,nclass,U,centroid,phi,scale)
% Fuzzy linear discriminant analysis
%   Input:
% data      : data matix        data(ndata,ndim)
% nclass    : no. classes
% U         : membership matrix U(ndata,nclass)
% centroid  : centroid of data  centroid(nclass,ndim)
% phi       : fuzzy exponent
% scale     : type of scaling applied to data 
%               0= no scaling, 
%               1= center data, mean = 0
%               2= center & standardise, mean = 0, sd = 1
%   Output:
% F     : fuzzy covariance matrix
% WC    : SSP within class
% BC    : SSP between classes
% wilks : fuzzy wilks criterion
% V     : eigenvector of canonical variates
% e     : eigenvalues
% vexp  : percent variance explained
% z     : projection of data into canonical axis
% zc    : projection of centroid into canonical axis
%
ndata=size(data,1);
ndim=size(data,2);

xmean=mean(data);
xstd=std(data);

% WC: SSP within class
% BC: SSP between classes
% F: fuzzy covariance matrix
BC=zeros(ndim,ndim);
WC=zeros(ndim,ndim);
for k=1:nclass,
    ufi=U(:,k).^phi;
    c1=data-repmat(centroid(k,:),ndata,1);
    c2=repmat(ufi,1,ndim).*c1;
    c3=c2'*c1;
    WC=WC+c3;
    F{k}=c3./sum(ufi);
        
    c4=centroid(k,:)-xmean;
    c5=sum(ufi).*c4;
    c6=c5'*c4;
    BC=BC+c6;        
end
clear c1 c2 c3 c4 c5 c6 ufi;

% Total SSP
TC=WC+BC
% Wilks Lambda
wilks=det(WC)/det(TC);
    
% canonical variates
WB=inv(WC)*BC;
[V,D] =eig(WB); % D = eigenvalue; V = eigenvector
e=diag(D);
tvar = sum(e);          % total variance
vexp = 100*e./tvar;     % varaince explained
   
% center data
if(scale==0),
    xmean=zeros(1,ndim);
    xstd=ones(1,ndim);
elseif(scale==1),
    xstd=ones(1,ndim);
end

cdata = (data - repmat(xmean,ndata,1))./ repmat(xstd,ndata,1);
centr = (centroid - repmat(xmean,nclass,1))./repmat(xstd,nclass,1);

% project data & centroid into canonical axis
z=cdata*V;
zc=centr*V;
  
% plot Canonical variate
figure;
hold on;
plot(z(:,1),z(:,2),'o');
plot(zc(:,1),zc(:,2),'r+');
hold on;

if(scale==2),
    % plot eigenvector
    for j=1:ndim,
        xp=[0; V(1,j)];
        yp=[0; V(2,j)];
        plot(xp,yp,'m-');
    end;
end;
hold off;
