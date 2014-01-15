function [pcfirst,eigvect,eigval,pcorder]=pca_of_data(data,np)
% extracting PCA components from data (along dimensions)

[m,n]=size(data);
datcov=cov(data');
[V,S] = eig(datcov);
eigvalue= diag(S);
[Vs,eigorder] = sort(-abs(eigvalue));
eigvalue  = eigvalue(eigorder);
Vs = V(:,eigorder);

for i = 1:np
    V(:,i) = Vs(:,i)/norm(Vs(:,i));
end

eigvect = V(:,1:np);           % first eigenvectors
pcfirst=eigvect'*data;
pcorder=eigorder(1:np);
eigval = eigvalue(1:np);
%eigval = abs(eigvalue(1:np)); % /sum(abs(eigvalue));
