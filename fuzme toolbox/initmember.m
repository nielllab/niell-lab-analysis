function U= initmember(scatter,nclass,ndata)

if(scatter<=0), scatter=0.001,end;
% initialise membership
U=1/nclass*ones(ndata,nclass);
U= U + scatter* rand(ndata,nclass);
cs=sum(U')';
U = U./cs(:,ones(nclass,1));
