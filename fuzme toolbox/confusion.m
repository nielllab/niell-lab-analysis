function ci = confusion(nclass,data,U)
% confusion index of Burrough and McDonnell (1998) 
% ci = confusion(nclass,data,U)
%   Input:
%   nclass=number of class
%   data= data matrix
%   U=membership matrix
us=U';
us=sort(us)';
ci=1-us(:,nclass)+us(:,nclass-1);
