function [dout,wide,dmin]=data_normalization(data,choise,oldwide,oldm)
% each column is normalized by its center, variance or max/min values
% All rights Reserved, March 2006, Feb. 2007. Kaijun WANG

[nrow,ncol] = size(data);
dout = data;
dmin = 0;
wide = ones(1,ncol);

 if oldwide(1) > 0
  wide = oldwide;
  dmin = oldm;
 elseif choise == 1          % normalization: zero mean & unit variance
  dmin = mean(data);
  wide = std(data,0,1);     % standard variance
  R = find(wide==0);
  wide(R) = dmin(R);
elseif choise == 2 || choise == 3 
  dmax = max(data);
  dmin = min(data);
  wide = dmax-dmin;        % standardization: [0,1] for each column
 end
 
if choise == 1 || choise == 2 || choise == 3
   data = data - repmat(dmin,nrow,1);
   dout = data./repmat(wide,nrow,1);
end
if choise == 3                     % standardization: [-1,1] for each column
   dout = dout+dout-1;
end