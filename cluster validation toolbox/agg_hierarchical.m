function [index, Cmerge] = agg_hierarchical(Distance, Method, k)
% agg_hierarchical Performs agglomerative clustering, the original file 
% DCAGG.m is renamed and altered by adding some codes, so that the 
% middle result stored in "Cmerge" may be returned. (Kaijun Wang July 2009)

% Cluster = DCAGG(Distance, Method, k) where Distance is
% square dissimilarity matrix, with Inf on leading diagonal,
% Method is one of 'single', 'complete' or 'centroid',
% and k is number of clusters. Cluster is a cell array showing
% which entities belong to which cluster.
%(C) David Corney (2000)   		D.Corney@cs.ucl.ac.uk

%If you want a nice hierarchical graph, try Jang's hierarchical
%clustering: http://neural.cs.nthu.edu.tw/jang/mlbook/demo/

[R,C] = size(Distance);
cluster = num2cell(1:C);		%Put every point is in its own cluster
k2 = k(2)-k(1)+1;
kn = k2;
Cmerge = cell(1,k2);
index = ones(R,k2);
k2 = R-k(2)+1;
m = R-k(1)+1;
Dmatrix = Distance;

for i = 2:m
   %Find closest clusters
   [MinRow, IdxDow] = min(Distance);
   [temp, MinJ] = min(MinRow);
   MinI = IdxDow(MinJ);
   
   if MinI > MinJ
      t=MinI;
      MinI=MinJ;
      MinJ=t;
   end
   
   %Merge cluster j into cluster i, then delete j
   C1 = cluster{MinI};
   C2 = cluster{MinJ};
   C = [C1 C2];
   cluster{MinI} = C;
   cluster(MinJ) = [];
   if i >= k2
     if strcmp(Method, 'single')
         md = min(min(Dmatrix(C1, C2)));
     elseif strcmp(Method, 'complete')
         md = max(max(Dmatrix(C1, C2)));
     else
         md = -1;
     end
      nc = length(C1);
      Cmerge{kn} = [md nc C];
      for j = 2:R-i+1
        index(cluster{j},kn) = j;
      end
      kn = kn-1;
   end
   %Calculate new Distance matrix
   switch Method
   case 'single'
      Distance(:, MinI) = min(Distance(:, MinI), Distance(:, MinJ)); 
      Distance(MinI, :) = min(Distance(MinI, :), Distance(MinJ, :)); 
   case 'complete'
      Distance(:, MinI) = max(Distance(:, MinI), Distance(:, MinJ)); 
      Distance(MinI, :) = max(Distance(MinI, :), Distance(MinJ, :)); 
   case 'centroid'
      Distance(:, MinI) = (Distance(:, MinI) + Distance(:, MinJ))/2; 
      Distance(MinI, :) = (Distance(MinI, :) + Distance(MinJ, :))/2;
   otherwise
      error('Unsupported Method in DCAGG!');
   end
   
   Distance(MinJ, :) = [];
   Distance(:, MinJ) = [];
   Distance(MinI, MinI) = inf;   
end