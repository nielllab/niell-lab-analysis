function [RM,RS,SPR,CD] = valid_hierarchical(data, labels, Cmerge, dtype, id)
% cluster validity indices for a hierarchical clustering algorithm

nc = size(labels,2);
RM = zeros(1,nc);
RS = zeros(1,nc);
CD = zeros(1,nc);
SPR = zeros(1,nc);

% computing R-Squared
for i = 1:nc
   if dtype == 1
     [St,Sw,Sb] = valid_sumsqures(data,labels(:,i),i+1);
   else
     [St,Sw,Sb] = valid_sumpearson(data,labels(:,i),i+1);
   end
   sst = trace(St);
   ssb = trace(Sb);
   RS(i) = ssb/sst;
end

if isempty(Cmerge)
   return
end

% computing RMSSTD, SPR, CD
if id
nc = length(Cmerge);
RM = zeros(1,nc);
for i = 1:nc
   C = Cmerge{i};
   CD(i) = C(1);
   m = C(2);
   C(1) = [];
   C(2) = [];
   nk = length(C);
   Clabel = [ones(m,1); 2*ones(nk-m,1)];
   if nk > 1
      R = data(C,:);
      Cm = mean(R);
      if dtype == 1
        R = R - repmat(Cm,nk,1);  %Cm(ones(nk,1),:)
        R = sum(R.^2);
        nk = length(R)-1;
        RM(i) = sqrt(sum(R)/nk);
        [St,Sw,Sb,Sintra,Sinter] = valid_sumsqures(data(C,:),Clabel,2);
        St = trace(St);
        Sw = trace(Sw);
      else
        R = similarity_pearsonC(R', Cm');
        nk = length(R)-1;
        RM(i) = sum(R)/nk;
        [St,Sw,Sb,Sintra,Sinter] = valid_sumpearson(data(C,:),Clabel,2);
      end
      SPR(i) = St-Sw;      % Semi-partial R-squared (SPR)
      if CD(i) == -1
         CD(i) = Sinter(3); % Distance between two clusters (CD)
      end
   end
end
RM = RM/max(RM);
CD = CD/max(CD);
SPR = SPR/max(SPR);
end
