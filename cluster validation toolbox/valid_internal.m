function validty = valid_internal(data, labels, nk, Dmatrix, Cm, dtype, id)
% calculating intenal validity indices

N = length(nk);
nrow = size(data,1);
if N == 1
  m = 1;
else
  m = 1:N;
end
ids = (id > 5 && id < 9 || id == 11 || id == 12 || id > 16 && id < 21);

  Sil = zeros(1,N);
  DB = zeros(1,N);
  CH = zeros(1,N);
  Dunn = zeros(1,N);
  RS = zeros(1,N);
  Cindex = zeros(1,N);
  KL = zeros(1,N);
  Ha = zeros(1,N);
  RM = zeros(1,N);
  Hom = zeros(1,N);
  Sep = zeros(1,N);
  wtertra = zeros(1,N);

if id == 5 || id > 16 && id < 21
  if dtype == 1
    R = 'euclidean';
  else
    R = 'correlation';
  end
  for k = m
    S = silhouette(data, labels(:,k), R);
    Sil(k) = mean(S);
  end
end

if ids
  for k = m
    [DB(k), CH(k), Dunn(k), KL(k), Ha(k), ST] = ...
        valid_internal_deviation(data,labels(:,k), dtype);
   disp(['Validation state: running at k = ' num2str(nk(k))]);
  end
end
if N > 1 && ids
  S = trace(ST);
  KL = [S KL];
  Ha = [S Ha];
  R = abs(KL(1:N)-KL(2:N+1));
  S = [R(2: end) R(end)];
  KL = R./S;
  KL(N) = KL(N-1);
  R = Ha(1:N)./Ha(2:N+1);
  Ha = (R-1).*(nrow-[nk(1)-1 nk(1:N-1)]-1); 
end

ids = (id == 9 || id == 13 || id > 16);
if ids
   R = (id > 9);
   [RM,RS,SPR,CD] = valid_hierarchical(data, labels, Cm, dtype, R);
end

ids = (id == 10 || id > 13);
if ids
  for k = m
     S = ind2cluster(labels(:,k));
     [Hom(k), Sep(k), Cindex(k), wtertra(k)] = ... %, Dunn(k), DB(k)
         valid_internal_intra(Dmatrix, S, dtype, id==10);
   disp(['Validation state: running at k = ' num2str(nk(k))]);
  end
end

if id == 13
  validty = [RM; RS; SPR; CD];
elseif id < 18
  %R = (1-sqrt(nk).*nk/nrow);
  %Hom = R.*(Hom-Sep);
  %if dtype == 1
  %   Hom = -Hom;
  %end
  validty = [Sil; DB; CH; Dunn; RS; Cindex; KL; Ha; RM; wtertra; Hom; Sep];
else
  validty = [Sil; DB; CH; Dunn; RS; Hom; Sep];
end
