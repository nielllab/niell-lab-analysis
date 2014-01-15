function valid_plot_index(validty, ks, B, dname)
% preparing for plotting indices

m = size(validty,1);
N = length(ks);
Q = {'b*','r-'};

for i=1:m
  R = i;
  if 0 %i == m
    R = R-1;
    Q = {'b^','m-'};
    dname{i} = [dname{i} '-' dname{i-1}];
  end
  R = mod(R-1,4)+1;
  if R == 1
     figure('color','white');
  end
  subplot(2,2,R); 
  
nk = ks;
S = validty(i,:);
[high, ko] = max(S);
[low, ke] = min(S);
[ko, nk] = valid_findk(S, B(i), 1, nk, N);

  margin = (high-low)*0.1;
  if i == m-1
     oldhigh = high;
     oldlow = low;
  elseif 0 %i == m
     high = max([high oldhigh]);
     low = min([low oldlow]);
     margin = (high-low)*0.1;
  end
  if margin == 0
     margin = 0.1;
     S(1:end) = NaN;
  end
  
  plot(nk,S,Q{2});  hold on;
  plot(nk,S,Q{1}); hold on;
  plot(nk(ko),S(ko),'ks','MarkerSize',11);
  xlim2([nk(1) nk(end)]);
  xlim([nk(1)-0.3 nk(end)+0.3]);
  ylim([low-margin high+margin]);
  
  xlabel('number of clusters(k)','FontSize',11,'FontWeight','demi');
  title([dname{i} ' index']);
end