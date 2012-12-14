% Plotting data by their class labels
[nr, nc] = size(data);
if nc > 2
  D = pca_of_data(data',2);
  D = D';
  xlabel('PC1','FontSize',11,'FontWeight','demi');
  ylabel('PC2','FontSize',11,'FontWeight','demi');
  title('Data are plotted in Principal Component space');
else
  D = data;
  xlabel('');
  ylabel('');
  title('Two dimensional data are plotted');
end

  if strcmp(id, 'bk')
     symb = ['k.';'ks';'k+';'ko';'k*';'k^';'k.';'k+';'ko';'k*';'ks';'k.'];
     %plot(D(:,1),D(:,2),'y.');
  elseif strcmp(id, 'co')
     symb = ['b.';'r.';'m.';'g.';'k.';'c.';'b+';'r*';'mo';'gs';'kd';'c^'];
     %plot(D(:,1),D(:,2),'y.');
  elseif strcmp(id, 'nb')
     symb = 1;
  end
  
k = max(labels); 
clas = cell(1,k);
for i = 1:k
    Q = find(labels==i);
    clas{i} = Q;
end

for i = 1:k
    id=clas{i};
  if symb == 1     
     text(D(id,1),D(id,2),['\fontsize{10}' int2str(i)]);
   else
     plot(D(id,1),D(id,2),symb(mod(i-1,12)+1,:));
  end
    hold on;
end

set(gca,'XTickMode','auto');
set(gca,'YTickMode','auto');
set(gca,'XTickLabelMode','auto');
set(gca,'YTickLabelMode','auto');
[fd,dt] = max(D(:,1));
[ft, dt] = min(D(:,1));
dt = 0.05*(fd-ft); 
xlim([ft-dt fd+dt]);
[fd,dt] = max(D(:,2));
[ft, dt] = min(D(:,2));
dt = 0.05*(fd-ft); 
ylim([ft-dt fd+dt]);
