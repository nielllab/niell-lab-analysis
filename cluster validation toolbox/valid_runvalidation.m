
S = [];
data = handles.data;
labels = handles.labels;
if id > 4 && (isempty(data) || size(data,1) ~= size(labels,1))
   Q = 'Validation state: Please load right data for computing validity indices !';
   set(handles.Outext2, 'String', Q);
   return;
end
ko = isempty(labels);
ke = sum(isnan(labels));
if ko || prod(ke)
   if ko
   Q = 'Validation state: No clustering sulution for validition !';
   else
   Q = 'Validation state: Class labels are "NaN" for validition !';
   end
   set(handles.Outext2, 'String', Q);
   return;
end

if id > 4
  Q = get(handles.Loadsolution, 'Checked');
  if strcmp(Q, 'on')
    distype = get(handles.popMeasure, 'Value');
  else
    distype = 2;
    if strcmp(handles.distype, 'Euclidean distance')
      distype = 1;
    end
  end
  
  if isempty(handles.Dmatrix)
    if distype == 1
       [Dmatrix, dmax] = similarity_euclid(data);
    else
       Dmatrix = 1-(1+similarity_pearson(data'))/2;
       dmax = 1;
    end
    handles.Dmatrix = Dmatrix;
    handles.dmax = dmax;
  else
    Dmatrix = handles.Dmatrix;
    dmax = handles.dmax;
  end
  if dmax ~= 1
    Dmatrix = Dmatrix/dmax;
  end
  
end
ke = handles.Cmerge;
truelabels = handles.truelabels;
k = handles.NC;
N = length(k);
ko = [];

Q = 'Validation state: Computation of validity indices is running';
Q = [Q ' (see Command window), please wait ...'];
set(handles.Outext2, 'String', Q);
pause(0.01);
disp(' ');
R = {[dname{id} ' index']};
kfind = [20 20 20 20 2 1 2 2 20 4 2 5 20 2 21 6 20 20];
%kfind = [20 20 20 20 2 1 2/3 2 3/20 1/4 2 5 20 2/3 2 6 20 20];
Q = get(handles.EstimateK, 'Checked');
if strcmp(Q, 'on') && N > 1
   Q = {'-b*', '-ro', '-m^', '-k+', '-gx', '-cs', '-yv', '-.mp', '-.bd'};
   dname{15} = 'Homogeneity';
   dname{16} = 'Separation';
else
   kfind = 21+kfind;
   if strcmp(dname{id}, 'All validity indices')
      id = 18;
      dname{10} = 'Homogeneity';
      dname{11} = 'Separation';
      dname{12} = 'All validity indices';
   end
   Q = {'b*', 'k*'};
end

if id < 5
   validty = valid_external(labels, truelabels);
   if isempty(validty)
      E = 'Validation state: Error class labels, not positive integers !';
      E = [E ' Please check pop-up menu "Position" or the data file.'];
      set(handles.Outext2, 'String', E);
      return;
   end
   if max(handles.truelabels) == 1
      E = 'Validation state: External index is not applicable !';
      set(handles.Outext2, 'String', E);
      return;
   end
   S = validty(id,:);
elseif id < 16
 Hs = strcmp(dname{id}, 'Homogeneity-Separation');
 Hs = Hs || strcmp(dname{id}, 'Homogeneity');
 if Hs && kfind(id) > 20
   validty = valid_internal(data, labels, k, Dmatrix, ke, distype, 21);
   S = validty([end-1 end],:); %[id-4 id-3]
   if distype ==1
      S = S*dmax;
   end
   R = {'Homogeneity', 'Separation'};
 elseif id == 13
   if ialg < 3 || ialg > 5
      set(handles.Outext2, 'String', ...
      'Validation state: RMSSTD group is only for Hierarchical results, please try another !');
      return;
   end
   S = valid_internal(data, labels, k, Dmatrix, ke, distype, id);
   R = {'RMSSTD', 'R-squared', 'Semi-partial R-squared', 'Centroid Distance'};
 else
   validty = valid_internal(data, labels, k, Dmatrix, ke, distype, id);
   S = validty(id-4,:);
 end
elseif id == 16
   Q = 'Validation state: none';
   set(handles.Outext2, 'String', Q);
   return;
   [S, ko] = valid_systemEvolution(k, labels, Dmatrix, distype);
   kfind(16) = -ko;
   disp(' ');
   R = {'Partition Energy', 'Merging Energy'};
    if id == 17
    [ke, ko] = valid_systemEvolution(k, labels, Dmatrix, distype);
    S = [S; ke];
    kfind(16) = -ko;
    end
elseif id >= 17
    S = valid_external(labels, truelabels);
    validty = valid_internal(data, labels, k, Dmatrix, ke, distype, id);
    S = [S; validty];
    R = {'All validity indices'};
end

axes(handles.axes1);
if single
   cla;
end

if id >= 17 && N > 1
   title('All validity indices in new figures');
   disp(dname(1:end-1));
   disp(S);
   ko = kfind;
   valid_plot_index(S, k, kfind, dname);
else

  if N > 1
     [ko, k] = valid_findk(S, kfind, id, k, N);
  end
  title('');
  if kfind(id) < 10
     title('Optimal k is indicated by a square symbol');
  end
  if id >= 17
     cla;
     ylim([0 1]); 
     title([R{1} ' at k = ' num2str(k)]);
     for i = 1:length(S)
       T = [dname{i} ' = ' num2str(S(i))];
       disp(T);
       text(k-0.9,1.02-0.08*i, T,'FontSize',11);
     end
  elseif single
     for j = 1:size(S,1)
       plot(k, S(j,:), Q{j});  hold on;
       disp([R{j} ' = ' num2str(S(j,:))]);
       if size(S,2) == 1
          text(k+0.1,S(j,:), num2str(S(j,:)),'FontSize',11);
       end
     end
     plot(k(ko),S(1,ko),'ks','MarkerSize',11);
  else
     Q = {'-b*','-ro','-m^','-k+','-gx','-cs','-yv','-.mp','-.bd','-.ro','-.m^','-.k+', ...
     '-.gx','-.cs','-.yv','-b*','-ro','-m^','-k+','-gx','-cs','-yv','-.mp','-.bd','-.ro','-.m^',};
     for j = 1:size(S,1)
       plot(k, S(j,:), Q{ialg});  hold on;
       disp([R{j} ' = ' num2str(S(j,:))]);
     end
     plot(k(ko),S(1,ko),'ks','MarkerSize',11);
     title(R); 
  end

  if id >= 17
     ylim([0 1]); 
     legend off;
  else
     if single
        legend(R,1);
     else
        T = [Salg ialg];
        Q = get(handles.popAlgorithm, 'String');
        high = find(T>10);
        if ~isempty(high)
          Q(T(high)) = handles.ialgFile;
        end
        legend(Q(T),1);
     end
  end
   ylim('auto');
   low = get(gca,'ylim');
   ymax = max([max(max(S)) low(2)]);
   ymin = min([min(min(S)) low(1)]);
   low(1) = ymin-0.05*abs(low(2)-low(1));
   low(2) = ymax+0.05*abs(low(2)-low(1));
   ylim(low);
   set(gca,'XTickLabelMode','auto');
   set(gca,'YTickLabelMode','auto');
   xlim2([k(1)-1 k(end)+1]);
   xlabel('number of clusters (k)','FontSize',11);
   ylabel('','FontSize',11);
end

if strcmp(R{1}, 'All validity indices')
   Q = dname;
else
   Q = R;
end

if N == 1&& single
  handles.indexname = Q;
  handles.indexvalue = S;
  handles.koptimal = k(ko);
elseif single
  handles.indexname = Q;
  handles.ialg = handles.algcluster;
  low = size(S,1);
  Svalue = NaN*ones(low,max(k));
  Svalue(:,k) = S;
  handles.ialgvalue = Svalue;
  if length(ko) == 1
     handles.koptimal = k(ko);
  else
     handles.koptimal = ko;
  end
else
  T = handles.ialgvalue;
  [nr, nc] = size(T);
  kmax = max([max(k) nc]);
  low = size(S,1);
  Svalue = NaN*ones(nr+low,kmax);
  Svalue(1:nr,1:nc) = T;
  Svalue(nr+1:nr+low,k) = S;
  handles.ialgvalue = Svalue;
  handles.indexname = Q;
  handles.koptimal = k(ko);
end
Q = handles.runstate;
if strcmp(R{1}, 'Homogeneity')
Q = ['Validation state: ' dname{id} ' values under ' handles.distype ' for ' Q];
else
Q = ['Validation state: ' R{1} ' values under ' handles.distype ' for ' Q];
end
handles.validstate = Q;
guidata(hObject,handles);

set(handles.Outext2, 'String', Q);