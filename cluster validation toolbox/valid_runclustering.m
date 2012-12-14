ialg = get(handles.popAlgorithm, 'Value');
if ialg > 20
  disp(' # # Clustering process is not applicable to this selection !'); 
  return
end

data = handles.data;
[nrow, dim] = size(data);

k = get(handles.edit1, 'String');
k = str2double(k);
k2 = str2double(k2);
if  k < 2 || k > nrow || int8(k) ~= k || k2 < k || k2 > nrow || int8(k2) ~= k2
   set(handles.Outext1, 'String', ' Error ! Input Number of clusters is wrong or too big.');
   set(handles.Outext2, 'String', 'Running state: Clustering stops on acount of errors');
   return;
end

if k2 > k
   Q = ['at k = ' num2str(k) ' - ' num2str(k2)];
   ko = max(handles.truelabels);
   if ko >= k && ko <= k2
      kp = ko-k+1;
   else
      kp = 1;
   end
else
   Q = ['at k = ' num2str(k)];
   ko = 1;
   kp = 1;
end

%icol = get(handles.LoadSmatrix, 'Checked');
if 1 %strcmp(icol, 'off')
   R = get(handles.Normaliz, 'Checked');
   icol = get(handles.Standz, 'Checked');
   if strcmp(R, 'on')
     data = data_normalization(data',1,0); % each column to 0 mean & variance 1
     data = data';
   elseif strcmp(icol, 'on')
     data = data_normalization(data',2,0);  % columns are standardized within [0,1]
     data = data';
   end
   Dist = [];
else
   Dist = handles.Smatrix;
end

icol = get(handles.popAlgorithm, 'String');
distype = get(handles.popMeasure, 'Value');
R = get(handles.popMeasure, 'String');
handles.distype = R{distype};
labels = ones(nrow,k2-k+1);
NC = [];
Cmerge = [];
dmax = 1;
if ialg > 2 && ialg < 6 || ialg > 7
   [labels,NC,Dist,dmax,Cmerge] = valid_clusteringAlgs(data,ialg,[k k2],distype,Dist,dmax);
else
  ki = 10;
  if ialg == 2
    R = {'Number of times to repeat the clustering for returning a best solution'};
    dlg_title = 'replicates';
    ki = inputdlg(R, dlg_title, 1, {'10'});
    if isempty(ki)
       return;
    else
       ki = str2double(ki{1});
    end
  end
for i = 1:k2-k+1
   [labels(:,i),NC,Dist,dmax] = valid_clusteringAlgs(data,ialg,[i+k-1 k2],distype,Dist,dmax,ki);
   R = ['Running state: ' icol{ialg}  ' Clustering is running at iteration = ' num2str(i)];
   set(handles.Outext1, 'String', R);
   pause(0.01);
end
end

dim = [];
if strcmp(get(handles.Standz, 'Checked'), 'on')
   dim = 'standardized';
elseif strcmp(get(handles.Normaliz, 'Checked'), 'on')
   dim = 'normalized';
end
dname = handles.dname;
R = ['Clustering state: ' icol{ialg} ' clustering on ' dim ' ' dname ' ' Q  ' runs over.'];
handles.runstate = [icol{ialg} ' results'];
if ~isempty(NC)
handles.NC = NC;
else
handles.NC = k:k2;
end
handles.labels = labels;
handles.Cmerge = Cmerge;
handles.Dmatrix = Dist;
handles.dmax = dmax;
handles.algcluster = ialg;
guidata(hObject,handles);

set(handles.Loadsolution, 'Checked', 'off');
set(handles.Outext1, 'String', R);
set(handles.Outext2, 'String', 'Validating state: none');
fprintf(['\n' R '\n']);

if 0
if k == k2
   R = [R ' Output class labels are in Command Window.'];
else
   R = [R ' Output class labels may be saved with "Save Results" in File menu'];
   R = [R ' or printed with "Print Class Lables" in View menu.'];
end
distype = ['===> ' icol{ialg} ' algorithm on data set ' dname '  runs over'];
disp(' ');
disp(distype);
if ko > k && ko <= k2
  R = num2str(ko);
else
  R = num2str(k);
end
fprintf(['    and gives the clustering result of class labels at k = ' R ' :\n']);
disp(num2str(labels(:,kp)'));
end