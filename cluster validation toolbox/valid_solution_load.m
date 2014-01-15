[file,filePath] = uigetfile('*.txt'); 
if isequal(file, 0)  
   return; 
end
dname = [filePath file];
try
  labels = load(dname);
catch
  set(handles.Outext1, 'String', 'Running state: incorrect solution file !');
  return
end

valid_clear_clustering;

[nrow, dim] = size(labels);
Q = labels(1,:);
labels = labels(2:nrow,:);
if Q(1) == 0
   handles.NC = Q(2:dim);
   handles.truelabels = labels(:,1);
   handles.labels = labels(:,2:dim);
else
   handles.NC = Q;
   handles.truelabels=ones(nrow-1,1);
   handles.labels = labels;
end

handles.algcluster = file;                    % 'solution file'
C = unique(handles.truelabels);
handles.koptimal = length(C);
handles.runstate = ['solution file - ' file];
guidata(hObject,handles);

R = ['Running state: Solution file ' file ' is loaded. Load its'];
R = [R ' corresponding data set (appears in popup menu "Data'];
R = [R ' set" if loaded) before computing internal indices.'];
set(handles.Outext1, 'String', R);
set(handles.Outext2, 'String', 'Validation state: none');
set(hObject, 'Checked', 'on');
set(handles.Loadata, 'Checked', 'off');
set(handles.Loademo, 'Checked', 'off');

%valid_clear_validation;
%valid_clear_plotting;