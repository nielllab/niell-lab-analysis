% loading a data set

Q = get(handles.Loadsolution, 'Checked');
if ~strcmp(Q, 'on')
  valid_clear_clustering;
end
  valid_clear_validation;
  valid_clear_plotting;

id = get(handles.ClusterColumn, 'Checked');
if strcmp(id, 'on')
   data = data';
end
[nrow, dim] = size(data);

id = get(handles.popPosition, 'Value');
if id == 1
   handles.truelabels = data(:,1);
   handles.data = data(:,2:dim);
elseif id == 2
   handles.truelabels = data(:,dim);
   handles.data = data(:,1:dim-1);
else
   handles.truelabels = ones(nrow,1);
   handles.data = data;
end
C = unique(handles.truelabels);
k = length(C);
handles.koptimal = k;
handles.dname = file; 
handles.datastate = R;
guidata(hObject,handles);

Q = handles.runstate;
if ~isempty(Q)
   R = [R ' Notice: ' Q ' exist.'];
else
  id = get(handles.ClusterColumn, 'Checked');
  if strcmp(id, 'on')
   R = [R ' Notice: columns will be clustered.'];
  else
   R = [R ' Notice: rows will be clustered (default).'];
  end
end
set(handles.Outext1, 'String', R);
set(handles.Outext2, 'String', 'Validation state: none');
