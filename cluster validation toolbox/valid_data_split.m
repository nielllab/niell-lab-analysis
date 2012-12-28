icol = get(handles.ClusterColumn, 'Checked');
if strcmp(icol, 'on')
   data = data';
end
[nrow, dim] = size(data);

icol = get(handles.popPosition, 'Value');
if icol == 1
   handles.truelabels=data(:,1);
   handles.data = data(:,2:dim);
elseif icol == 2
   handles.truelabels=data(:,dim);
   handles.data = data(:,1:dim-1);
else
   handles.truelabels=ones(nrow,1);
   handles.data = data;
end
handles.koptimal = max(handles.truelabels);
guidata(hObject,handles);