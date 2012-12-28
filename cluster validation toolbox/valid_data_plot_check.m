id_return = 0;
if isempty(handles.data)
   R = 'Running state: Please load data first for using this function !';
   set(handles.Outext1, 'String', R);
   id_return = 1;
   return;
end
id = get(handles.Plotnewfigure, 'Checked');
if strcmp(id, 'on')
   figure('color','white');
else
   axes(handles.axes1);
end