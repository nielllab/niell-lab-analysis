
set(hObject, 'Checked', 'on'); 
set(handles.EstimateK, 'Checked', 'off');
set(handles.Validityanalys, 'Checked', 'on');
set(handles.EstimateK, 'Checked', 'off');
dname = {'Rand', 'Adjusted Rand', 'Jaccard', 'FM', 'Silhouette', ...
                'Davies-Bouldin', 'Calinski-Harabasz', 'Dunn', 'R-Squared', ...
                 'Homogeneity-Separation', 'All validity indices'};
set(handles.popIndex, 'String', dname);
set(handles.popIndex, 'Value', 5);
Q = get(handles.Validone, 'Checked');
if strcmp(Q, 'on')
   valid_clear_validation;
   axes(handles.axes1);
   cla;
end
set(handles.Validone, 'Checked', 'off');