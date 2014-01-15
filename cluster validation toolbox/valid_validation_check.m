if isempty(handles.runstate)
   Q = ['Validation state: Please run clustering or load a solution file first !'];
   set(handles.Outext2, 'String', Q);
   return;
end
ialg = handles.algcluster;
id = get(handles.popIndex, 'Value');
dname = get(handles.popIndex, 'String');
single = get(handles.Validmutiple, 'Checked');
if strcmp(dname{id}, 'All validity indices') && strcmp(single, 'on')
   Q = ['Validation state: Please click "Single-Algorithm Validation" '];
   Q = [Q 'in Validation menu to to use it !'];
   set(handles.Outext2, 'String', Q);
   return;
end

if strcmp(single, 'on')
   Q = handles.indexname;
 if ~isempty(Q)
   if strcmp(Q{1}, 'Homogeneity') || strcmp(Q{1}, 'Separation')
      Q = {'Homogeneity-Separation'};
   end
   R = strcmp(dname{id}, Q{1}) || strcmp([dname{id} ' index'], Q{1});
   if ~R || strcmp(dname{id}, 'All validity indices')
      R = ['Validation state: New selection is not ' Q{1} ', please turn to it,'];
      R = [R ' or press "Clear Validation" to start a new work']; %& "Clear Plotting"
      set(handles.Outext2, 'String', R);
      return
   end
 end
   single = 0;
   Salg = handles.ialg;
   if ischar(ialg) || isempty(find(Salg==ialg))
      if ischar(ialg)
        if isempty(Salg) || max(Salg) < 21
          Q{1} = ialg;
          handles.ialgFile = Q;
          ialg = 21;
        else
          Q = handles.ialgFile;
          R = max(Salg)-20+1;
          Q{R} = ialg;
          handles.ialgFile = Q;
          ialg = R+20;
        end
      end
      valid_runvalidation;
      
      if isempty(S)
        return;
      end
      Salg = [Salg ialg];
      handles.ialg = Salg;
      guidata(hObject,handles);
   else
      Q = get(handles.popAlgorithm, 'String');
      R = ['Validation state: Validation results for ' Q{ialg} ' exist already, please'];
      R = [R ' select a new algorithm, or press "Clear Validation" to start a new work'];
      set(handles.Outext2, 'String', R);
   end

else
   single = 1;
   valid_runvalidation;
end
