labels = handles.labels;
truelabel = handles.truelabels;
if isempty(labels) || isempty(truelabel)
   R = 'Validation state: Please run clustering or load a solution file first !';
  set(handles.Outext2, 'String', R);
  return;
end
C = unique(handles.truelabels);
ko = length(C);
NC = handles.NC;
k = find(NC == ko);
if ~isempty(k)
  fprintf('# clustering solution (class labels) compared with true labels (if any):');
  %E = valid_errorate(labels(:,k), truelabel);
else
  k = 1;
end
  R = {['Which k= ' num2str(NC)]};
  dlg_title = 'Selecting a NC';
  k = inputdlg(R, dlg_title, 1, {num2str(NC(k))});
  if ~isempty(k)
    k = str2double(k{1});
  end
  if isempty(k) || k < NC(1) || k > NC(end)
    fprintf('# Error Rate is not applicable in this case !');
    return
  end
  k = find(NC == k);
  E = valid_errorate(labels(:,k), truelabel);
Q = 'Validation state: Error rate of the clustering solution (or class labels at k= ';
Q = [Q num2str(NC(k)) ') compared with true labels (at k= ' num2str(ko) ') : '];
Q = [Q num2str(E) '%'];
set(handles.Outext2, 'String', Q);
