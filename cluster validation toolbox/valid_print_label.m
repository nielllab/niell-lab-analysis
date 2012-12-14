labels = handles.labels;
if isempty(labels)
  labels = handles.truelabels;
  if isempty(labels)
     return;
  end
  k = max(handles.truelabels);
else
  k = handles.NC;
end
R = length(k);
if R > 1
  R = ['Print result at which k (' num2str(k(1)) '-' num2str(k(R)) ') ?'];
  R = {R, 'End number of clusters'};
  dlg_title = 'Select k';
  ko = inputdlg(R, dlg_title, 1, {'Start number of clusters', ''});
  if isempty(ko)
     return;
  end
  ko = [str2double(ko{1}) str2double(ko{2})];
  if ko(1) == 0
    labels = handles.truelabels;
    ko = max(handles.truelabels);
    kp = 1;
    R = ['## The true class labels at k = ' num2str(ko) ' :\n'];
  else
  if ko(2) < ko(1) || isnan(ko(2))
     ko(2) = ko(1);
  end
  ko = ko(1):ko(2);
  kp = find(k == ko(1)):find(k == ko(end));
  R = num2str(ko);
  R = ['## The clustering result (class labels) at k = ' R ' :\n'];
  end
else
  ko = k;
  kp = 1;
  R = num2str(ko);
  R = ['## The class labels at True k = ' R ' :\n'];
end
fprintf(R);
disp(num2str(labels(:,kp)'));
