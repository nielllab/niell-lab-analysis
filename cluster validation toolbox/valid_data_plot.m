
labels = handles.labels;
ko = isempty(labels);
ke = sum(isnan(labels));
if  ko || prod(ke)
   R = '# Label plotting state: no clustering result, true labels (if any) are used !';
   %set(handles.Outext1, 'String', R);
   disp(R);
   labels = handles.truelabels;
  if isempty(labels)
     return;
  end
   title('Data are plotted in class labels');
end

data = handles.data;
ko = handles.koptimal;
if length(ko) > 1
ko = [];
end
ns = size(labels,2);
if ns > 1
  k = handles.NC;
  ko = num2str(ko);
  R = [num2str(k(1)) '-' num2str(k(end))];
  R = ['Select one of k = ' R];
  dlg_title = 'Selection';
  ko = inputdlg(R, dlg_title, 1, {ko});
  if isempty(ko)
     return;
  end
  ko = str2double(ko{1});
  kp = find(k == ko);
  title('Solutions in class labels & Data in colors if shown');
end

if isempty(ko) ||  ns == 1
   kp = 1;
end
labels = labels(:,kp);
id ='nb';
plot_data_bylabels;
box on;