
k = handles.NC;
Q = {'-b*','-ro','-m^','-k+','-gx','-cs','-yv','-.mp', '-.bd', ...
   '-.ro','-.m^','-.k+','-.gx','-.cs','-.yv','-b*','-ro','-m^','-k+','-gx'};
Sname = handles.indexname;
if strcmp(Sname{1}, 'Hartigan index') || strcmp(Sname{1}, 'Hartigan')
   k = k-1;
end
S = handles.indexvalue;
  
if length(k) == 1 && ~isempty(S)
  m = length(Sname);
  if m > 3
     ylim([0 1]); 
     xlim2([k(1)-1 k(end)+1]);
     for i = 1:length(S)
       R = [Sname{i} ' = ' num2str(S(i))];
       text(k-0.9,1.02-0.08*i, R,'FontSize',11);
     end
     title(Sname{i+1});
     box on;
  else
    for j = 1:m
       plot(k, S(j,:), Q{j});  hold on;
       text(k+0.1,S(j,:), num2str(S(j,:)),'FontSize',11);
    end
     legend(Sname,1);
  end
  
else
  Salg = handles.ialg;
  Svalue = handles.ialgvalue;
  if isempty(Svalue)
     return
  end
  
  single = get(handles.Validone, 'Checked');
  if strcmp(single, 'on')
     dname = handles.indexname;
     Salg = 1:length(dname);
  else
     dname = get(handles.popAlgorithm, 'String');
  end
  
  if length(k) > 1
    m = size(Svalue,2);
    if m > length(k)
      % k = 1:m;
    end
    if ~strcmp(dname(end), 'All validity indices')
      ymin = num2str(k(1));
      ymax = num2str(k(end));
      R = ['Redrawing at which k (' ymin  '-'  ymax ') ?'];
      R = {R, 'End number of clusters'};
      dlg_title = 'Select k range';
      ko = inputdlg(R, dlg_title, 1, {ymin, ymax});
      if ~isempty(ko) || ~strcmp(dname(end), 'All validity indices')
        k = str2double(ko{1}):str2double(ko{2});
      end
    end
  end
  
  ko = handles.koptimal;
  if strcmp(dname(end), 'All validity indices')
     title('All validity indices in new figures');
     valid_plot_index(Svalue(:,k), k, ko, dname);
     return;
  end
  
  m = length(Salg);
  n = round(size(Svalue,1)/m);
  for i = 1:n
    for j = 1:m
      md = i+n*(j-1);
      plot(k, Svalue(md,k), Q{j});  hold on;
      %disp([dname{Salg(j)} ' = ' num2str(Svalue(md,:))]);
    end
  end
  plot(ko,Svalue(1,ko),'ks','MarkerSize',11);
  
  S = strcmp(Sname{1}, 'Homogeneity');
  if S
     Q = get(handles.Validmutiple, 'Checked');
     S = strcmp(Q, 'on');
     Q = strcmp(handles.distype, 'Euclidean distance');
     if S && Q
        Sname = {'Separation', 'Homogeneity'};
     end
  end
  
  if strcmp(Sname{1}, 'RMSSTD')
    title('');
  else
    title(Sname);
  end
  R = find(Salg>20);
  if ~isempty(R)
    Q = Salg(R);
    dname(Q) = handles.ialgFile;
  end
  legend(dname{Salg},1);
  xlim2([k(1)-1 k(end)+1]);
  ylim('auto');
  low = get(gca,'ylim');
  ymax = max(max(Svalue(:,k)));
  ymin = min(min(Svalue(:,k)));
  low(1) = ymin-0.08*abs(low(2)-low(1));
  low(2) = ymax+0.08*abs(low(2)-low(1));
  ylim(low);
  set(gca,'XTickLabelMode','auto');
  set(gca,'YTickLabelMode','auto');
  xlabel('number of clusters(k)','FontSize',11);
  ylabel('','FontSize',11);
end
