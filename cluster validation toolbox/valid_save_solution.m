[file,filePath] = uiputfile('test.mat'); 
if file == 0
   return
end
NCrange = handles.NC;
labelsInNCrange = handles.labels;
id = handles.algcluster;
algorithmName = get(handles.popAlgorithm, 'String');
ResultTitle = ['Clustering & Validation results at k = ' num2str(NCrange)];

estimatedNC = [];
validityIndexName = handles.indexname;
validityIndexValue = handles.ialgvalue;
if ~isempty(validityIndexValue)
  id = handles.ialg;
  estimatedNC = handles.koptimal;
  if length(estimatedNC) > 1
     estimatedNC = [];
  end
end
  R = find(id>20);
  if ~isempty(R)
    Q = id(R);
    algorithmName(Q) = handles.ialgFile;
  end  
algorithmName = algorithmName(id);
save([filePath file], 'ResultTitle', 'algorithmName', 'labelsInNCrange', ...
   'validityIndexName', 'validityIndexValue', 'estimatedNC', 'NCrange');
