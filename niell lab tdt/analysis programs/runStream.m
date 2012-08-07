[fname pname] = uiputfile('','data folder');

snipName =[fname 'snip.mat']

analyzeStream(snipName,pname);
cluster_tetrode_stream(snipName,pname);
