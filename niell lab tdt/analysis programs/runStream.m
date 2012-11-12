[fname pname] = uiputfile('','data folder');

snipName =[fname 'snip.mat']

%choose_pca = input('manually choose pca (0=no / 1 = yes ?')
analyzeStream(snipName,pname);

cluster_tetrode_stream(snipName,pname,1);
