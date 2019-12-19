%PathList_New
n= 0;
[outfile outpath] = uiputfile('*.txt','text file for output');

tic
for iDir = 1:length(pname);
   iDir
   files = dir([ pname{iDir} '*eye*.avi']);   
    for f = 1:length(files)  
        sprintf('dir %d/%d   file %d/%d',iDir,length(pname),f, length(files))
        if isempty(strfind(files(f).name,'DeInter'))
            n = n+1;
            name = fullfile(files(f).folder,files(f).name)
            movieFiles{n} = deInterlaceVids(name);
        end
    end
end

toc

filePh = fopen(fullfile(outpath,outfile),'w');
fprintf(filePh,'[r''%s'',\n',movieFiles{1});
for f = 2:n-1
    fprintf(filePh,'r''%s'',\n',movieFiles{f});
end
fprintf(filePh,'r''%s'']\n',movieFiles{n});
fclose(filePh);


    