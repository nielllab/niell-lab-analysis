

[fname, pname] = uigetfile('*.mat','cluster data');
clustfile = fullfile(pname,fname);

[fname, pname] = uigetfile('*.mat','analysis data');
afile = fullfile(pname,fname);

[fname, pname] = uigetfile('c:\data\movie files','white noise movie file');
wnfile = fullfile(pname,fname);

[fname, pname] = uigetfile('c:\data\movie files','flash spots movie file');
flfile = fullfile(pname,fname);

[fname, pname] = uigetfile('c:\data\movie files','mv spots movie file');
mvfile = fullfile(pname,fname);


[fname, pname] = uiputfile('*.mat','pdf file');
pdfFile = fullfile(pname,fname);

load(clustfile,'Block_Name');

for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
barblock = input('which block for  bars (0=none) ? ');

for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end

driftblock = input('which block for drift ? ');

for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
flblock = input('which block for flash spots ? ');

for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
mvblock = input('which block for mv spots (0=none)? ');

for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
contrablock = input('which block for contra movie ? ');


for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
ipsiblock = input('which block for  ipsi movie (0=none) ? ');

for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
moveblock = input('which block for movement analysis (0=none) ? ');




if barblock~=0
    bars_analysis(clustfile,afile,pdfFile,Block_Name,barblock)   %%%% moving spots
    close all
end

if driftblock~=0
    drift_analysis(clustfile,afile,pdfFile,Block_Name,driftblock);
close all
fclose all
end

noise_analysis(clustfile,afile,pdfFile,wnfile,Block_Name,contrablock,1, 1);  %%% white noise contra
close all
fclose all

if ipsiblock~=0
    noise_analysis(clustfile,afile,pdfFile,wnfile,Block_Name,ipsiblock,1, 2)  %%% white noise ipsi
    close all
    fclose all
end

noise_analysis(clustfile,afile,pdfFile,flfile,Block_Name,flblock,2)  %%% flashing spots
close all
fclose all

if mvblock~=0
    noise_analysis(clustfile,afile,pdfFile,mvfile,Block_Name,mvblock,3)   %%%% moving spots
close all
fclose all
end

if moveblock~=0
    pptgAnalysis(clustfile,afile,pdfFile,Block_Name,moveblock,0,60)
 
    fclose all
end
