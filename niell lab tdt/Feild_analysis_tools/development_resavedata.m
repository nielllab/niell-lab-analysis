% this just shows you how I resave it per session, have to change the directories.
cd F:\Jennifer_Development\Jen_analysis_development\Data_2_17_17
outputDir = 'F:/Jennifer_Development/Jen_analysis_development/Data_2_17_17/Bar';
mkdir(outputDir);
load('JLH_Development.mat')

for k = 1:size(exptdata,1)
  for j = 1:size(exptdata,2)
    data = exptdata{k,j};
    if isempty(data), continue,end
    isInSession = false;
    for iUnit = 1:length(unitdata)
      isInSession(iUnit) = unitdata{iUnit}.expnum==k & unitdata{iUnit}.GT==j;
    end
    unitdataSession = unitdata(isInSession);
    
    % save the data to the tank
    filename = fullfile(outputDir, exptdata{k,j}.analysis_file);
     [pathstr,name,ext] = fileparts(filename);
      save(name, 'unitdataSession', 'data');
  end
end