% this just shows you how I resave it per session, have to change the directories.
%cd D:\Jen_analysis\analysis_martin\data_1_25_16
cd F:\New_extraction
outputDir = 'F:/New_extraction/barAnalysis2';%
%outputDir = 'D:/Jen_analysis/analysis_martin/data_1_25_16/Hoy_bar_dataB';

mkdir(outputDir);%
load('JLH_NMDA_KO_bar_7_14_17.mat')

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
    save(filename, 'unitdataSession', 'data');
  end
end