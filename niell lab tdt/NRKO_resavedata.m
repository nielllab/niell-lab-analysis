% this just shows you how I resave it per session, have to change the directories.
cd D:\Jen_analysis\NR5A_Pinping\Jen_NR5A_analysis_files\analysis_files
outputDir = 'D:/Jen_analysis/NR5A_Pinping/Jen_NR5A_analysis_files/analysis_files/Martin_analysis/Data/Hoy/';
mkdir(outputDir);
load('JLH_NMDA_KO_drift_test1.mat')

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