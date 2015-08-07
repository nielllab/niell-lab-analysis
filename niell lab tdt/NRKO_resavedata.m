% this just shows you how I resave it per session, have to change the directories.
cd ~/Downloads
outputDir = '~/Data/Hoy/';
mkdir(outputDir);
load('Hoy_adult_awake_NR5A_NR2B_KO_and_control_bars.mat')

for k = 1:size(exptdata,1)
  for j = 1:size(exptdata,2)
    data = exptdata{k,j};
    if isempty(data), continue,end
    isInSession = false;
    for iUnit = 1:length(unitdata)
      isInSession(iUnit) = unitdata{iUnit}.expnum==k & unitdata{iUnit}.cond==j;
    end
    unitdataSession = unitdata(isInSession);
    
    % save the data to the tank
    filename = fullfile(outputDir, exptdata{k,j}.tank);
    save(filename, 'unitdataSession', 'data');
  end
end