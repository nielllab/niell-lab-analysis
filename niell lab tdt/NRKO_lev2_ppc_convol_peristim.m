function [] = NRKO_lev2_ppc_convol_peristim

% take the globals from the info script that has been ran at the start
global info
global outputDir
input   = 'lev1_sts_convol_peristim';
output  = 'lev2_ppc_convol_peristim';

% loop over the various files
nDirs = length(info);

for iDir =1:nDirs
  load(fullfile(outputDir, input, info(iDir).dataname, input));
  
  % different ppc versions, angle, rayleigh test
  methods  = {'ppc0', 'ppc1', 'ppc2', 'ang', 'ral'};%usually used ppc1 most conservative/sesnsitive 
  nMethods = length(methods); 
  clear stat stat_tfr   
  for iMethod = 1:nMethods
    nUnits = length(sts.label);
    unitCnt = 0;
    for iUnit = 1:nUnits      
      unitCnt = unitCnt + 1;
      if size(sts.fourierspctrm{iUnit},1)<10, continue,end % reject less than 10 spikes in total
      cfg = [];
      cfg.spikechannel = iUnit;
      cfg.method = methods{iMethod};      
      cfg.avgoverchan = 'unweighted'; % average the phases prior to computing phase locking
      % we need to add something for rejecting same spike-field combinations here.
      % but for that I need more info % compute per CSD 
      stat(iMethod).static(unitCnt) = ft_spiketriggeredspectrum_stat(cfg, sts);
     end
  end  
  
  % save the STS
  filename = fullfile(outputDir, output, info(iDir).dataname, output);
  mkdir(fullfile(outputDir, output, info(iDir).dataname));
  if exist ('stat', 'var')
  save(filename, 'stat')
  end
end