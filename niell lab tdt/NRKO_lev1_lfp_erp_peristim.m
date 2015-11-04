function [] = NRKO_lev1_lfp_erp_peristim

% NRKO_LEV1_ERP produces spike triggered spectra of the LFP

% take the globals from the info script that has been ran at the start
global info
global outputDir
input   = 'lev0_lfp_peristim';
output  = 'lev0_lfp_erp_peristim';

% loop over the various files
nDirs = length(info);
for iDir = 1:nDirs
  
  load(fullfile(outputDir, input, info(iDir).dataname, input));
    
  cfg = [];
  erp = ft_timelockanalysis(cfg, dataTrials);
  
  % save the STS
  filename = fullfile(outputDir, output, info(iDir).dataname, output);
  mkdir(fullfile(outputDir, output, info(iDir).dataname));
  save(filename, 'erp')
   
end

