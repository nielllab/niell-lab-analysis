function [] = NRKO_lev1_sts_convol_peristim

% NRKO_LEV1_ERP produces spike triggered spectra of the LFP

% take the globals from the info script that has been ran at the start
global info
global outputDir
input   = 'lev0_lfp_peristim';
output  = 'lev0_lfp_tfr_peristim';

% loop over the various files
nDirs = length(info);
for iDir = 1:nDirs
  
  load(fullfile(outputDir, input, info(iDir).dataname, input));
     
  cfg           = [];
  cfg.foi       = 2:2:160;
  cfg.t_ftimwin = 5./cfg.foi;
  cfg.toi       = -0.1:0.05:3;     
  cfg.method    = 'mtmconvol';
  cfg.taper     = 'hanning';
  tfr = ft_freqanalysis(cfg, dataTrials);
   
  % save the STS
  filename = fullfile(outputDir, output, info(iDir).dataname, output);
  mkdir(fullfile(outputDir, output, info(iDir).dataname));
  save(filename, 'sts')   
end

