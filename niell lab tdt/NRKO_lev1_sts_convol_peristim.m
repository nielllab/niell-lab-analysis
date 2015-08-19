function [] = NRKO_lev1_sts_convol_peristim

% NRKO_LEV1_STS_PERISTIM produces spike triggered spectra of the LFP

% take the globals from the info script that has been ran at the start
global info
global outputDir
input1  = 'lev0_lfp_peristim';
input2  = 'lev0_spike_peristim';
output  = 'lev1_sts_convol_peristim';

% loop over the various files
nDirs = length(info);
for iDir = 1:nDirs
  
  load(fullfile(outputDir, input1, info(iDir).dataname, input1));
  load(fullfile(outputDir, input2, info(iDir).dataname, input2));
    
  % compute the spike triggered spectra
  cfg = [];
  cfg.method    = 'mtmconvol';% play with window size (200ms for freq for example)
  cfg.foi       = [20 65 100]; % can use a different set of frequencies, e.g. 2:1:160, create any resolution you want
  cfg.t_ftimwin = 11./cfg.foi; % gamma is quite narrowband, so need sufficient cycles/freq usually used 9 or 7cycles/freq play with resolution
  cfg.t_ftimwin(cfg.t_ftimwin>3) = 3; % use maximum 2 seconds
  cfg.taper     = 'hanning'; % could use other settings, e.g.:
  sts = ft_spiketriggeredspectrum(cfg, dataTrials, spikeTrials); %#ok<NASGU>
      
  % save the STS
  filename = fullfile(outputDir, output, info(iDir).dataname, output);
  mkdir(fullfile(outputDir, output, info(iDir).dataname));
  save(filename, 'sts')
   
end

