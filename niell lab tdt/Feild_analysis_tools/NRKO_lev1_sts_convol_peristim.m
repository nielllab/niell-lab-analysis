function [] = NRKO_lev1_sts_convol_peristim(overwrite,dirsel)
% NRKO_LEV1_STS_PERISTIM produces spike triggered spectra of the LFP

% take the globals from the info script that has been ran at the start
global info
global outputDir
input1  = 'lev0_lfp_peristim';
input2  = 'lev0_spike_peristim';
output  = 'lev1_sts_convol_peristim_new';

% loop over the various files
nDirs = length(info);
%if nargin<1, dirsel = 1:nDirs; end

if nargin<2
  dirsel = 1:nDirs;
end
if nargin<1, 
  overwrite = 0;
  warning('not overwriting the old data')
end

for iDir = dirsel(:)'
  
    iDir

   filename = fullfile(outputDir, output, info(iDir).dataname, output);
  if exist([filename '.mat'], 'file') && overwrite==0, 
      fprintf('skipping directory %d because has already been done\n', iDir), 
      continue
  end
  
    load(fullfile(outputDir, input1, info(iDir).dataname, input1));
  load(fullfile(outputDir, input2, info(iDir).dataname, input2));
    
    
  if isfield(info(iDir),'stimType')
      hasbars = strcmp(info(iDir).stimType,'bars');
  else
      hasbars = strcmp(info(iDir).stimtype,'bars');
  end
  
  % remove the ERP if it does have drifting gratings
  if ~hasbars 
      cfg = [];
      cfg.vartrllength = 2;
      erp = ft_timelockanalysis(cfg, dataTrials);
      for k = 1:length(dataTrials.trial)
          try
              dataTrials.trial{k} = dataTrials.trial{k} - erp.avg;
          catch
              disp('trial has a different length'); 
          end
      end
  end
  
      
%  cfg = [];
%    cfg.hpfilter = 'yes';
%    cfg.hpfreq   = [30];
%    cfg.hpinstabilityfix = 'reduce'
%    dataFilt = ft_preprocessing(cfg, dataTrials);
%   
%      cfg = [];
%      cfg.keeptrials = 'no';
%      cfg.timwin     = [-0.1 0.1];
%      sta = ft_spiketriggeredspectrum_average2(cfg, dataTrials,spikeTrials);
%           
%     for k = 1:9
%         figure, plot(sta.time, squeeze(sta.avg(1,k,:)))
%     end
  % compute the spike triggered spectra
  cfg = [];
  cfg.channel = dataTrials.label(1:4:end);
  cfg.method    = 'mtmconvol';% play with window size (200ms for freq for example)
  cfg.foi       = [4:2:140]; % can use a different set of frequencies, e.g. 2:1:160, create any resolution you want
  cfg.t_ftimwin = 11./cfg.foi; % gamma is quite narrowband, so need sufficient cycles/freq usually used 9 or 7cycles/freq play with resolution
  cfg.t_ftimwin(cfg.t_ftimwin>1.5) = 1.5; % use maximum 2 seconds
  cfg.taper     = 'hanning'; % could use other settings, e.g.:
  sts_long = ft_spiketriggeredspectrum(cfg, dataTrials, spikeTrials); %#ok<NASGU>
%   cfg = [];
%   cfg.method = 'ppc1';
%   cfg.trials = spikeTrials.trialinfo(:,2)==1;
%   cfg.latency = [1 2]
%   pc = ft_spiketriggeredspectrum_stat(cfg, sts_long);
  % compute the spike triggered averages of the data
   % compute the spike triggered spectra
%   cfg = [];
%   cfg.channel = dataTrials.label(1:4:end);
%   cfg.method    = 'mtmconvol';% play with window size (200ms for freq for example)
%   cfg.foi       = [4:2:140]; % can use a different set of frequencies, e.g. 2:1:160, create any resolution you want
%   cfg.t_ftimwin = 7./cfg.foi; % gamma is quite narrowband, so need sufficient cycles/freq usually used 9 or 7cycles/freq play with resolution
%   cfg.t_ftimwin(cfg.t_ftimwin>1.5) = 1.5; % use maximum 2 seconds
%   cfg.taper     = 'hanning'; % could use other settings, e.g.:
%   sts_short = ft_spiketriggeredspectrum(cfg, dataTrials, spikeTrials); %#ok<NASGU>
% 
%   % compute the spike field coherence of the data
%     % compute the spike triggered averages of the data
%    % compute the spike triggered spectra
%   cfg = [];
%   cfg.channel = dataTrials.label(1:4:end);
%   cfg.method    = 'mtmconvol';% play with window size (200ms for freq for example)
%   cfg.foi       = [4:2:140]; % can use a different set of frequencies, e.g. 2:1:160, create any resolution you want
%   cfg.t_ftimwin = 3./cfg.foi; % gamma is quite narrowband, so need sufficient cycles/freq usually used 9 or 7cycles/freq play with resolution
%   cfg.t_ftimwin(cfg.t_ftimwin>1.5) = 1.5; % use maximum 2 seconds
%   cfg.taper     = 'hanning'; % could use other settings, e.g.:
%   sts_very_short = ft_spiketriggeredspectrum(cfg, dataTrials, spikeTrials); %#ok<NASGU>

  % save the STS
  filename = fullfile(outputDir, output, info(iDir).dataname, output);
  mkdir(fullfile(outputDir, output, info(iDir).dataname));
  %save(filename, 'sts_long','sts_short', 'sts_very_short','-v7.3')
   save(filename, 'sts_long','-v7.3')
   
end
exit
