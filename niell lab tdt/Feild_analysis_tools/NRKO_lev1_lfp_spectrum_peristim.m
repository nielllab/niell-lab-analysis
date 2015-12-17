function [] = NRKO_lev1_sts_convol_peristim(overwrite, dirsel)

% NRKO_LEV1_ERP produces spike triggered spectra of the LFP

% take the globals from the info script that has been ran at the start
global info
global outputDir
input   = 'lev0_lfp_peristim';
output  = 'lev1_lfp_spectrum_peristim';

% loop over the various files
nDirs = length(info);
if nargin<2
    dirsel = 1:nDirs;
end
for iDir = dirsel
  iDir
  clear dataTrials dataStim 
  filename = fullfile(outputDir, output, info(iDir).dataname, output);
  if exist([filename '.mat'],'file') && overwrite==0, 
      fprintf('skipping %d\n', iDir); continue,
  end
  load(fullfile(outputDir, input, info(iDir).dataname, input));
     
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
  clear fr
  for iLatency = 1:2
      iLatency
      cfg = [];
      t0 = []; t_end = [];
      for k = 1:length(dataTrials.trial)
          t0(k) = dataTrials.time{k}(1);
          t_end(k) = dataTrials.time{k}(end);
      end
      stimDur = max(t_end);
      
      if hasbars
          if iLatency==1
              cfg.latency = [1 2]; % remove the onset transient from the data
          elseif iLatency==2
              cfg.latency = [min(t0)+0.2 0.5];
          end
      else
          if iLatency==1
              cfg.latency = [0.15 stimDur]; % remove the onset transient from the data
          elseif iLatency==2
              cfg.latency = [min(t0)+0.2 0];
          end
      end
      if (cfg.latency(2)-cfg.latency(1))<0.25, continue,end   
      cfg2.trials = 2:length(dataTrials.trial)-1;
      dataStim = ft_selectdata(cfg2, dataTrials);
      dataStim  = ft_selectdata(cfg, dataStim);
      for k = 1
          k
          cfg = [];  
          if k==1
              cfg.length = 0.25;
              cfg.overlap = 0.5;              
          elseif k==2
              cfg.length = 0.4;
              cfg.overlap = 0.5;
          end
          dataStim = ft_redefinetrial(cfg, dataStim);

          cfg           = [];
          cfg.foi       = 0:2:160;
          cfg.method    = 'mtmfft';
          cfg.taper     = 'hanning';
          cfg.keeptrials = 'yes';
          cfg.output     = 'fourier';
          cfg.channel = dataTrials.label(1:4:end);          
          fr(k,iLatency) = ft_freqanalysis(cfg, dataStim);
      end
  end

  % save the STS
  filename = fullfile(outputDir, output, info(iDir).dataname, output);
  mkdir(fullfile(outputDir, output, info(iDir).dataname));
  save(filename, 'fr')   
end

