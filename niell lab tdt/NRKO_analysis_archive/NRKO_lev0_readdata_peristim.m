function [] = Hoy_lev0_readdata_peristim

% HOY_LEV0_READDATA_PERISTIM reads in the files as provided by JH and produces Fieldtrip
% structures around the stimulus presentations

% take the globals from the info script that has been ran at the start
global info
global outputDir
output1 = 'lev0_lfp_peristim';
output2 = 'lev0_spike_peristim';

% loop over the various files
nDirs = length(info);
for iDir = 1:nDirs
  
  % load the data
  load(info(iDir).datafile);
    
  % create a new structure DAT that fits with Fieldtrip
  dat = [];
  

  dat.time{1}   = data.lfpT{1}{1}(:)'; % if this is 1,1 needs to be checked
  nChannels     = length(data.lfpT{1});
  for iChannel  = 1:nChannels
    dat.trial{1}(iChannel,:) = double(data.lfpV{1}{iChannel}); % convert to double
    dat.label{iChannel} = ['CSC_' num2str(data.lfpChans(iChannel))];
  end
  fsample = nanmean(1./diff(dat.time{1})); % sampling rate

  % create stimulus trials based on the timestamps that were given
  cfg = [];
  trl = data.stimEpocs{1}(2,:); % this variable will contain the trial on and offsets in samples
  nTrials  = length(trl);
  samples  = [];
  
  for iTrial = 1:nTrials
    samples(iTrial,1) = nearest(dat.time{1}, trl(iTrial));
  end
  baseDur = median(diff(data.stimEpocs{1}(2,:))) - info(iDir).stimDur;
  cfg.trl(:,1) = samples - round(fsample*baseDur);
  cfg.trl(:,2) = samples + round(fsample*info(iDir).stimDur);
  cfg.trl(:,3) = round(-0.1*fsample);
  cfg.trl(:,4) = data.stimEpocs{1}(1,:); % this will contain the trial information

  % remove trials at beginning or end that didn't have the full stimulus interval available.
  toDel = any(cfg.trl(:,1)<0,2) |  any(cfg.trl(:,2)>length(dat.time{1}),2); 
  cfg.trl(toDel,:) = [];  

  % construct trials in the LFP
  dataTrials = ft_redefinetrial(cfg, dat);

  % read in the unit data, create a spike structure
  nCells = length(unitdataSession);
  for iCell = 1:nCells
    spike.timestamp{iCell} = unitdataSession{iCell}.spikes{1};
    spike.label{iCell} = ['spk_' num2str(iCell)];
    cnds(iCell) = unitdataSession{iCell}.GT;
  end
 
  % make trials in the spike structure
  cfg = [];
  trl = data.stimEpocs{1}(2,:);
  cfg.trl(:,1) = trl - baseDur;
  cfg.trl(:,2) = trl + 3.0535;
  cfg.trl(:,3) = -baseDur;
  cfg.trl(:,4) = data.stimEpocs{1}(1,:);
  cfg.timestampspersecond = 1;
  cfg.trl(toDel,:) = [];
  spikeTrials = ft_spike_maketrials(cfg, spike);
  
  % save the LFPs.
  filename1 = fullfile(outputDir, output1, info(iDir).dataname, output1);
  mkdir(fullfile(outputDir, output1, info(iDir).dataname));
  save(filename1, 'dataTrials')
  
  filename2 = fullfile(outputDir, output2, info(iDir).dataname, output2);
  mkdir(fullfile(outputDir, output2, info(iDir).dataname));  
  save(filename2, 'spikeTrials', 'spike')    
end

