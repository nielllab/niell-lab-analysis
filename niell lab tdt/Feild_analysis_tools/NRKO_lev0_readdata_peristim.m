function [] = NRKO_lev0_readdata_peristim(overwrite,dirsel)

% HOY_LEV0_READDATA_PERISTIM reads in the files as provided by JH and produces Fieldtrip
% structures around the stimulus presentations
% it also detects the movement transitions in the data and stores this information

% take the globals from the info script that has been ran at the start
global info
global outputDir
output1 = 'lev0_lfp_peristim';
output2 = 'lev0_spike_peristim';

% loop over the various files
nDirs = length(info);
if nargin<2
  dirsel = 1:nDirs;
end
if nargin<1, 
  overwrite = 0;
  warning('not overwriting the old data')
end
for iDir = dirsel
  
  iDir 
  clear dat
  % check if directory has been processed already
  filename1 = fullfile(outputDir, output1, info(iDir).dataname, output1);
  mkdir(fullfile(outputDir, output1, info(iDir).dataname));
  if exist([filename1 '.mat'], 'file') && overwrite==0, 
    fprintf('not processing file %d, %s because already exists\n', iDir, info(iDir).dataname);
    continue
  end
    
  % load the data
  load(info(iDir).datafile);
    
  % create a new structure DAT that fits with Fieldtrip
  dat = []; 
  indxChan = ~cellfun(@isempty,data.lfpT); 
  indxChan=1
  dat.time{1}   = data.lfpT{indxChan}'; 
  nChannels     = length(data.lfpT);
 
  for iChannel  = 1:nChannels
      lfpChans=[1:1:64];
    dat.trial{1}(iChannel,:) = double(data.lfpV{iChannel}); % convert to double
    dat.label{iChannel} = ['CSC_' num2str(lfpChans(iChannel))];
    %dat.label{iChannel} = ['CSC_' num2str(data.lfpChans(iChannel))];
  end
  fsample = nanmean(1./diff(dat.time{1})); % sampling rate

  % filter out the line noise; I (MV) checked the bandwidths and this seems sufficient
  % theoretically AC freq should be in 59.95-60.05 band.
  cfg = [];
  cfg.bsfilter = 'yes';
  cfg.bsfreq = [59.93 60.07];
  cfg.bsinstabilityfix = 'reduce';
  dat = ft_preprocessing(cfg, dat);
      
%   cfg = [];
%   cfg.length = 100;
%   cfg.overlap = 0.5;
%   data_long = ft_redefinetrial(cfg, dat);
%   
%   cfg = [];
%   cfg.method = 'mtmfft';
%   cfg.taper  = 'hanning';
%   cfg.output = 'pow'
%   ft = ft_freqanalysis(cfg, dataTrials);
  % detect the state point transitions on the data based on 1 cm threshold
  % and check for each trial in which state transition the trial fell
  threshold = 1.4;
  isInMove = data.mouseV{indxChan}>threshold;
  moveOn   = find(diff(isInMove)==1) + 1;
  moveOff  = find(diff(isInMove)==-1) + 1;
  
  % check if first datapoint is a move on or a move off
  moveOnTs = data.mouseT{indxChan}(moveOn);
  moveOffTs = data.mouseT{indxChan}(moveOff);
  if data.mouseV{indxChan}(1)>threshold
    moveOnTs = [data.lfpT{indxChan}(1) moveOnTs];
  else
    moveOffTs = [data.lfpT{indxChan}(1) moveOffTs];
  end    
  if moveOffTs(end)>moveOnTs(end), 
      moveOnTs = [moveOnTs data.lfpT{indxChan}(end)]; 
  else
      moveOffTs(end)<moveOnTs(end), moveOffTs = [moveOffTs data.lfpT{indxChan}(end)]; 
  end
  
  transitionPoints = [moveOffTs moveOnTs];
  transitionType   = [zeros(1,length(moveOffTs)) ones(1,length(moveOnTs))];
  [val,indx]       = sort(transitionPoints);
  transitionPoints = transitionPoints(indx);
  transitionType   = transitionType(indx);
  
  % create stimulus trials based on the timestamps that were given
  cfg = [];
  trl = data.stimEpocs{indxChan}(2,:); % this variable will contain the trial on and offsets in samples
  nTrials  = length(trl);
  samples  = [];
  
  % for each trial determine if it is in movement or sitting period, and also note how far
  [N,B] = histc(trl,transitionPoints);
  remove = B==0;
  if sum(remove>0), error('some spikes fall outside the lfp time'); end
  trl(remove) = []; B(remove) = [];
  stimulusTransitionType = transitionType(B);
  timeSinceChange        = trl-transitionPoints(B);
  timeToChange           = transitionPoints(B+1) - trl;
  
  % compute the average velocity per trial during stim and baseline as well
  for iTrial = 1:nTrials
    samples(iTrial,1) = nearest(dat.time{1}, trl(iTrial));
  end
  baseDur = median(diff(data.stimEpocs{indxChan}(2,:))) - info(iDir).stimDur;
  cfg.trl(:,1) = samples - round(fsample*baseDur);
  cfg.trl(:,2) = samples + round(fsample*info(iDir).stimDur);
  cfg.trl(:,3) = round(-baseDur*fsample);
  cfg.trl(:,4) = data.stimEpocs{indxChan}(1,:); % this will contain the trial information
  
  % remove the first trial if it was incomplete and correct baseline if <0  
  trialsToDel = [];
  if cfg.trl(1,1)<1, 
    cfg.trl(1,1) = 1; 
    cfg.trl(1,3) = cfg.trl(1,1)-samples(1); 
    if trl(1)<dat.time{1}, 
      trialsToDel = [trialsToDel 1];
    end
  end
  
  % delete the last trial if it was incomplete.
  if cfg.trl(end,2)>length(data.lfpT{indxChan}), 
    trialsToDel = [trialsToDel size(cfg.trl,1)];
  end  
  
  % put the trial info in there as well
  velocity = zeros(1,nTrials); % preallocate
  for iTrial = 1:nTrials
    indxStart = nearest(data.mouseT{indxChan},trl(iTrial));
    indxEnd   = nearest(data.mouseT{indxChan},trl(iTrial) + info(iDir).stimDur);
    velocity(iTrial) = nanmean(data.mouseV{indxChan}(indxStart:indxEnd));
  end
  
  cfg.trl(:,5) = stimulusTransitionType;
  cfg.trl(:,6) = timeSinceChange;
  cfg.trl(:,7) = timeToChange;
  cfg.trl(:,8) = velocity;
  
  % these two we don't know and have to set to NaN
  cfg.trl(1,6)   = NaN;
  cfg.trl(end,7) = NaN;  
  cfg.trl(trialsToDel,:) = [];  
  
  % construct trials in the LFP
  dataTrials = ft_redefinetrial(cfg, dat);

  % read in the unit data, create a spike structure
  nCells = length(unitdataSession);
  clear spike
  for iCell = 1:nCells
    spike.timestamp{iCell} = unitdataSession{iCell}.spikes{1};
    spike.label{iCell} = ['chan_' num2str(unitdataSession{iCell}.ch) '_cluster_' num2str(unitdataSession{iCell}.clust) ];
    cnds(iCell) = unitdataSession{iCell}.GT;
  end
  
  % make trials in the spike structure
  cfg = [];
  trl = data.stimEpocs{indxChan}(2,:);
  cfg.trl(:,1) = trl - baseDur;
  cfg.trl(:,2) = trl + info(iDir).stimDur;
  cfg.trl(:,3) = -baseDur;
  cfg.trl(:,4) = data.stimEpocs{indxChan}(1,:);
  cfg.trl(:,5) = stimulusTransitionType;
  cfg.trl(:,6) = timeSinceChange;
  cfg.trl(:,7) = timeToChange;
  cfg.trl(:,8) = velocity;  
  cfg.trl(1,6)   = NaN;
  cfg.trl(end,7) = NaN;
 
  cfg.timestampspersecond = 1;
  cfg.trl(trialsToDel,:) = [];
  spikeTrials = ft_spike_maketrials(cfg, spike); 
  if length(spikeTrials.label)~=length(unitdataSession)
      error('number of cells dont match')
  end
  trialinfo = {'move1sit0', 'time_since_change', 'time_to_change', 'velocity'};
%   cfg = [];
%   cfg.spikechannel = [1:5];
%   cfg.binsize = 0.025;
%   psth = ft_spike_psth(cfg, spikeTrials);
%    cfg = [];
%   cfg.spikechannel = [1:5];
%  
%   figure
%   cfg.topplotfunc = 'line';
%   ft_spike_plot_raster(cfg,spikeTrials,psth)
%   cfg = [];
%   jpsth = ft_spike_jpsth(cfg, psth);
%   figure
%   cfg = [];
%   cfg.channelcmb = {psth.label{1}, psth.label{2}};
%   ft_spike_plot_jpsth(cfg, jpsth)
  
  save(filename1, 'dataTrials', 'trialinfo')
  
  filename2 = fullfile(outputDir, output2, info(iDir).dataname, output2);
  mkdir(fullfile(outputDir, output2, info(iDir).dataname));  
  save(filename2, 'spikeTrials', 'spike', 'trialinfo')    
end
exit
