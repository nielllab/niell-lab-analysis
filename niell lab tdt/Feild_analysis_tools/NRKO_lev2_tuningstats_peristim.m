function [] = NRKO_lev1_sts_convol_peristim(dirsel, overwrite)

% NRKO_LEV1_STS_PERISTIM produces spike triggered spectra of the LFP

% take the globals from the info script that has been ran at the start
global info
global outputDir
input  = 'lev1_sdf_peristim';
input2  = 'lev0_spike_peristim';
output  = 'lev2_tuning_peristim';

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
     
  load(fullfile(outputDir, input, info(iDir).dataname, input));    
    load(fullfile(outputDir, input2, info(iDir).dataname, input2));    
  
  % compute measures of orientatoin tuning
  hasbar = strcmp(info(iDir).stimType, 'bars');
  if hasbar
    angles = linspace(0,360,17); angles(end) = [];
    angles = repmat(angles(:), [1 2]);
    WB     = [ones(size(angles,1),1) zeros(size(angles,1),1)];
    angles = angles'*pi/180; 
    angles = mod(angles(:),pi);
    WB     = WB';
    WB     = WB(:);

    sl        = rate(2).trialinfo(:,1)<33;
    angTrials = angles(rate(2).trialinfo(sl,1))*2;
    sm        = rate(2).trial(sl,:)'*exp(1i*angTrials);
    norm      = sum(rate(2).trial(sl,:)',2);
    osi       = abs(sm./norm); % compute the orientation selectivity index

    % see if the firing rate of the cell increases
    if ~isempty(rate(1).avg)
        modIndx = (rate(2).avg-rate(1).avg)./(rate(2).avg+rate(1).avg);
    else
        modIndx = NaN;
    end
  else
    angles = linspace(0,360,13); angles(end) = [];
    angles = repmat(angles(:), [1 6]);
    angles = angles'*pi/180; 
    angles = mod(angles(:),pi);
    
    sl        = rate(2).trialinfo(:,1)<33;
    angTrials = angles(rate(2).trialinfo(sl,1))*2;
    sm        = rate(2).trial(sl,:)'*exp(1i*angTrials);
    norm      = sum(rate(2).trial(sl,:)',2);
    osi       = abs(sm./norm); % compute the orientation selectivity index

    % see if the firing rate of the cell increases
    if ~isempty(rate(1).avg)
        modIndx = (rate(2).avg-rate(1).avg)./(rate(2).avg+rate(1).avg);
    else
        modIndx = NaN;
    end
  end
    
    % compute some locomotion signal 
    statelabels{1} = 'move';
    statelabels{2} = 'sit';
    statelabels{3} = 'early_sit'; % <10 s.
    statelabels{4} = 'late sit';
    mnRate = [];
    for iState   = 1:length(statelabels)
     for iLatency = 1:length(latencies)
          stimDur     = median(spikeTrials.trialinfo(:,2));              
          if iState==1
            cfg.trials = rate(2).trialinfo(:,2) == 1 & rate(2).trialinfo(:,4)>stimDur;
          elseif iState==2
            cfg.trials = rate(2).trialinfo(:,2) == 0 & rate(2).trialinfo(:,4)>(stimDur+2) & rate(2).trialinfo(:,3)>2;
          elseif iState==3
            cfg.trials = rate(2).trialinfo(:,2) == 0 & rate(2).trialinfo(:,4)>(stimDur+2) & rate(2).trialinfo(:,3)>2 & rate(2).trialinfo(:,3)<10;
          elseif iState==4
            cfg.trials = rate(2).trialinfo(:,2) == 0 & rate(2).trialinfo(:,4)>(stimDur+2) & rate(2).trialinfo(:,3)>20;
          end
          if ~isempty(rate(iLatency).avg)
              mnRate(iLatency,iState,:) = nanmean(rate(iLatency).trial(cfg.trials,:),1);
          else
              mnRate(iLatency,iState,:) = NaN(1,1,length(rate(2).label));
          end
     end         
    end                       
    Stat.osi = osi;
    Stat.mnRate = mnRate;
    Stat.modIndx_stim = modIndx;
    Stat.modIndx_locomotion = squeeze((mnRate(:,1,:) - mnRate(:,2,:))./(mnRate(:,1,:) + mnRate(:,2,:)));

    filename = fullfile(outputDir, output, info(iDir).dataname, output);
    mkdir(fullfile(outputDir, output, info(iDir).dataname));
    save(filename, 'Stat','-v7.3')
end

