function [] = NRKO_lev1_sts_convol_peristim(dirsel, overwrite)

% NRKO_LEV1_STS_PERISTIM produces spike triggered spectra of the LFP

% take the globals from the info script that has been ran at the start
global info
global outputDir
input  = 'lev0_spike_peristim';
output  = 'lev1_sdf_peristim';

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
  
  binsizes = [0.01 0.025 0.05];
  
 % compute the psth
 clear psth
 for iBinsizes = 1:length(binsizes)
      cfg.binsize = 0.025;
      psth(iBinsizes) = ft_spike_psth(cfg, spikeTrials);
 end
  
  if strcmp(info(iDir).stimType,'drift')          
      if min(spikeTrials.trialtime(:,1))<-0.2
          latencies{1} = [min(spikeTrials.trialtime(:,1))+0.2 0];
      else
          latencies{1} = [min(spikeTrials.trialtime(:,1)) 0];
      end
      latencies{2} = [0.06 max(spikeTrials.trialtime(:,2))];
  else
      latencies{1} = [min(spikeTrials.trialtime(:,1)) 0.5];
      latencies{2} = [1 2];
  end
  
  clear rate
  for iLatency = 1:2
      cfg = [];
      cfg.latency = latencies{iLatency};
      if cfg.latency(2)<cfg.latency(1), continue,end
      rate(iLatency) = ft_spike_rate(cfg, spikeTrials);
  end
  
  filename = fullfile(outputDir, output, info(iDir).dataname, output);
  mkdir(fullfile(outputDir, output, info(iDir).dataname));
  save(filename, 'latencies', 'rate','psth', '-v7.3')
end

