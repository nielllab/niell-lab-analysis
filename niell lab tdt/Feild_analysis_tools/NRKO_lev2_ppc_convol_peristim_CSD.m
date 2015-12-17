function [] = NRKO_lev2_ppc_convol_peristim(overwrite, dirsel)

% take the globals from the info script that has been ran at the start
global info
global outputDir
input   = 'lev1_sts_convol_peristim_new';
input2  = 'lev0_spike_peristim';
output  = 'lev2_ppc_convol_peristim_CSD';

% loop over the various files
nDirs = length(info);

if nargin<2
    dirsel = 1:nDirs;
end
    
for iDir = dirsel(:)'
  clear sts_long sts_short sts_very_short
  filename = fullfile(outputDir, output, info(iDir).dataname, output);
  if exist([filename,'.mat'],'file'), fprintf('skipping %d',iDir); continue,end
  
  load(fullfile(outputDir, input, info(iDir).dataname, input));  
  load(fullfile(outputDir, input2, info(iDir).dataname, input2));

  iDir
  
  % different ppc versions, angle, rayleigh test
  methods  = {'ppc1', 'ang', 'ral'};%usually used ppc1 most conservative/sesnsitive 
  nMethods = length(methods); 
  clear stat stat_tfr   
  statelabels{1} = 'move';
  statelabels{2} = 'sit';
  try
      if strcmp(info(iDir).stimType,'drift')                          
          latencies{1} = [min(spikeTrials.trialtime(:,1))+0.2 0];
          latencies{2} = [0.06 max(spikeTrials.trialtime(:,2))];
      else
          latencies{1} = [min(spikeTrials.trialtime(:,1)) 0.5];
          latencies{2} = [1 2];
      end
  catch
     if strcmp(info(iDir).stimtype,'drift')                          
          latencies{1} = [min(spikeTrials.trialtime(:,1))+0.2 0];
          latencies{2} = [0.06 max(spikeTrials.trialtime(:,2))];
      else
          latencies{1} = [min(spikeTrials.trialtime(:,1)) 0.5];
          latencies{2} = [1 2];
      end
  end             
                    
  % compute the CSDs, just the second derivative, difference of two
  % channels minus the difference between the other two channels
  for k = 2:7
      for iUnit = 1:length(sts_long.fourierspctrm)
          try
            sts_long.fourierspctrm{iUnit}(:,k,:) = (sts_long.fourierspctrm{iUnit}(:,k+1,:)-sts_long.fourierspctrm{iUnit}(:,k,:))-(sts_long.fourierspctrm{iUnit}(:,k,:)-sts_long.fourierspctrm{iUnit}(:,k-1,:))
            sts_short.fourierspctrm{iUnit}(:,k,:) = (sts_short.fourierspctrm{iUnit}(:,k+1,:)-sts_short.fourierspctrm{iUnit}(:,k,:))-(sts_short.fourierspctrm{iUnit}(:,k,:)-sts_short.fourierspctrm{iUnit}(:,k-1,:))
            sts_very_short.fourierspctrm{iUnit}(:,k,:) = ( sts_very_short.fourierspctrm{iUnit}(:,k+1,:)- sts_very_short.fourierspctrm{iUnit}(:,k,:))-( sts_very_short.fourierspctrm{iUnit}(:,k,:)- sts_very_short.fourierspctrm{iUnit}(:,k-1,:))
          catch
              disp('error with CSD, no spikes possible')              
              iUnit
          end
      end
  end
  
  for k = 9:15
      for iUnit = 1:length(sts_long.fourierspctrm)
          try
            sts_long.fourierspctrm{iUnit}(:,k,:) = (sts_long.fourierspctrm{iUnit}(:,k+1,:)-sts_long.fourierspctrm{iUnit}(:,k,:))-(sts_long.fourierspctrm{iUnit}(:,k,:)-sts_long.fourierspctrm{iUnit}(:,k-1,:))
            sts_short.fourierspctrm{iUnit}(:,k,:) = (sts_short.fourierspctrm{iUnit}(:,k+1,:)-sts_short.fourierspctrm{iUnit}(:,k,:))-(sts_short.fourierspctrm{iUnit}(:,k,:)-sts_short.fourierspctrm{iUnit}(:,k-1,:))
            sts_very_short.fourierspctrm{iUnit}(:,k,:) = ( sts_very_short.fourierspctrm{iUnit}(:,k+1,:)- sts_very_short.fourierspctrm{iUnit}(:,k,:))-( sts_very_short.fourierspctrm{iUnit}(:,k,:)- sts_very_short.fourierspctrm{iUnit}(:,k-1,:))
          catch
              disp('error with CSD')
              iUnit
          end
      end
  end
      
  for iTaper = 1:2
      for iState   = 1:length(statelabels)
          for iLatency = 1:length(latencies)
              for iMethod = 1%:nMethods
                nUnits = length(sts_long.label);
                unitCnt = 0;
                for iUnit = 1:nUnits      
                  unitCnt = unitCnt + 1;
                  if size(sts_long.fourierspctrm{iUnit},1)<10, continue,end % reject less than 10 spikes in total          
                  cfg = [];
                  cfg.spikechannel = iUnit;
                  cfg.method = methods{iMethod};      
                  cfg.latency = latencies{iLatency};
                  stimDur     = median(spikeTrials.trialinfo(:,2));
                  if iState==1
                    cfg.trials = spikeTrials.trialinfo(:,2) == 1 & spikeTrials.trialinfo(:,4)>stimDur;
                  elseif iState==2
                    cfg.trials = spikeTrials.trialinfo(:,2) == 0 & spikeTrials.trialinfo(:,4)>(stimDur+2) & spikeTrials.trialinfo(:,3)>2;
                  elseif iState==3
                    cfg.trials = spikeTrials.trialinfo(:,2) == 0 & spikeTrials.trialinfo(:,4)>(stimDur+2) & spikeTrials.trialinfo(:,3)>2 & spikeTrials.trialinfo(:,3)<10;
                  elseif iState==4
                    cfg.trials = spikeTrials.trialinfo(:,2) == 0 & spikeTrials.trialinfo(:,4)>(stimDur+2) & spikeTrials.trialinfo(:,3)>20;
                  end
                    
                  try
                      iDir
                      if iTaper==1
                          stat(iMethod,iState,iLatency,iTaper).static(unitCnt) = ft_spiketriggeredspectrum_stat(cfg, sts_long);
                      elseif iTaper==2
                          stat(iMethod,iState,iLatency,iTaper).static(unitCnt) = ft_spiketriggeredspectrum_stat(cfg, sts_short);
                      else
                          stat(iMethod,iState,iLatency,iTaper).static(unitCnt) = ft_spiketriggeredspectrum_stat(cfg, sts_very_short);
                      end              
                  catch
                      S = lasterror;
                      S.stack
                      fprintf('not enough spikes in unit %d\n', iUnit);
                  end
                end      
              end
          end
      end  
  end
  
  % save the STS
  filename = fullfile(outputDir, output, info(iDir).dataname, output);
  mkdir(fullfile(outputDir, output, info(iDir).dataname));
  %if exist ('stat', 'var')
  save(filename, 'stat', 'statelabels', 'latencies', 'methods')
  %end
end