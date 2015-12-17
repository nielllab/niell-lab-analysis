function [] = NRKO_lev2_ppc_convol_peristim(overwrite, dirsel)

%% take the globals from the info script that has been ran at the start
global info
global outputDir
input   = 'lev1_lfp_spectrum_peristim';
input2  = 'lev0_lfp_peristim';
output  = 'lev2_lfp_spectrum_peristim';

% loop over the various files
nDirs = length(info);

% if nargin<2, dirsel = 1:nDirs; end
% if nargin<1, overwrite = 0; end
%  cfg = [];
%  cfg.timwin = [-0.025 0.025];
%  [sdf,sdfdata] = ft_spikedensity(cfg, spikeTrials);
%  
%  cfg = [];
%  erp = ft_timelockanalysis(cfg, dataTrials);

% psth = ft_spike_psth(cfg,spikeTrials);
%  for k = 1:25
%      cfg = [];
%      cfg.spikechannel = spikeTrials.label{k}
%      cfg.topplotfunc = 'line'
%      figure, ft_spike_plot_raster(cfg,spikeTrials,sdf)
%  end
    
if nargin<2
    dirsel = 1:nDirs;
end
%%    
for iDir = dirsel(:)'
  load(fullfile(outputDir, input, info(iDir).dataname, input));  
  load(fullfile(outputDir, input2, info(iDir).dataname, input2));
  keyboard
  iDir
  
  % different ppc versions, angle, rayleigh test
   clear stat stat_tfr   
  statelabels{1} = 'move';
  statelabels{2} = 'sit';
  statelabels{3} = 'early_sit'; % <10 s.
  statelabels{4} = 'late sit';
  
  for iTaper = 1:2      
      for iState   = 1:length(statelabels)
          for iLatency = 1:3%length(latencies)
              for iMethod = 1:nMethods
                nChans = length(fr(1,1).label);
                unitCnt = 0;
                for iChan = 1:nChans
                  unitCnt = unitCnt + 1;
                  cfg = [];
                  cfg.latency = latencies{iLatency};
                  keyboard
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
                      keyboard
%                       if iTaper==1
%                           stat(iMethod,iState,iLatency,iTaper).static(unitCnt) = ft_spiketriggeredspectrum_stat(cfg, sts_long);
%                       elseif iTaper==2
%                           stat(iMethod,iState,iLatency,iTaper).static(unitCnt) = ft_spiketriggeredspectrum_stat(cfg, sts_short);
%                       else
%                           stat(iMethod,iState,iLatency,iTaper).static(unitCnt) = ft_spiketriggeredspectrum_stat(cfg, sts_very_short);
%                       end              
                  catch
%                      S = lasterror;
%                      S.stack
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