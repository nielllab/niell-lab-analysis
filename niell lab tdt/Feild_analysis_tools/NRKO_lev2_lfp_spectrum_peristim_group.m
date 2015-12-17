function [] = NRKO_lev3_ppc_convol_group_peristim

% HOY_LEV0_READDATA_PERISTIM reads in the files as provided by JH and produces Fieldtrip
% structures around the stimulus presentations

%% take the globals from the info script that has been ran at the start
global info
global outputDir
input2   = 'lev0_lfp_peristim';
input   =  'lev1_lfp_spectrum_peristim';
output  = 'lev2_lfp_spectrum_group_peristim';

% loop over the various files
nDirs = length(info);
latencies = {'stim', 'prestim'};
statelabels{1} = 'move';
statelabels{2} = 'sit';
statelabels{3} = 'early_sit'; % <10 s.
statelabels{4} = 'late sit';
tapers = {'short'};
layers = {'superficial', 'granular', 'deep'};
nStates  = length(statelabels);
nLatencies = length(latencies);
%%
for iDir = 1:length(info)
  iDir
  load(fullfile(outputDir, input, info(iDir).dataname, input));  
  inp2 = load(info(iDir).datafile);
  dirinf(iDir) = info(iDir).genotype;
  animalid{iDir} = inp2.data.tank;
    for iState = 1:nStates
        for iLatency = 1:nLatencies
            for iTaper = 1:length(tapers)
                try
                    struc = fr(1,iLatency);
                catch
                    continue,
                end
                if isempty(struc), continue,end
                for iLayer = 1:3
                    chans = fr(1,1).label;
                    lay = [];
                    nums = [];
                    for iChan = 1:length(chans)
                        indx = strfind(chans{iChan}, '_');
                        chanNum = str2num(chans{iChan}(indx+1:end));
                        indx = inp2.data.layer{1}(:,1)==chanNum;
                        if ~isempty(indx)
                            lay(iChan) = median(inp2.data.layer{1}(indx,2));
                        else
                            lay(iChan) = NaN;
                        end
                        nums(iChan)= chanNum;
                    end

                    % set first two to 1 and the third to 2
                    lay(isnan(lay) & (nums==1 | nums==5 | nums==33 | nums==37)) = 1;                                                        
                    if isnan(lay(3))
                        lay(3) = 2; 
                    end
                    if isnan(lay(11))
                        lay(11) = 2;
                    end

                    % fill in the info from the other shank
                    indx = find(isnan(lay));                            
                    if ~isempty(indx)
                       for k = 1:length(indx)
                           if indx(k)<9
                               lay(indx(k)) = lay(indx(k)+8);
                           else
                               lay(indx(k)) = lay(indx(k)-8);
                           end
                       end
                    end                     

                    % check if we are missing superficial channels                            
                    indxLast1 = find(lay(1:8)==3|lay(1:8)==2,1,'last'); % determine the last L23, shank 1
                    indxLast2 = find(lay(9:16)==3|lay(9:16)==2,1,'last');                                                                                                               
                    indx = find(isnan(lay));
                    if ~isempty(indx)
                        for k = 1:length(indx)
                            if ~isempty(indxLast1) & indx(k)<indxLast1
                                lay(indx(k)) = 2;
                            elseif ~isempty(indxLast2) & (indx(k)>8 & indx(k)<(indxLast2+8))
                                lay(indx(k)) = 2;
                            end
                        end
                    end

                    
                     stimDur     = info(iDir).stimDur;
                      if iState==1
                        cfg.trials = struc.trialinfo(:,2) == 1 & struc.trialinfo(:,4)>stimDur;
                      elseif iState==2
                        cfg.trials = struc.trialinfo(:,2) == 0 & struc.trialinfo(:,4)>(stimDur+2) & struc.trialinfo(:,3)>2;
                      elseif iState==3
                        cfg.trials = struc.trialinfo(:,2) == 0 & struc.trialinfo(:,4)>(stimDur+2) & struc.trialinfo(:,3)>2 ...
                            & struc.trialinfo(:,3)<10;
                      elseif iState==4
                        cfg.trials = struc.trialinfo(:,2) == 0 & struc.trialinfo(:,4)>(stimDur+2) & struc.trialinfo(:,3)>20;
                      end

                    if iLayer==1
                        if any(lay>1 & lay<4)
                            powAllCat(iState,iLatency,iLayer).powAll(iDir,:) = nanmean(nanmean(abs(struc.fourierspctrm(cfg.trials,lay>1 & lay<4,:)),1),2);
                        end
                    elseif iLayer==2
                        if any(lay==4)
                            powAllCat(iState,iLatency,iLayer).powAll(iDir,:) = nanmean(nanmean(abs(struc.fourierspctrm(cfg.trials,lay==4,:)),1),2);
                        end
                    elseif iLayer==3
                        if any(lay==5)
                            powAllCat(iState,iLatency,iLayer).powAll(iDir,:) = nanmean(nanmean(abs(struc.fourierspctrm(cfg.trials,lay==5,:)),1),2);
                        end                
                    end                                                                
            end
        end
    end
  end
end
%%
G_stat.dirinfo = dirinf;
G_stat.animalinfo = animalid;
G_stat.powAllCat = powAllCat;
G_stat.dimord = 'state_latency_layer';    
G_stat.latencies = latencies;
G_stat.state = statelabels;
G_stat.layers = layers;
filename = fullfile(outputDir,output,output);
mkdir(fullfile(outputDir, output))
save(filename,'G_stat');
%G_stat.unitinfo(G_stat.unitinfo(:,4)==2,4) = 3;
% %%
% for iLatency= 1:2
%     cors = {'b', 'r'};
%     iState = 1
%     iStim = 2
%     for iLayer = 4%3:5
%         for iType = 0
%             figure
%             for k = 0:1
%                 try                
%                     sl =  G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==k & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
%                     mn = squeeze(nanmean(G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).ppcAll(sl,:,:),1));
%                     sm = squeeze(nanstd(G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).ppcAll(sl,:,:),1))./sqrt(sum(sl));
%                     freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).freq;
%                     errorbar(freq, mn,sm,cors{k+1})
%                     hold on
%                     ylim([0 0.04])
%                     title([num2str(iLayer) '_' num2str(iType)])
%                 catch
%                 end
%             end
%         end
%     end
% end
% %%
% G_stat.unitinfo(G_stat.unitinfo(:,4)==2,4) = 3;
% cors = {'b', 'r'};
% iStim = 2
% for iLayer = 3:5
%     for iType = 0:2
%         figure
%         for k = 0:1
%             try                
%                 sl = G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==k & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
%                 mn = squeeze(nanmean(G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).ppcAll(sl,:,:),1));
%                 sm = squeeze(nanstd(G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).ppcAll(sl,:,:),1))./sqrt(sum(sl));
%                 freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).freq;
%                 errorbar(freq, mn,sm,cors{k+1})
%                 hold on
%                 ylim([0 0.03])
%                 title([num2str(iLayer) '_' num2str(iType)])
%             catch
%             end
%         end
%     end
% end
% 
% %%
% iStim = 2;
% cors = {'b', 'g'};
% for iType = 0:2
%     figure
%     for k = 0:1
%         try
%             if k==0
%                 sl = G_stat.unitinfo(:,4)==4 & G_stat.unitinfo(:,2)==k & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==0 & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
%             else
%                 sl = G_stat.unitinfo(:,2)==k & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==0 & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
%             end
%             mn = squeeze(nanmean(G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).ppcAll(sl,:,:),1));
%             sm = squeeze(nanstd(G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).ppcAll(sl,:,:),1))./sqrt(sum(sl));
%             freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).freq;
%             errorbar(freq, mn,sm,cors{k+1})
%             hold on
%             ylim([0 0.03])
%         catch
%         end
%     end
% end
% %%
% % nLatencies (base - stim)
% % run sit
% % layer unit
% % layer LFP
% % genotype
% % interneuron pyramidal
% % l4 cre or no cre
% % 2 x 2 x 3 x 3 x 3 x 2 = 216 conditions 
% %2 x 2 x 3 x 3 x 3 x 2 = 216 conditions
% 
% 
% 
% 








