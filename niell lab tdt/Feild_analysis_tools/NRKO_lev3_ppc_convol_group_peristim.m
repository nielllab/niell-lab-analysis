function [] = NRKO_lev3_ppc_convol_group_peristim

% HOY_LEV0_READDATA_PERISTIM reads in the files as provided by JH and produces Fieldtrip
% structures around the stimulus presentations

%% take the globals from the info script that has been ran at the start
global info
global outputDir
input   = 'lev2_ppc_convol_peristim';
input2  = 'lev0_spike_peristim';
output  = 'lev3_ppc_convol_group_peristim';

% loop over the various files
nDirs = length(info);
methods   = {'ppc1'};
tapers    = {'long'};
latencies = {'prestim', 'stim'};
statelabels{1} = 'move';
statelabels{2} = 'sit';
statelabels{3} = 'early_sit'; % <10 s.
statelabels{4} = 'late sit';
layers = {'all', 'superficial', 'granular', 'deep', 'all_diff', 'superficial_diff', 'granular_diff'};
nMethods = length(methods); 
nStates  = length(statelabels);
nLatencies = length(latencies);
nTapers    = length(tapers);
[dofAllCat, ppcAllCat] = deal(cell(1,nMethods));
genotype = []; animalid = [];
dimord = 'method_state_latency_taper';
%%
totUnits = 0;
clear ppcAllCat
for iDir = 1:nDirs
  iDir
  clear dat data
  load(fullfile(outputDir, input, info(iDir).dataname, input));  
  inp2 = load(info(iDir).datafile);
  
  for iMethod = 1:nMethods
    for iState = 1:nStates
        for iLatency = 1:nLatencies
            for iTaper = 1:nTapers
                nUnits = length(inp2.unitdataSession);
                if size(stat,2)<iState, continue,end 
                pcStruc = stat(iMethod, iState, iLatency, iTaper).static;
                for iUnit = 1:nUnits
                    unitCnt =  totUnits + iUnit;                   
                    unitinfo(unitCnt,1) = inp2.unitdataSession{iUnit}.GT*0+info(iDir).genotype;
                    %unitinfo(unitCnt,2) = inp2.unitdataSession{iUnit}.pinp;
                    unitinfo(unitCnt,2) = inp2.unitdataSession{iUnit}.responsive;
                    unitinfo(unitCnt,3) = inp2.unitdataSession{iUnit}.layer;
                    unitinfo(unitCnt,4) = inp2.unitdataSession{iUnit}.inhibitory;             
                    try
                        if strcmp(info(iDir).stimtype,'drift')
                            unitinfo(unitCnt,6) = 1;  
                        else
                            unitinfo(unitCnt,6) = 2;
                        end
                    catch
                        if strcmp(info(iDir).stimType,'drift')
                            unitinfo(unitCnt,6) = 1;  
                        else
                            unitinfo(unitCnt,6) = 2;
                        end
                    end
                    
                    animalinfo{unitCnt} = inp2.data.tank;
                                        
                    if iMethod==1 && length(pcStruc)>=iUnit
                        coh = pcStruc(iUnit).ppc1;
                        nspikes = pcStruc(iUnit).nspikes;                                      
                    else
                        coh = []; 
                        nspikes = 0; 
                    end
                    if ~isempty(coh)
                        % loop over the layers
                        for iLayer = 1:7
                            % identify and reject the same channe
                            unitChan = pcStruc(iUnit).labelcmb{1,1};
                            indx = strfind(unitChan,'_');
                            unitChanNum = unitChan(indx(1)+1:indx(2)-1);
                            sameChan = ~cellfun(@isempty,strfind(pcStruc(iUnit).labelcmb(:,2),unitChanNum));
                           
                            % extract the layer for each electrode
                            chans = pcStruc(iUnit).labelcmb(:,2);
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
                            if unitChanNum<33
                                sameProbe = nums<33;
                            else
                                sameProbe = nums>=33;
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
                            sameChan = sameChan(:); sameProbe = sameProbe(:); lay = lay(:); 
                            
                            if iLayer==1
                                if any(~sameChan & sameProbe)
                                    ppcAllCat(iMethod,iState,iLatency,iTaper,iLayer).ppcAll(unitCnt,:) = nanmean(pcStruc(iUnit).ppc1(~sameChan & sameProbe,:),1);
                                end
                            elseif iLayer==2
                                if any(~sameChan & lay>1 & lay<4 & sameProbe)
                                    ppcAllCat(iMethod,iState,iLatency,iTaper,iLayer).ppcAll(unitCnt,:) = nanmean(pcStruc(iUnit).ppc1(~sameChan & lay>1 & lay<4 & sameProbe,:),1);
                                end
                            elseif iLayer==3
                                if any(~sameChan & lay==4 & sameProbe)
                                    ppcAllCat(iMethod,iState,iLatency,iTaper,iLayer).ppcAll(unitCnt,:) = nanmean(pcStruc(iUnit).ppc1(~sameChan & lay==4 & sameProbe,:),1);
                                end
                            elseif iLayer==4
                                if any(~sameChan & lay==5)
                                    ppcAllCat(iMethod,iState,iLatency,iTaper,iLayer).ppcAll(unitCnt,:) = nanmean(pcStruc(iUnit).ppc1(~sameChan & lay==5,:),1);
                                end
                            elseif iLayer==5
                                if any(~sameChan & ~sameProbe)
                                    ppcAllCat(iMethod,iState,iLatency,iTaper,iLayer).ppcAll(unitCnt,:) = nanmean(pcStruc(iUnit).ppc1(~sameChan & ~sameProbe,:),1);
                                end
                            elseif iLayer==6
                                if any(~sameChan & lay>1 & lay<4 & ~sameProbe)
                                    ppcAllCat(iMethod,iState,iLatency,iTaper,iLayer).ppcAll(unitCnt,:) = nanmean(pcStruc(iUnit).ppc1(~sameChan & lay>1 & lay<4 & ~sameProbe,:),1);
                                end
                            elseif iLayer==7
                                if any(~sameChan & lay==4 & ~sameProbe)
                                    ppcAllCat(iMethod,iState,iLatency,iTaper,iLayer).ppcAll(unitCnt,:) = nanmean(pcStruc(iUnit).ppc1(~sameChan & lay==4 & ~sameProbe,:),1);
                                end
                            end                            
                                
                            ppcAllCat(iMethod,iState,iLatency,iTaper,iLayer).freq = pcStruc(iUnit).freq;
                            ppcAllCat(iMethod,iState,iLatency,iTaper,iLayer).nSpikes(unitCnt) = max(nspikes(:));                       
                        end
                    else
                        for iLayer = 1:7
                            ppcAllCat(iMethod,iState,iLatency,iTaper,iLayer).nSpikes(unitCnt) = 0;                                                 
                        end
                    end
                end
            end
        end
    end
  end
  totUnits = totUnits + nUnits;
  size(ppcAllCat)
end
%%
G_stat.unitinfo = unitinfo;
G_stat.animalinfo = animalinfo;
G_stat.ppcAllCat = ppcAllCat;
G_stat.dimord = 'method_state_latency_taper_layer';    
G_stat.methods = methods;
G_stat.tapers = tapers;
G_stat.latencies = latencies;
G_stat.state = statelabels;
G_stat.layers = layers;
%G_stat.unitinfo_dim = {'wt_N2B_N2A', 'pinp_nopinp', 'responsive', 'layer', 'inhibitory', 'drift1_bar2'};
G_stat.unitinfo_dim = {'stage', 'responsive', 'layer', 'inhibitory', 'stimulus type'};

filename = fullfile(outputDir,output,output);
mkdir(fullfile(outputDir,output))
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








