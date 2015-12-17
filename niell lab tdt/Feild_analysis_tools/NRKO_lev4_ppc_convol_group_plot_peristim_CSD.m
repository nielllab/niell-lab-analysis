function [] = NRKO_lev3_ppc_convol_group_peristim

global info
global outputDir
input   = 'lev2_ppc_convol_peristim';
input2  = 'lev0_spike_peristim';
output  = 'lev3_ppc_convol_group_peristim';
filename = fullfile(outputDir,output,output);
mkdir(fullfile(outputDir,output))
load(filename,'G_stat');
%%

G_stat.unitinfo(G_stat.unitinfo(:,4)==2,4) = 3;
%% wild type
% first figure: relative to different layers, with running and with sitting
iType = 1
iMethod = 1; % phsae locking
iTaper = 1;
for iStim = 1:2
    for iState = 1:2
            figure
            cnt = 0;
            for iLayer = 3:5
                for k = 0:1
                   for jLayer = [2 3 4]

                        cors = {'k', 'g'};
                        cnt=cnt+1;
                        subplot(6,3,cnt)
                        for iLatency= 1:2
                            hold on
                            try                                            
                                sl =  G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==k & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                dat = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
                                dat(dat==0) = NaN;
                                dat(~isfinite(dat)) = NaN;

                                mn = squeeze(nanmean(dat,1));
                                sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                                freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).freq;
                                errorbar(freq, mn,sm,cors{iLatency})
                                hold on
                                ylim([0 0.05])
                                latencies = {'prestim', 'stim'};
                                stims = {'drift', 'bar'};
                                layers = {'nan', 'nan','superficial', 'granular', 'deep'};
                                nt = {'pyr', 'int'};
                                H = title(sprintf('unit: %s LFP: %s %s %s %s df=%d',  layers{iLayer}, G_stat.layers{jLayer}, G_stat.state{iState}, stims{iStim}, nt{k+1},sum(sl)));
                                set(H,'FontSize', 6);
                            end                        
                        end
                    end
                end
        end
    end
end
%% SIT VS RUN WITHIN EACH PLOT
iType = 1
iMethod = 1; % phsae locking
iTaper = 1;
for iStim = 1:2
    for iLatency = 1:2
            figure
            cnt = 0;
            for iLayer = 3:5
                for k = 0:1
                   for jLayer = [2 3 4]

                        cors = {'g', 'r', 'c', 'k'};
                        cnt=cnt+1;
                        subplot(6,3,cnt)
                        for iState= [1:2]
                            hold on
                            try                                            
                                sl =  G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==k & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                dat = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
                                dat(dat==0) = NaN;
                                dat(~isfinite(dat)) = NaN;

                                mn = squeeze(nanmean(dat,1));
                                sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                                freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).freq;
                                errorbar(freq, mn,sm,cors{iState})
                                hold on
                                ylim([0 0.05])
                                latencies = {'prestim', 'stim'};
                                stims = {'drift', 'bar'};
                                layers = {'nan', 'nan','superficial', 'granular', 'deep'};
%                                title(sprintf('Layer unit: %s \n  Layer LFP: %s',  layers{iLayer}, G_stat.layers{jLayer}));      nt = {'pyr', 'int'};
                                H = title(sprintf('unit: %s LFP: %s %s %s %s df=%d',  layers{iLayer}, G_stat.layers{jLayer}, latencies{iLatency}, stims{iStim}, nt{k+1},sum(sl)));
                                set(H,'FontSize', 6);                 
                            end                        
                        end
                    end
                end
        end
    end
end
%
%% wild type
% INTERNEURONS VS PYRAMIDS
iType = 1
iMethod = 1; % phsae locking
iTaper = 1;
for iStim = 1:2
    for iState = 1:2
            figure
            cnt = 0;
            for iLayer = 3:5
                for iLatency = 1:2
                   for jLayer = [2 3 4]

                        cors = {'b', 'r'};
                        cnt=cnt+1;
                        subplot(6,3,cnt)
                        for k = 0:1
                            hold on
                            try                                            
                                sl =  G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==k & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                dat = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
                                dat(dat==0) = NaN;
                                dat(~isfinite(dat)) = NaN;

                                mn = squeeze(nanmean(dat,1));
                                sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                                freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).freq;
                                errorbar(freq, mn,sm,cors{k+1})
                                hold on
                                ylim([0 0.05])
                                latencies = {'prestim', 'stim'};
                                stims = {'drift', 'bar'};
                                layers = {'nan', 'nan','superficial', 'granular', 'deep'};
                                %title(sprintf('Layer unit: %s \n  Layer LFP: %s',  layers{iLayer}, G_stat.layers{jLayer}));
                                    H = title(sprintf('unit: %s LFP: %s %s %s %s df=%d',  layers{iLayer}, G_stat.layers{jLayer}, latencies{iLatency}, stims{iStim}, nt{k+1},sum(sl)));
                                set(H,'FontSize', 6);                 
                  
                            end                        
                        end
                    end
                end
        end
    end
end
%%
iType = 1
iMethod = 1; % phsae locking
iTaper = 1;
for iLatency = 1:2
    for iState = 1:2
            figure
            cnt = 0;
            for iLayer = 3:5
                for k = 0:1
                   for jLayer = [2 3 4]

                        cors = {'k', 'g'};
                        cnt=cnt+1;
                        subplot(6,3,cnt)
                        for iStim= 1:2
                            hold on
                            try                                            
                                sl =  G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==k & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                dat = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
                                dat(dat==0) = NaN;
                                dat(~isfinite(dat)) = NaN;

                                mn = squeeze(nanmean(dat,1));
                                sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                                freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).freq;
                                errorbar(freq, mn,sm,cors{iStim})
                                hold on
                                ylim([0 0.05])
                                latencies = {'prestim', 'stim'};
                                stims = {'drift', 'bar'};
                                layers = {'nan', 'nan','superficial', 'granular', 'deep'};
                                nt = {'pyr', 'int'};
                                H = title(sprintf('unit: %s  LFP: %s, %s, %s, %s',  layers{iLayer}, G_stat.layers{jLayer}, latencies{iLatency}, G_stat.state{iState}, nt{k+1}));
                                set(H,'FontSize', 6)
                            end                        
                        end
                    end
                end
        end
    end
end
%%

iType = 1
iMethod = 1; % phsae locking
iTaper = 1;
for iLatency = 1:2
    for iState = 1:2
            figure
            cnt = 0;
            for iLayer = 3:5
                for k = 0:1
                   for jLayer = [2 3 4]

                        cors = {'k', 'g'};
                        cnt=cnt+1;
                        subplot(6,3,cnt)
                        for iStim= 1:2
                            hold on
                            try                                            
                                sl =  G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==k & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                dat = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
                                dat(dat==0) = NaN;
                                dat(~isfinite(dat)) = NaN;

                                [A,B,C] = unique(G_stat.animalinfo);
                                cfg.semiweighted_adaptive = 'no';
                                for iFrequency = 1:size(dat,2)                                    
                                    stt = Lur_statistics(cfg, dat(:,iFrequency), C(sl), G_stat.unitinfo(sl,1));
                                    mn(iFrequency) = stt.means_semiweighted;
                                    sm(iFrequency) = stt.sems_semiweighted;
                                end
                                %mn = squeeze(nanmean(dat,1));
                                %sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                                freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).freq;
                                errorbar(freq, mn,sm,cors{iStim})
                                hold on
                                ylim([0 0.05])
                                latencies = {'prestim', 'stim'};
                                stims = {'drift', 'bar'};
                                layers = {'nan', 'nan','superficial', 'granular', 'deep'};
                                nt = {'pyr', 'int'};
                                H = title(sprintf('unit: %s  LFP: %s, %s, %s, %s',  layers{iLayer}, G_stat.layers{jLayer}, latencies{iLatency}, G_stat.state{iState}, nt{k+1}));
                                set(H,'FontSize', 6)
                            end                        
                        end
                    end
                end
        end
    end
end
%% plot the different layer 4 cells
for iType = [1 0 2]
    iMethod = 1; % phsae locking
    iTaper = 1;
    figure
    cnt = 0;
    for iStim = 1:2
     for iState = 1:2
        for iLatency = 1:2
             for jLayer = [2 3 4]
                  cnt = cnt + 1;
                            cors = {'k', 'g'};
                            subplot(8,3,cnt)
                            for iPing= 0:1
                                hold on
                                try                                            
                                    sl =  G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==0 & G_stat.unitinfo(:,2)==iPing & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                    dat = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
                                    dat(dat==0) = NaN;
                                    dat(~isfinite(dat)) = NaN;

                                    mn = squeeze(nanmean(dat,1));
                                    sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                                    freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).freq;
                                    errorbar(freq, mn,sm,cors{iPing+1})
                                    hold on
                                    ylim([0 0.05])
                                    latencies = {'prestim', 'stim'};
                                    stims = {'drift', 'bar'};
                                    layers = {'nan', 'nan','superficial', 'granular', 'deep'};
                                    H = title(sprintf('stim=%s,st=%s,lat=%s \n  Layer LFP: %s', stims{iStim},G_stat.state{iState},latencies{iLatency},G_stat.layers{jLayer}));
                                    set(H,'FontSize',6)
                                end                        
                            end
                end
            end
        end
    end
end
%%

% first figure: relative to different layers, with running and with sitting
iType = 1
iMethod = 1; % phsae locking
iTaper = 1;
for iStim = 1:2
   for iState = 1:2
       for iLatency = 1:2
            figure
            cnt = 0;
            for k = 0:1       
                for iLayer = 3:5
                   for jLayer = [2 3 4]

                        cors = {'r', 'k', 'g'};
                        cnt=cnt+1;
                        subplot(6,3,cnt)
                        for iType = [1 0 2]
                            hold on
                            try                                            
                                sl =  G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==k & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                dat = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
                                dat(dat==0) = NaN;
                                dat(~isfinite(dat)) = NaN;

                                [A,B,C] = unique(G_stat.animalinfo);
                                cfg.semiweighted_adaptive = 'yes';
                                for iFrequency = 1:size(dat,2)                                    
                                    stt = Lur_statistics(cfg, dat(:,iFrequency), C(sl), G_stat.unitinfo(sl,1));
                                    mn(iFrequency) = stt.means_semiweighted;
                                    sm(iFrequency) = stt.sems_semiweighted;
                                end
%                                mn = squeeze(nanmean(dat,1));
 %                               sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                                freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).freq;
                                errorbar(freq, mn,sm,cors{iType+1})
                                hold on
                                ylim([0 0.05])
                                latencies = {'prestim', 'stim'};
                                stims = {'drift', 'bar'};
                                layers = {'nan', 'nan','superficial', 'granular', 'deep'};
                                H = title(sprintf('Layer unit: %s   Layer LFP: %s',  layers{iLayer}, G_stat.layers{jLayer}));
                                if cnt==1
                                    H = title(sprintf('%s, %s, %s \n Layer unit: %s   Layer LFP: %s',  stims{iStim}, G_stat.state{iState}, latencies{iLatency}, layers{iLayer}, G_stat.layers{jLayer}));
                                end
                                xlim([0 140])
                                set(H,'FontSize',6)
                            end                        
                        end
                    end
                end
            end
        end
    end
end
%%




