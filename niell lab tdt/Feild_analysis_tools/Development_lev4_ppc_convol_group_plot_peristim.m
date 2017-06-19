function [] = NRKO_lev3_ppc_convol_group_peristim

global info
global outputDir
input   = 'lev2_ppc_convol_peristim';
input2  = 'lev0_spike_peristim';
output  = 'lev3_ppc_convol_group_peristim';
filename = fullfile(outputDir,output,output);
mkdir(fullfile(outputDir,output))
load(filename,'G_stat');
G_stat.unitinfo(G_stat.unitinfo(:,4)==2,4) = 3;

output  = 'lev3_tuningstats_group_peristim';
filename = fullfile(outputDir,output,output);
mkdir(fullfile(outputDir,output))
load(filename);
%
rate_stat.unitinfo(rate_stat.unitinfo(:,4)==2,4) = 3;
%%
% iType = 1;
% iMethod = 1; % phsae locking
% iTaper = 1;
% for iStim = 1:2
%     for iLatency = 1:2
%             figure
%             cnt = 0;
%             for iLayer = 3:5
%                 for k = 0:1
%                    for jLayer = [2 3 4]
% 
%                         cors = {'k', 'g'};
%                         cnt=cnt+1;
%                         subplot(6,3,cnt)
%                                                 
%                         hold on
%                         try 
%                            sl =  G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==k & G_stat.unitinfo(:,1)==iType ...
%                            & G_stat.ppcAllCat(iMethod,1,iLatency,iTaper,jLayer).nSpikes(:)>50 & G_stat.ppcAllCat(iMethod,2,iLatency,iTaper,jLayer).nSpikes(:)>50;
% 
%                             dat = G_stat.ppcAllCat(iMethod,1,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
%                             dat(dat==0) = NaN;
%                             dat(~isfinite(dat)) = NaN;
% 
%                             dat2 = G_stat.ppcAllCat(iMethod,2,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
%                             dat2(dat2==0) = NaN;
%                             dat2(~isfinite(dat2)) = NaN;
%                             d  = dat;
%                             freq = G_stat.ppcAllCat(iMethod,1,iTaper,jLayer).freq;
%                                                      
%                             kk = nearest(freq,62);                                  
%                             x = d(:,kk); x= x(:);
%                             y = rate_stat.loco_mod(sl,iLatency); y= y(:);
%                             datsel = ~isnan(x) & ~isnan(y);
%                             [rho(kk),pval] = corr(x(datsel),y(datsel),'Type', 'Spearman');
%                             %keyboard
%                             plot(x,y,'ko')                                
%                             latencies = {'prestim', 'stim'};
%                             stims = {'drift', 'bar'};
%                             layers = {'nan', 'nan','superficial', 'granular', 'deep'};
%                             nt = {'pyr', 'int'};
%                             H = title(sprintf('unit: %s LFP: %s %s %s %s df=%d',  layers{iLayer}, G_stat.layers{jLayer}, latencies{iLatency}, stims{iStim}, nt{k+1},sum(sl)));
%                             set(H,'FontSize', 6);
%                         end                        
%                     end
%                 end
%         end
%     end
% end
% %%
% 
% iType = 1;
% iMethod = 1; % phsae locking
% iTaper = 1;
% for iStim = 1:2
%     for iLatency = 1:2
%             figure
%             cnt = 0;
%             for iLayer = 3:5
%                 for k = 0:1
%                    for jLayer = [2 3 4]
% 
%                         cors = {'k', 'g'};
%                         cnt=cnt+1;
%                         subplot(6,3,cnt)
%                                                 
%                         hold on
%                         try 
%                            sl =  G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==k & G_stat.unitinfo(:,1)==iType ...
%                            & G_stat.ppcAllCat(iMethod,1,iLatency,iTaper,jLayer).nSpikes(:)>50 & G_stat.ppcAllCat(iMethod,2,iLatency,iTaper,jLayer).nSpikes(:)>50;
% 
%                             dat = G_stat.ppcAllCat(iMethod,1,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
%                             dat(dat==0) = NaN;
%                             dat(~isfinite(dat)) = NaN;
% 
%                             dat2 = G_stat.ppcAllCat(iMethod,2,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
%                             dat2(dat2==0) = NaN;
%                             dat2(~isfinite(dat2)) = NaN;
%                             d  = dat;
%                             freq = G_stat.ppcAllCat(iMethod,1,iTaper,jLayer).freq;
%                                                      
%                             kk = nearest(freq,62);                                  
%                             x = d(:,kk); x= x(:);
%                             y = rate_stat.stim_mod(sl); y= y(:);
%                             datsel = ~isnan(x) & ~isnan(y);
%                             [rho(kk),pval] = corr(x(datsel),y(datsel),'Type', 'Spearman');
%                             %keyboard
%                             plot(x,y,'ko')                                
%                             latencies = {'prestim', 'stim'};
%                             stims = {'drift', 'bar'};
%                             layers = {'nan', 'nan','superficial', 'granular', 'deep'};
%                             nt = {'pyr', 'int'};
%                             H = title(sprintf('unit: %s LFP: %s %s %s %s df=%d',  layers{iLayer}, G_stat.layers{jLayer}, latencies{iLatency}, stims{iStim}, nt{k+1},sum(sl)));
%                             set(H,'FontSize', 6);
%                         end                        
%                     end
%                 end
%         end
%     end
%end
%%

%%

iType = 1;%1= adults
iMethod = 1; % phsae locking
iTaper = 1;
%for iStim = 1:2
    for iState = 1:2
            figure
            cnt = 0;
            for iLayer = 3:5
                for k = 0:1
                   for jLayer = [3 4 5]

                        cors = {'k', 'g'};
                        cnt=cnt+1;
                        subplot(6,3,cnt)
                        for iLatency= 1:2
                            hold on
                            try    
                                if k==0
                                    sl =  rate_stat.stim_mod(:)>0.1 & G_stat.unitinfo(:,3)==iLayer & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                elseif k==1
                                    sl =  rate_stat.stim_mod(:)<0 & G_stat.unitinfo(:,3)==iLayer & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                end
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
                                nt = {'stimmodulated', 'suppressed by stim'};
                                H = title(sprintf('unit: %s LFP: %s %s %s %s df=%d',  layers{iLayer}, G_stat.layers{jLayer}, G_stat.state{iState}, nt{k+1},sum(sl)));
                                set(H,'FontSize', 6);
                            end                        
                        end
                    end
                end
        end
    end
%end
%%
iType =1;
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
                                if k==0
                                    sl =  rate_stat.osi(:)>0.3 & G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==0 & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                elseif k==1
                                    sl =  rate_stat.osi(:)<0.1 & G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==0 & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                end
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
                                nt = {'ori_modulated', 'ori_nonmodulated'};
                                H = title(sprintf('unit: %s LFP: %s %s %s %s df=%d',  layers{iLayer}, G_stat.layers{jLayer}, G_stat.state{iState}, stims{iStim}, nt{k+1},sum(sl)));
                                set(H,'FontSize', 6);
                            end                        
                        end
                    end
                end
        end
    end
end
%%

iType = 1;
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
                                if k==0
                                    sl =  (rate_stat.loco_mod(:,2))>0.35 & G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==0 & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                elseif k==1
                                    sl =  (rate_stat.loco_mod(:,2))<0.2 & G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==0 & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                end
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
                                nt = {'loco_modulated', 'loco_nonmodulated'};
                                H = title(sprintf('unit: %s LFP: %s %s %s %s df=%d',  layers{iLayer}, G_stat.layers{jLayer}, G_stat.state{iState}, stims{iStim}, nt{k+1},sum(sl)));
                                set(H,'FontSize', 6);
                            end                        
                        end
                    end
                end
        end
    end
end
%%
iType = 1;
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
%% PRE STIM VS STIM; WILD TYPES
% first figure: relative to different layers, with running and with sitting
iType = 1;
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
% % EARLY VS LATE QUIESCENCE
% iType = 1
% iMethod = 1; % phsae locking
% iTaper = 1;
% for iStim = 1:2
%     for iLatency = 1:2
%             figure
%             cnt = 0;
%             for iLayer = 3:5
%                 for k = 0:1
%                    for jLayer = [2 3 4]
% 
%                         cors = {'g', 'r', 'g', 'r'};
%                         cnt=cnt+1;
%                         subplot(6,3,cnt)
%                         for iState= [3:4]
%                             hold on
%                             try                                            
%                                 sl =  G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==k & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
%                                 dat = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
%                                 dat(dat==0) = NaN;
%                                 dat(~isfinite(dat)) = NaN;
% 
%                                 mn = squeeze(nanmean(dat,1));
%                                 sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
%                                 freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).freq;
%                                 errorbar(freq, mn,sm,cors{iState})
%                                 hold on
%                                 ylim([0 0.05])
%                                 latencies = {'prestim', 'stim'};
%                                 stims = {'drift', 'bar'};
%                                 layers = {'nan', 'nan','superficial', 'granular', 'deep'};
%                                title(sprintf('Layer unit: %s \n  Layer LFP: %s',  layers{iLayer}, G_stat.layers{jLayer}));      nt = {'pyr', 'int'};
%                                 H = title(sprintf('unit: %s LFP: %s %s %s %s df=%d',  layers{iLayer}, G_stat.layers{jLayer}, latencies{iLatency}, stims{iStim}, nt{k+1},sum(sl)));
%                                 set(H,'FontSize', 6);                 
%                             end                        
%                         end
%                     end
%                 end
%         end
%     end
% end
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
                                    H = title(sprintf('unit: %s LFP: %s %s %s %s df=%d',  layers{iLayer}, G_stat.layers{jLayer}, latencies{iLatency}, stims{iStim}, G_stat.state{iState},sum(sl)));
                                set(H,'FontSize', 6);                 
                  
                            end                        
                        end
                    end
                end
        end
    end
end
%% plot the different layer 4 cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
gt = {'N2BKO','wt','N2AKO'};
for iType = [0:2]
    iMethod = 1; % phsae locking
    iTaper = 1;
    figure
    cnt = 0;
    for iStim = 1:2
     for iState = 1:2
        for iLatency = 1:2
             for jLayer = [1 2 3 4]
                  cnt = cnt + 1;
                            cors = {'k', 'g'};
                            subplot(8,4,cnt)
                            df = 0;
                            for iPing= 0:1
                                hold on
                                try               
                                    if iPing==1
                                        sl = rate_stat.osi(:)>0.45 & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==0 & G_stat.unitinfo(:,2)==iPing & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                    elseif iPing==0
                                        sl =  rate_stat.osi(:)>0.55 & G_stat.unitinfo(:,4)==4 & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==0 & G_stat.unitinfo(:,2)==iPing & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                    end                                        
                                    dat = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
                                    dat(dat==0) = NaN;
                                    dat(~isfinite(dat)) = NaN;
                                    df(iPing+1) = sum(sl);
                                    mn = squeeze(nanmean(dat,1));
                                    sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                                    freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).freq;
                                    errorbar(freq, mn,sm,cors{iPing+1})
                                    hold on
                                    ylim([0 0.05])
                                    xlim([0 140])
                                    latencies = {'prestim', 'stim'};
                                    stims = {'drift', 'bar'};
                                    layers = {'nan', 'nan','superficial', 'granular', 'deep'};
                                    H = title(sprintf('%s %s %s %s %s df1=%d,df2=%d',  gt{iType+1},G_stat.layers{jLayer}, latencies{iLatency}, stims{iStim}, G_stat.state{iState},df(1), df(2)));
                                    set(H,'FontSize', 6);                                   
                                end                        
                            end
                end
            end
        end
    end
end
%%

%first figure: relative to different layers, with running and with sitting
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
                        df = [];
                        cnt=cnt+1;
                        subplot(6,3,cnt)
                        for iType = [1 0 2]
                            hold on
                            try                                            
                                sl =  G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,2)==k & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                dat = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
                                dat(dat==0) = NaN;
                                dat(~isfinite(dat)) = NaN;

                                [A,B,C] = unique(G_stat.animalinfo(sl));
                                cfg.semiweighted_adaptive = 'yes';
                                for iFrequency = 1:size(dat,2)                                    
                                    stt = Lur_statistics(cfg, dat(:,iFrequency), C, G_stat.unitinfo(sl,1));
                                    mn(iFrequency) = stt.means_semiweighted;
                                    sm(iFrequency) = stt.sems_semiweighted;
                                end
                                sttAll(iType).stt = stt;
                                df(iType + 1)  = length(A);
                                mn = squeeze(nanmean(dat,1));
                                sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                                freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).freq;
                                errorbar(freq, mn,sm,cors{iType+1})
                                hold on
                                ylim([0 0.05])
                                latencies = {'prestim', 'stim'};
                                stims = {'drift', 'bar'};
                                layers = {'nan', 'nan','superficial', 'granular', 'deep'};
                                H = title(sprintf('Layer unit: %s   Layer LFP: %s',  layers{iLayer}, G_stat.layers{jLayer}));
                                if cnt==1
                                    H = title(sprintf('%s, %s, %s \n unit: %s  LFP: %s, NRB(red)=%d,WT(black)=%d,NRA(green)=%d',  stims{iStim}, G_stat.state{iState}, latencies{iLatency}, layers{iLayer}, G_stat.layers{jLayer},df(1),df(2),df(3)));
                                end
                                xlim([0 140])
                                set(H,'FontSize',6)
                                mnAll(iType+1,:) = mn;
                                smAll(iType+1,:) = sm;
                                
                            end                        
                        end
                    end
                end
            end
        end
    end
end
%%

%% STIM VS PRE STIM all genotypes all cells that are orientation selective
iType = 1;
iMethod = 1; % phsae locking
iTaper = 1;
for iStim = 1:2
    for iState = 1:2
            figure
            cnt = 0;
            for iLayer = 3:5
                for iType = 0:2
                   for jLayer = [2 3 4]

                        cors = {'k', 'g'};
                        cnt=cnt+1;
                        subplot(9,3,cnt)
                        for iLatency= 1:2
                            hold on
                            try                                            
                                sl =  rate_stat.osi(:)>0.55  & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==0 & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
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
                                H = title(sprintf('unit: %s LFP: %s %s %s %s df=%d', layers{iLayer},  G_stat.layers{jLayer}, G_stat.state{iState}, stims{iStim}, gt{iType+1},sum(sl)));
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
for iStim = 1:2
    for iLatency = 1:2
            figure
            cnt = 0;
            for iLayer = 3:5
                for iType = 0:2
                   for jLayer = [2 3 4]

                        cors = {'g', 'r'};
                        cnt=cnt+1;
                        subplot(9,3,cnt)
                        for iState = 1:2
                            hold on
                            try                                            
                                sl =  G_stat.unitinfo(:,4)==iLayer & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,5)==0 & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                dat = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).ppcAll(sl,:,:);
                                dat(dat==0) = NaN;
                                dat(~isfinite(dat)) = NaN;

                                mn = squeeze(nanmean(dat,1));
                                sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                                freq = G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper,jLayer).freq;
                                errorbar(freq, mn,sm,cors{iState})
                                hold on
                                ylim([0 0.07])
                                latencies = {'prestim', 'stim'};
                                stims = {'drift', 'bar'};
                                layers = {'nan', 'nan','superficial', 'granular', 'deep'};
                                nt = {'pyr', 'int'};
                                H = title(sprintf('unit: %s LFP: %s %s %s %s df=%d',  layers{iLayer}, G_stat.layers{jLayer},latencies{iLatency}, stims{iStim}, gt{iType+1},sum(sl)));
                                set(H,'FontSize', 6);
                            end                        
                        end
                    end
                end
        end
    end
end
%%
    iMethod = 1; % phsae locking
    iTaper = 1;
    for iStim = 1:2
        for iState = 1:2
                figure
                cnt = 0;
                   for iPing = 0:1            
                        for iType = [1 0 2]               
                            for jLayer = [2 3 4]

                            cors = {'k', 'g'};
                            cnt=cnt+1;
                            subplot(6,3,cnt)
                            for iLatency= 1:2
                                hold on
                                try                                    
                                    if iPing==0
                                        sl =  G_stat.unitinfo(:,4)==4 & G_stat.unitinfo(:,5)==0 & G_stat.unitinfo(:,6)==iStim & G_stat.unitinfo(:,2)==iPing & G_stat.unitinfo(:,1)==iType & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                    elseif iPing==1
                                        sl =  G_stat.unitinfo(:,6)==iStim  & G_stat.unitinfo(:,1)==iType & G_stat.unitinfo(:,2)==iPing & G_stat.ppcAllCat(iMethod,iState,iLatency,iTaper).nSpikes(:)>50;
                                    end
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
                                    ping = {'noping', 'ping'}
                                    H = title(sprintf('unit: %s LFP: %s %s %s %s %s df=%d',  layers{iLayer}, G_stat.layers{jLayer}, G_stat.state{iState}, stims{iStim}, ping{iPing+1},gt{iType+1},sum(sl)));
                                    set(H,'FontSize', 6);
                                end                        
                            end
                       end
                  end
              end
        end
    end





