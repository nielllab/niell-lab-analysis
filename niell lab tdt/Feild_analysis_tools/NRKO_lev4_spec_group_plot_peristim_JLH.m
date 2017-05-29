function [] = NRKO_lev3_ppc_convol_group_peristim

global info
global outputDir
input   = 'lev2_ppc_convol_peristim';
input2  = 'lev0_spike_peristim';
output  = 'lev2_lfp_spectrum_group_peristim';
filename = fullfile(outputDir,output,output);
mkdir(fullfile(outputDir,output))
load(filename,'G_stat');
%%
for iDir = 1:length(info)
     try
        if strcmp(info(iDir).stimtype,'drift')
            unitinfo(iDir) = 1;  
        else
            unitinfo(iDir) = 2;
        end
    catch
        if strcmp(info(iDir).stimType,'drift')
            unitinfo(iDir) = 1;  
        else
            unitinfo(iDir)= 2;
        end
    end
    
end
%% STIM VS preSTIM; WILD TYPES
% first figure: relative to different layers, with running and with sitting
iStim = 2;
iMethod = 1; 
iTaper = 1;
 figure
cnt = 0;
for iType = 1:3         
    for iState = 1:2
         for iLayer = 1:3

            cors = {'g', 'k'};
            cnt=cnt+1;
            subplot(6,3,cnt)
            for iLatency= 1:2
                hold on
                try                                            
                    sl =  unitinfo==iStim & G_stat.dirinfo==iType-1;%1=wt 0=N2B, 2=N2A
                    dat = log10(G_stat.powAllCat(iState,iLatency,iLayer).powAll(sl,:,:));
                    dat(dat==0) = NaN;
                    dat(~isfinite(dat)) = NaN;

                    mn = squeeze(nanmean(dat,1));
                    sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                    freq = 4:4:160;
                    errorbar(log10(freq), mn,sm,cors{iLatency})
                    hold on
                    latencies = {'stim', 'prestim'};
                    Geno={'2B', 'wt', '2A'};
                    stims = {'drift', 'bar'};
                    layers = {'superficial', 'granular', 'deep'};
                    %nt = {'pyr', 'int'};
                    xlim(log10([4 160]))
                    ylim([0 2])
                    H = title(sprintf('%s %s %s %s',  layers{iLayer}, G_stat.state{iState}, stims{iStim},Geno{iType}, sum(sl)));
                    set(H,'FontSize', 6);
                  
                end         
            end
        end
    end
end

%%
iStim = 2; %iStim of 2 is bar stimulus
iMethod = 1; 
iTaper = 1;
 figure
cnt = 0;
iLatency=2;%iLatency of 1 is with stimulus
%for iType = 1:3         
    for iState = 1:2
         for iLayer = 1:2

            cors = {'g', 'k','r'};%2b, wt, 2a
            cnt=cnt+1;
            subplot(2,2,cnt)
            for iType= 1:3
                hold on
                try                                            
%                     sl =  unitinfo==iStim & G_stat.dirinfo==iType-1;
%                     dat = log10(G_stat.powAllCat(iState,iLatency,iLayer).powAll(sl,:,:));
%                     dat(dat==0) = NaN;
%                     dat(~isfinite(dat)) = NaN;
% 
%                     mn = squeeze(nanmean(dat,1));
%                     sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
%                     freq = 4:4:160;
%                     shadedErrorBar(log10(freq), mn,sm,cors{iType})
%                     hold on
%                     %latencies = {'stim', 'prestim'};
%                     Geno={'2B', 'wt', '2A'};
%                     stims = {'drift', 'bar'};
%                     layers = {'superficial', 'granular', 'deep'};
%                     %nt = {'pyr', 'int'};
%                     xlim(log10([4 160]))
%                     ylim([0 2])
%                     H = title(sprintf('%s %s %s',  layers{iLayer}, G_stat.state{iState}, stims{iStim}, sum(sl)));
%                     set(H,'FontSize', 6);
%                     
                    %for a linear scale
                     sl =  unitinfo==iStim & G_stat.dirinfo==iType-1;
                    dat = G_stat.powAllCat(iState,iLatency,iLayer).powAll(sl,:,:);
                    dat(dat==0) = NaN;
                    dat(~isfinite(dat)) = NaN;

                    mn = squeeze(nanmean(dat,1));
                    sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                    mLowG(iType,1)=nanmedian(mn(10:22));
                    stdLowG(iType,2)=semedian(mn(10:22));
                    freq = 8:4:120;
                    shadedErrorBar((freq), mn(2:30),sm(2:30),cors{iType})
                    hold on
                    %latencies = {'stim', 'prestim'};
                    Geno={'2B', 'wt', '2A'};
                    stims = {'drift', 'bar'};
                    layers = {'superficial', 'granular', 'deep'};
                    %nt = {'pyr', 'int'};
                    xlim([15 75])
                    ylim([0 30])
                    H = title(sprintf('%s %s %s',  layers{iLayer}, G_stat.state{iState}, stims{iStim}, sum(sl)));
                    set(H,'FontSize', 6);
%                  keyboard
                end         
            end
        end
    end
%end

%% plot median power at low gama band (20Hz-44Hz) and high gamma (46-65Hz)
iStim = 2; %iStim of 2 is bar stimulus
iMethod = 1; 
iTaper = 1;
cnt = 0;
iLatency=2;%iLatency of 1 is with stimulus
iState=1;%iState=1 is moving
iLayer=1;%iLayer=1 is superficial iLayer=2 is granular and iLayer=3 is subgranular
%for iType = 1:3         
   % for iState = 1:2
    %     for iLayer = 1:2
            figure
            cors = {'g', 'k','r'};%2b, wt, 2a
            cnt=cnt+1;
            %subplot(2,2,cnt)
            for iType= 1:3
                hold on
                try                                            
%                     sl =  unitinfo==iStim & G_stat.dirinfo==iType-1;
%                     dat = log10(G_stat.powAllCat(iState,iLatency,iLayer).powAll(sl,:,:));
%                     dat(dat==0) = NaN;
%                     dat(~isfinite(dat)) = NaN;
% 
%                     mn = squeeze(nanmean(dat,1));
%                     sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
%                     freq = 4:4:160;
%                     shadedErrorBar(log10(freq), mn,sm,cors{iType})
%                     hold on
%                     %latencies = {'stim', 'prestim'};
%                     Geno={'2B', 'wt', '2A'};
%                     stims = {'drift', 'bar'};
%                     layers = {'superficial', 'granular', 'deep'};
%                     %nt = {'pyr', 'int'};
%                     xlim(log10([4 160]))
%                     ylim([0 2])
%                     H = title(sprintf('%s %s %s',  layers{iLayer}, G_stat.state{iState}, stims{iStim}, sum(sl)));
%                     set(H,'FontSize', 6);
%                     
                    %for a linear scale
                     sl =  unitinfo==iStim & G_stat.dirinfo==iType-1;
                    dat = G_stat.powAllCat(iState,iLatency,iLayer).powAll(sl,:,:);
                    dat(dat==0) = NaN;
                    dat(~isfinite(dat)) = NaN;

                    mn = squeeze(nanmean(dat,1));
                    sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                   
                    mLowG(iType,1)=nanmedian(mn(6:13));
                    mLowG(iType,2)=semedian(mn(6:13));
                    
                    mHighG(iType,1)=nanmedian(mn(16:23));
                    mHighG(iType,2)=semedian(mn(16:23));
                    
                    freq = 8:4:120;
                    %figure
                    shadedErrorBar((freq), mn(2:30),sm(2:30),cors{iType})
                    hold on
                    %latencies = {'stim', 'prestim'};
                    Geno={'2B', 'wt', '2A'};
                    stims = {'drift', 'bar'};
                    layers = {'superficial', 'granular', 'deep'};
                    %nt = {'pyr', 'int'};
%                     xlim([15 75])
%                     ylim([0 30])
                    H = title(sprintf('%s %s %s',  layers{iLayer}, G_stat.state{iState}, stims{iStim}, sum(sl)));
                    set(H,'FontSize', 6);
                 %keyboard
                end         
            end
       % end
    %end
%end

% rearrange matrix since GTs are out of order, N2B is iType=1, WT=iType==2
% and N2A is iType=3
low(1,:)=mLowG(2,:);
low(2,:)=mLowG(3,:);
low(3,:)=mLowG(1,:);

high(1,:)=mHighG(2,:);
high(2,:)=mHighG(3,:);
high(3,:)=mHighG(1,:);

figure
barweb(low(:,1),low(:,2));

figure
barweb(high(:,1),high(:,2));
%%

iType = 2;
iMethod = 1; % phsae locking
iTaper = 1;
 figure
cnt = 0;
for iStim = 1:2         
    for iState = 1:2
         for iLayer = 1:3

            cors = {'g', 'k'};
            cnt=cnt+1;
            subplot(4,3,cnt)
            for iLatency= [2 1]
                hold on
                try                                            
                    sl =  unitinfo==iStim & G_stat.dirinfo==iType;
                    dat = G_stat.powAllCat(iState,iLatency,iLayer).powAll(sl,:,:);
                    dat(dat==0) = NaN;
                    dat(~isfinite(dat)) = NaN;
                    freq = 4:4:160;                                        
                    if iLatency==2
                        norm = repmat(nanmean(dat,2),[1 length(freq)]);
                    end
                    dat = log10(dat./norm);
                    
                    mn = squeeze(nanmean(dat,1));
                    sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                    errorbar(log10(freq), mn,sm,cors{iLatency})
                    hold on
                    latencies = {'stim', 'prestim'};
                    stims = {'drift', 'bar'};
                    layers = {'superficial', 'granular', 'deep'};
                    nt = {'pyr', 'int'};
                    xlim(log10([4 160]))
                    H = title(sprintf('%s %s %s',  layers{iLayer}, G_stat.state{iState}, stims{iStim}, sum(sl)));
                    set(H,'FontSize', 6);
                end         
            end
        end
    end
end
%%
iType = 1;
iMethod = 1; % phsae locking
iTaper = 1;
 figure
cnt = 0;
for iStim = 1:2         
    for iLatency = 1:2
         for iLayer = 1:3

            cors = {'g', 'r', 'c', 'm'};
            cnt=cnt+1;
            subplot(4,3,cnt)
            for iState= [2:-1:1]
                hold on
                try                                            
                    sl =  unitinfo==iStim & G_stat.dirinfo==iType;
                    dat = G_stat.powAllCat(iState,iLatency,iLayer).powAll(sl,:,:);
                    dat(dat==0) = NaN;
                    dat(~isfinite(dat)) = NaN;
                    freq = 4:4:160;                                        
                    if iState==2
                        norm = repmat(nanmean(dat,2),[1 length(freq)]);
                    end
                    dat = log10(dat./norm);
                    
                    mn = squeeze(nanmean(dat,1));
                    sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                    errorbar(log10(freq), mn,sm,cors{iState})
                    hold on
                    latencies = {'stim', 'prestim'};
                    stims = {'drift', 'bar'};
                    layers = {'superficial', 'granular', 'deep'};
                    nt = {'pyr', 'int'};
                    xlim(log10([4 160]))
                    H = title(sprintf('%s %s %s',  layers{iLayer}, latencies{iLatency}, stims{iStim}, sum(sl)));
                    set(H,'FontSize', 6);
                end         
            end
        end
    end
end
%% GT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STIM VS PRESTIM; genotypes
gt = {'N2B', 'WT', 'N2A'}
for iType = 0:2
    iMethod = 1; % phsae locking
    iTaper = 1;
     figure
    cnt = 0;
    for iStim = 1:2         
        for iState = 1:2
             for iLayer = 1:3

                cors = {'g', 'k'};
                cnt=cnt+1;
                subplot(4,3,cnt)
                for iLatency= [2 1]
                    hold on
                    try                                            
                        sl =  unitinfo==iStim & G_stat.dirinfo==iType;
                        dat = G_stat.powAllCat(iState,iLatency,iLayer).powAll(sl,:,:);
                        dat(dat==0) = NaN;
                        dat(~isfinite(dat)) = NaN;
                        freq = 4:4:160;                                        
                        if iLatency==2
                            norm = repmat(nanmean(dat,2),[1 length(freq)]);
                        end
                        dat = log10(dat./norm);

                        mn = squeeze(nanmean(dat,1));
                        sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                        errorbar(log10(freq), mn,sm,cors{iLatency})
                        hold on
                        latencies = {'prestim', 'stim'};
                        stims = {'drift', 'bar'};
                        layers = {'superficial', 'granular', 'deep'};
                        nt = {'pyr', 'int'};
                        xlim(log10([4 160]))
                        H = title(sprintf('%s %s %s %s',  layers{iLayer}, G_stat.state{iState}, stims{iStim}, gt{iType+1}, sum(sl)));
                        set(H,'FontSize', 6);
                    end         
                end
            end
        end
    end
end
%%
gt = {'N2B', 'WT', 'N2A'}
for iType = 0:2
    iMethod = 1; % phsae locking
    iTaper = 1;
     figure
    cnt = 0;
    for iStim = 1:2         
        for iLatency = 1:2
             for iLayer = 1:3

                cors = {'g', 'r'};
                cnt=cnt+1;
                subplot(4,3,cnt)
                for iState= [2 1]
                    hold on
                    try                                            
                        sl =  unitinfo==iStim & G_stat.dirinfo==iType;
                        dat = G_stat.powAllCat(iState,iLatency,iLayer).powAll(sl,:,:);
                        dat(dat==0) = NaN;
                        dat(~isfinite(dat)) = NaN;
                        freq = 4:4:160;                                        
                        if iState==2
                            norm = repmat(nanmean(dat,2),[1 length(freq)]);
                        end
                        dat = log10(dat./norm);

                        mn = squeeze(nanmean(dat,1));
                        sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                        errorbar(log10(freq), mn,sm,cors{iState})
                        hold on
                        latencies = {'prestim', 'stim'};
                        stims = {'drift', 'bar'};
                        layers = {'superficial', 'granular', 'deep'};
                        nt = {'pyr', 'int'};
                        xlim(log10([4 160]))
                        H = title(sprintf('%s %s %s %s',  layers{iLayer}, latencies{iLatency}, stims{iStim}, gt{iType+1}, sum(sl)));
                        set(H,'FontSize', 6);
                    end         
                end
            end
        end
    end
end

%%

gt = {'N2B', 'WT', 'N2A'}
iLayer = 2;
iMethod = 1; % phsae locking
iTaper = 1;
 figure
cnt = 0;
for iStim = 1:2         
    for iLatency = 1:2
         for iType = 0:2

            cors = {'g', 'r'};
            cnt=cnt+1;
            subplot(4,3,cnt)
            for iState= [2:-1:1]
                hold on
                try                                            
                    sl =  unitinfo==iStim & G_stat.dirinfo==iType;
                    dat = G_stat.powAllCat(iState,iLatency,iLayer).powAll(sl,:,:);
                    dat(dat==0) = NaN;
                    dat(~isfinite(dat)) = NaN;
                    freq = 4:4:160;                                        
                    if iState==2
                        norm = repmat(nanmean(dat,2),[1 length(freq)]);
                    end
                    dat = log10(dat./norm);
                    
                    mn = squeeze(nanmean(dat,1));
                    sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                    errorbar(log10(freq), mn,sm,cors{iState})
                    hold on
                    latencies = {'stim', 'prestim'};
                    stims = {'drift', 'bar'};
                    layers = {'superficial', 'granular', 'deep'};
                    nt = {'pyr', 'int'};
                    Genotype={'N2B','wt','N2A'};
                    xlim(log10([4 160]))
                    H = title(sprintf('%s %s %s',  Genotype{iType+1}, latencies{iLatency}, stims{iStim}, sum(sl)));
                    set(H,'FontSize', 6);
                end         
            end
        end
    end
end 












