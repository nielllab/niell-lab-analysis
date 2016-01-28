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
iType = 0;
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
            for iLatency= 1:2
                hold on
                try                                            
                    sl =  unitinfo==iStim & G_stat.dirinfo==iType;
                    dat = log10(G_stat.powAllCat(iState,iLatency,iLayer).powAll(sl,:,:));
                    dat(dat==0) = NaN;
                    dat(~isfinite(dat)) = NaN;

                    mn = squeeze(nanmean(dat,1));
                    sm = squeeze(nanstd(dat,1))./sqrt(sum(sl));
                    freq = 4:4:160;
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












