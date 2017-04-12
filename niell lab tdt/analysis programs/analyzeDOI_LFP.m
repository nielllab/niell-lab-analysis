%load('compileAll_withLFP_121916.mat')

close all
clear all

batchDOIephys_filtered; %%% load batch file

%%% select the sessions you want based on filters
use =  find(strcmp({files.notes},'good data'))%useSess = use;
sprintf('%d selected sessions',length(use))

saline=1; doi=2; ht=3; ketanserin=4; ketandoi=5; mglur2=6; mglur2doi=7; lisuride=8;
treatLabel = {'saline','DOI','5HT'};
close all
clear darkSpikes
clear stLFPall stLFPallMean stLFPallStd nsAll
ntotal=0;

prepostLabel= {'pre','post'};
%%% calculate spike-triggered LFP

neighborsite = [2 3 4 5 4 5 6 7 10 11 12 13 12 13 14 15];

for i = 1:length(use)-1
%for i = 1:3
    
    i
    afile = [pathname '\' files(use(i)).dir '\' files(use(i)).analysisfile '.mat'];
    clustfile = [pathname '\' files(use(i)).dir '\' files(use(i)).clusterfile '.mat'] ;
    load(afile,'cells');
    
    [inh mid] = getWaveform(clustfile,afile,0);
    nc = length(inh); cellrange = (ntotal+1):(ntotal+nc);
    nc
    cellrange;
    inhAll(cellrange)=inh;
    
    if strcmp(files(use(i)).treatment,'Saline'), treatment(cellrange)=saline; end;
    if strcmp(files(use(i)).treatment,'DOI'), treatment(cellrange)=doi; end;
    if strcmp(files(use(i)).treatment,'5HT'), treatment(cellrange)=ht; end;
    
    
    
    if ~isempty(files(use(i)).blockDrift{1}) & ~isempty(files(use(i)).blockDrift{2})
        %  try
        %%% get wn response
        clear st_lfp ns st_lfpMean st_lfpStd
        
        %%% loop over pre and post
        for prepost = 1:2
            
            %%% load in spikes, LFP, and running speed
            spikes = getSpikes(clustfile,afile,files(use(i)).blockDrift{prepost},0);
            lfpraw = getLFPraw(clustfile,afile,files(use(i)).blockDrift{prepost},0);
            spd = getSpeed(clustfile,afile,files(use(i)).blockDrift{prepost},0);
            
            lfp = median(lfpraw.data,2);  %%% take median across all sites (low noise, but maybe lose local effects?)
            lfpsitesRaw = double(lfpraw.data);       %%% data for each "tetrode"
            dt =mean(diff(lfpraw.t))
            if ~isempty(dt)
                Fs = 1/dt;
            end
            fpass = [1 15];
            [b,a] = butter(4,fpass*2/Fs);
            clear lfpsites
            for ch = 1:16
                lfpsites(:,ch) = filtfilt(b,a,lfpsitesRaw(:,ch));
            end
            
            %%% loop over all cells
            for j = 1:length(cellrange)
                %   j
                
                %%% interpolate speed at each of the spike times
                spikespeed = interp1(spd.t,spd.v,spikes.sp{j});
                
                %%% find site for this cell
                site = ceil(cells(j,1)/4);
                
                %%% calculate spike-triggered LFP for stationary and moving
                for move =1:2
                    s = spikes.sp{j};
                    darkSpikes{cellrange(j),prepost} = s;
                    
                    %%% collect all spikes (with buffer at beginning/end for window, at appropriate speed
                    if move==1
                        s = s(s>1.5 & s<max(lfpraw.t)-1.5 &spikespeed<0.5);
                    else
                        s = s(s>1.5 & s<max(lfpraw.t)-1.5 &spikespeed>1);
                    end
                    
                    clear lfpsegs
                    %%% for each spike, grab the lfp segment around it
                    for n = 1:length(s);
                        ind = round(s(n)/dt);
                        lfpsegs(n,:) = lfpsites(ind-768:ind+768,neighborsite(site));
                    end
                    
                    
                    %%% average lfp segments
                    if length(s)>0
                        st_lfp(j,:,prepost,move) = median(lfpsegs,1);
                        st_lfpMean(j,:,prepost,move) = mean(lfpsegs,1);
                        st_lfpStd(j,:,prepost,move) = std(lfpsegs,[],1);
                        ns(j,prepost,move) = length(s);
                    else
                        st_lfp(j,:,prepost,move) = NaN;
                        st_lfpMean(j,:,prepost,move) = NaN;
                        st_lfpStd(j,:,prepost,move) = NaN;
                        ns(j,prepost,move) = NaN;
                    end
                end
            end
        end
        
        %%% plot figures for each cell (pre/post, move/stationary)
        for j = 1:length(cellrange);
            
            thisLFP = st_lfpMean(j,:,:,:);
            range = [min(thisLFP(:)) max(thisLFP(:))];
            figure
            for prepost=1:2
                for mv = 1:2
                    subplot(2,2,2*(prepost-1)+mv)
                    plot((1:size(lfpsegs,2))*(lfpraw.t(2)-lfpraw.t(1)),squeeze(st_lfpMean(j,:,prepost,mv)));
                    ylim(range); xlim([0 2])
                    title(sprintf('n=%d %s mv%d %s inh%d',ns(j,prepost,mv),files(use(i)).treatment,mv,prepostLabel{prepost},inh(j)))
                    drawnow
                end
            end
        end
        
        stLFPall(cellrange,:,:,:) = st_lfp;
        size(stLFPall)
        
        stLFPallMean(cellrange,:,:,:) = st_lfpMean;
        stLFPallStd(cellrange,:,:,:) = st_lfpStd;
        
        nsAll(cellrange,:,:) = ns;
    end
    %     figure
    %     plot((1:length(lfpsegs))*(lfpraw.t(2)-lfpraw.t(1)),squeeze(mean(st_lfp,1)));
    
    %     catch
    %     sprintf('couldnt do %s', files(use(i)).expt)
    %
    %     end
    ntotal= ntotal+nc;
    close all
end
figure
imagesc(squeeze(stLFPall(:,:,1,2)),[-20 20])

%%% average over treatments, and plot
dt =mean(diff(lfpraw.t))
treatment = treatment(1:size(stLFPall,1));
nsmin = min(min(nsAll,[],3),[],2);
for treat = 1:2
    figure
    meanAll = squeeze(nanmean(stLFPall(treatment==treat & nsmin'>50 & inhAll(1:length(nsmin))>=0,:,:,:),1));
    
    for prepost=1:2
        for mv = 1:2
            subplot(2,2,0*(prepost-1)+mv)
            plot((1:length(lfpsegs))*(2/length(lfpsegs)),meanAll(:,prepost,mv));
            hold on
            ylim([-25 40]); xlim([0.5 1.5])
            drawnow
            title(sprintf('prepost %d mv %d',prepost,mv))
        end
    end
    set(gcf,'Name',(treatLabel{treat}));
end


%%% read in drifting grating spectra, if they aren't already
clear LFPallChDrift
if ~exist('LFPallChDrift','var')
    for i = 1:length(use)
        i
        afile = [pathname '\' files(use(i)).dir '\' files(use(i)).analysisfile '.mat'];
        clustfile = [pathname '\' files(use(i)).dir '\' files(use(i)).clusterfile '.mat'] ;
        
        clear layer
        load(afile,'layer');
        if exist('layer','var')
            layerAll(cellrange) = layer;
        else
            layerInfo = getLayer(clustfile,afile,files(use(i)).tip1,files(use(i)).tip2,files(use(i)).angle, 0); %(needs histo information, but will give layers for all sites)
            layerSites(i,:) = layerInfo.sites;
        end
        
        if strcmp(files(use(i)).treatment,'Saline'), sessionTreatment(i)=saline; end;
        if strcmp(files(use(i)).treatment,'DOI'), sessionTreatment(i)=doi; end;
        if strcmp(files(use(i)).treatment,'5HT'),sessionTreatment(i)=ht; end;
        
        
        if ~isempty(files(use(i)).blockDrift{1}) & ~isempty(files(use(i)).blockDrift{2})
            i
            for prepost=1:2
                lfpMoveDrift = getLfpMovement(clustfile,afile,files(use(i)).blockDrift{prepost},0);
                
                if size(lfpMoveDrift.meanSpect,1)==16
                    LFPallDrift(i,:,:,prepost) =squeeze(nanmedian(lfpMoveDrift.meanSpect, 1));
                    LFPallChDrift(i,:,:,:,prepost) = lfpMoveDrift.meanSpect;
                    LFPfreqDrift(i,:) = lfpMoveDrift.freq;
                    display('good lfp')
                    lfpMoveDrift
                else
                    display('lfp wrong size')
                end
            end
        end
    end
end




%%% average LFP spectra across layers, for each treatment



titles = {'pre stop','post stop','pre move','post move'};

freqs = median(LFPfreqDrift,1);
range = [0 10^4];
midlayers = layerSites(:,3:4:end);


% %%% dark
% badrecs = [23 27 28 41]; use(badrecs)=0;
% %%%wn
% badrecs = [ 11 27 28 38]; use(badrecs)=0;


treat = {'saline','doi'};
for i = 1:size(LFPallChDrift,1);
    figure
    if sessionTreatment(i)<=2
        set(gcf,'Name',treat{sessionTreatment(i)});
    end
    
    for layer = 0:6;
        lfpLayer(i,layer+1,:,:,:) = squeeze(nanmean(LFPallChDrift(i,midlayers(i,:)==layer,:,:,:),2));
    end
    for j = 1:4
        subplot(2,2,j);
        imagesc(squeeze(lfpLayer(i,:,:,ceil(j/2),mod(j-1,2)+1)),range);
        title(titles{j});
    end
end


used = 1;
for t = 1:2
    
    meanLfpLayer = squeeze(nanmedian(lfpLayer(used==1 & sessionTreatment==t,:,:,:,:),1));
    figure; set(gcf,'Name',treat{t});
    for j = 1:4
        subplot(2,2,j);
        imagesc(squeeze(meanLfpLayer(:,:,ceil(j/2),mod(j-1,2)+1)),range/2);
        title(titles{j});
    end
    
    figure; set(gcf,'Name',treat{t});
    for j = 1:4
        subplot(2,2,j);
        d = mean(squeeze(meanLfpLayer(4:5,:,ceil(j/2),mod(j-1,2)+1)),1);
        specdata(j,:) = interp1([1:116 124:140],d([1:116 124:140]),1:140);
        plot(d); ylim([0 5000])
        title(titles{j});
    end
    
    figure
    plot(freqs(1:140),specdata([3 1],:),'Linewidth',2); ylim([0 6000])
    xlabel('freq (hz)'); set(gca,'Ytick',0:1000:6000);set(gca,'YtickLabel',{'0','1','2','3','4','5','6'});
    set(gca,'Fontsize',20); xlim([0 70]); legend('run','stop'); ylabel('power a.u.'); title(['pre ' treatLabel{t}])
    
    
    figure
    plot(freqs(1:140),specdata([4 2],:),'Linewidth',2); ylim([0 6000])
    xlabel('freq (hz)'); set(gca,'Ytick',0:1000:6000);set(gca,'YtickLabel',{'0','1','2','3','4','5','6'});
    set(gca,'Fontsize',20); xlim([0 70]); legend('run','stop'); ylabel('power a.u.');  title(['post ' treatLabel{t}])
    
    figure
    plot(freqs(1:140),specdata([3 4],:),'Linewidth',2); ylim([0 6000])
    xlabel('freq (hz)'); set(gca,'Ytick',0:1000:6000);set(gca,'YtickLabel',{'0','1','2','3','4','5','6'});
    set(gca,'Fontsize',20); xlim([0 70]); legend('pre',['post ' treatLabel{t}]); ylabel('power a.u.'); title('running only')
    
    figure
    plot(freqs(1:140),specdata([1 2],:),'Linewidth',2); ylim([0 6000])
    xlabel('freq (hz)'); set(gca,'Ytick',0:1000:6000);set(gca,'YtickLabel',{'0','1','2','3','4','5','6'});
    set(gca,'Fontsize',20); xlim([0 70]); legend('pre',['post ' treatLabel{t}]); ylabel('power a.u.'); title('stop only')
    
    
    
end



