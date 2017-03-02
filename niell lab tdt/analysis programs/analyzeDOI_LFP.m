%load('compileAll_withLFP_121916.mat')
close all
clear darkSpikes
clear stLFPall stLFPallMean stLFPallStd nsAll
ntotal=0;
for i = 1:length(use)
    
    i
    afile = [pathname '\' files(use(i)).dir '\' files(use(i)).analysisfile '.mat'];
    clustfile = [pathname '\' files(use(i)).dir '\' files(use(i)).clusterfile '.mat'] ;
    load(afile,'cells');
    
    
    [inh mid] = getWaveform(clustfile,afile,0);
    nc = length(inh); cellrange = (ntotal+1):(ntotal+nc);
    nc
    cellrange
    
    %%%if ~isempty(files(use(i)).blockDark{1}) & ~isempty(files(use(i)).blockDark{2})
    try
    %%% get wn response
    clear st_lfp ns st_lfpMean st_lfpStd
    for prepost = 1:2
        spikes = getSpikes(clustfile,afile,files(use(i)).blockWn{prepost},0);
        lfpraw = getLFPraw(clustfile,afile,files(use(i)).blockWn{prepost},0);
        spd = getSpeed(clustfile,afile,files(use(i)).blockWn{prepost},0);
        
        lfp = median(lfpraw.data,2);
        lfpsites = lfpraw.data;
        dt =mean(diff(lfpraw.t))
        for j = 1:length(cellrange)
            %   j
            
            spikespeed = interp1(spd.t,spd.v,spikes.sp{j});
            site = ceil(j/4);
            
            for move =1:2
                s = spikes.sp{j};
                darkSpikes{cellrange(j),prepost} = s;
                if move==1
                    s = s(s>1.5 & s<max(lfpraw.t)-1.5 &spikespeed<0.5);
                else
                    s = s(s>1.5 & s<max(lfpraw.t)-1.5 &spikespeed>0.5);
                end
                clear lfpsegs
                for n = 1:length(s);
                    
                    %[err ind] = min(abs(lfpraw.t - s(n)));
                    ind = round(s(n)/dt);
                    lfpsegs(n,:) = lfpsites(ind-768:ind+768,site);
                end
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
    for j = 1:length(cellrange);
        
        %             figure
        %                 plot((1:length(lfpsegs))*(lfpraw.t(2)-lfpraw.t(1)),squeeze(st_lfp(j,:,:)));
        %                 title(sprintf('pre =%d post=%d %s',ns(j,1),ns(j,2),files(use(i)).treatment))
        %                 drawnow
        
        thisLFP = st_lfpMean(j,:,:,:);
        range = [min(thisLFP(:)) max(thisLFP(:))];
        figure
        for prepost=1:2
            for mv = 1:2
                subplot(2,2,2*(prepost-1)+mv)
                plot((1:length(lfpsegs))*(lfpraw.t(2)-lfpraw.t(1)),squeeze(st_lfpMean(j,:,prepost,mv)));
                ylim(range); xlim([0 2])
                title(sprintf('n=%d %s mv%d prepost%d',ns(j,prepost,mv),files(use(i)).treatment,mv,prepost))
                drawnow
            end
        end
    end
    
    stLFPall(cellrange,:,:,:) = st_lfp;
    size(stLFPall)
    
    stLFPallMean(cellrange,:,:,:) = st_lfpMean;
    stLFPallStd(cellrange,:,:,:) = st_lfpStd;
    
    nsAll(cellrange,:,:) = ns;
%     figure
%     plot((1:length(lfpsegs))*(lfpraw.t(2)-lfpraw.t(1)),squeeze(mean(st_lfp,1)));
    
    catch
    sprintf('couldnt do %s', files(use(i)).expt)
    
    end
    ntotal= ntotal+nc;
    close all
end
figure
imagesc(squeeze(stLFPall(:,:,1,2)),[-20 20])

treatment = treatment(1:size(stLFPall,1));

for treat = 1:2
    figure
    meanAll = squeeze(nanmean(stLFPall(treatment==treat,:,:,:),1));

        for prepost=1:2
            for mv = 1:2
                subplot(2,2,2*(prepost-1)+mv)
                plot(meanAll(:,prepost,mv));
              ylim([-10 10])
                drawnow
            end
        end
end

        keyboard
        

titles = {'pre stop','post stop','pre move','post move'};

freqs = lfpMoveDark.freq;
range = [0 10^4];

midlayers = layerSites(:,3:4:end);

use = ones(1,size(LFPallCh,1));

%%% dark
badrecs = [23 27 28 41]; use(badrecs)=0;

%%%wn
badrecs = [ 11 27 28 38]; use(badrecs)=0;


treat = {'saline','doi'};
for i = 1:size(LFPallChDark,1);
    figure
    if sessionTreatment(i)<=2
        set(gcf,'Name',treat{sessionTreatment(i)});
    end
    
    for layer = 0:6;
        lfpLayer(i,layer+1,:,:,:) = squeeze(nanmean(LFPallCh(i,midlayers(i,:)==layer,:,:,:),2));
    end
    for j = 1:4
        subplot(2,2,j);
        imagesc(squeeze(lfpLayer(i,:,:,ceil(j/2),mod(j-1,2)+1)),range);
        title(titles{j});
    end
end


meanLFPclean = meanLfpLayer;
meanLFPclean(:,117:124,:,:)=NaN;

for t = 1:2
    
    meanLfpLayer = squeeze(nanmedian(lfpLayer(use==1 & sessionTreatment==t,:,:,:,:),1));
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
end

figure
plot(freqs(1:140),specdata([3 1],:),'Linewidth',2); ylim([0 6000])
xlabel('freq (hz)'); set(gca,'Ytick',0:1000:6000);set(gca,'YtickLabel',{'0','1','2','3','4','5','6'});
set(gca,'Fontsize',20); xlim([0 70]); legend('run','stop'); ylabel('power a.u.'); title('pre DOI')


figure
plot(freqs(1:140),specdata([4 2],:),'Linewidth',2); ylim([0 6000])
xlabel('freq (hz)'); set(gca,'Ytick',0:1000:6000);set(gca,'YtickLabel',{'0','1','2','3','4','5','6'});
set(gca,'Fontsize',20); xlim([0 70]); legend('run','stop'); ylabel('power a.u.'); title('post DOI')

figure
plot(freqs(1:140),specdata([3 4],:),'Linewidth',2); ylim([0 6000])
xlabel('freq (hz)'); set(gca,'Ytick',0:1000:6000);set(gca,'YtickLabel',{'0','1','2','3','4','5','6'});
set(gca,'Fontsize',20); xlim([0 70]); legend('pre','post'); ylabel('power a.u.'); title('running only')

figure
plot(freqs(1:140),specdata([1 2],:),'Linewidth',2); ylim([0 6000])
xlabel('freq (hz)'); set(gca,'Ytick',0:1000:6000);set(gca,'YtickLabel',{'0','1','2','3','4','5','6'});
set(gca,'Fontsize',20); xlim([0 70]); legend('pre','post'); ylabel('power a.u.'); title('stop only')


