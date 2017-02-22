%load('compileAll_withLFP_121916.mat')

clear darkSpikes
n=0;
for i = 1:length(use)
    
    
    afile = [pathname '\' files(use(i)).dir '\' files(use(i)).analysisfile '.mat'];
    clustfile = [pathname '\' files(use(i)).dir '\' files(use(i)).clusterfile '.mat'] ;
    
    
 
    [inh mid] = getWaveform(clustfile,afile,0);
    nc = length(inh); cellrange = n+1:n+nc;
    
    %%%if ~isempty(files(use(i)).blockDark{1}) & ~isempty(files(use(i)).blockDark{2})
 %   try
        %%% get wn response
         clear st_lfp
         for prepost = 1:2
            spikes = getSpikes(clustfile,afile,files(use(i)).blockWn{prepost},0);
            lfpraw = getLFPraw(clustfile,afile,files(use(i)).blockWn{prepost},0);
            lfp = median(lfpraw.data,2);
           
            for j = 1:length(cellrange)
                j
                s = spikes.sp{j};
                darkSpikes{cellrange(j),prepost} = s;
                s = s(s>1.5 & s<max(lfpraw.t)-1.5);
                clear lfpsegs
                for n = 1:length(s);
                 
                  [err ind] = min(abs(lfpraw.t - s(n)));
                    lfpsegs(n,:) = lfp(ind-768:ind+768);
                end
                st_lfp(j,:,prepost) = median(lfpsegs,1);
                ns(j,prepost) = length(s);
            end
 
            
        end
        for j = 1:length(cellrange);
          
            figure
                plot((1:length(lfpsegs))*(lfpraw.t(2)-lfpraw.t(1)),squeeze(st_lfp(j,:,:)));
                title(sprintf('pre =%d post=%d',ns(j,1),ns(j,2)))
        end
        
        figure
       plot((1:length(lfpsegs))*(lfpraw.t(2)-lfpraw.t(1)),squeeze(mean(st_lfp,1)));
        keyboard
   % catch
        for j = 1:length(cellrange)
            darkSpikes{cellrange(j),1} = [];
            darkSpikes{cellrange(j),1} = [];
            
        end
   % end
    n= n+nc;
end


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
        plot(mean(squeeze(meanLfpLayer(4:5,:,ceil(j/2),mod(j-1,2)+1)),1)); ylim([0 4000])
        title(titles{j});
    end
end
