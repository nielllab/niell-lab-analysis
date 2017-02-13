%load('compileAll_withLFP_121916.mat')

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
