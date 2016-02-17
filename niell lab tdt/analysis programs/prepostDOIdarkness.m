clear all; close all
[f p] = uigetfile('*.mat','pre data');
load(fullfile(p,f));

[afname, apname] = uigetfile('*.mat','analysis data');
darkpname = apname;
afile = fullfile(apname,afname);
load(afile);

preSpikes = blockSpike;
preRunT = tsampDark;
preRunV = vsmoothDark;

[f p] = uigetfile('*.mat','post data');
load(fullfile(p,f));

psfilename = 'D:\Angie_analysis\analysisPS.ps';
if exist(psfilename,'file')==2;delete(psfilename);end %%% 

postSpikes = blockSpike;
postRunT = tsampDark;
postRunV = vsmoothDark;


dur = 600;
dt = 1;
%histbins = dt/2:dt:dur

for c = 1:length(preSpikes);
    figure
    for cond = 1:2;
        if cond==1
            col = 'b';
            sp = preSpikes{c};
            preISI  = diff(sp(1:end-1));
            postISI = diff(sp(2:end))
            subplot(2,2,4)
            loglog(preISI,postISI,'b.'); xlabel 'pre ISI'; ylabel 'postISI'; axis square
            hold on
        else
            col = 'r';
            sp = postSpikes{c};
            preISI  = diff(sp(1:end-1));
            postISI = diff(sp(2:end))
            subplot(2,2,4)
            loglog(preISI,postISI,'r.');
        end
        R(c,:,cond) = (hist(sp(sp<dur),dt/2:dt:dur));
        cv2(c,:,cond) = mean(2*abs(preISI-postISI)./(preISI+postISI));
    end
    
    subplot(2,2,1)
    plot(squeeze(R(c,:,:))); xlabel 'time(s)'; ylabel 'sp/sec';
    
    subplot(2,2,2)
    bar(squeeze(cv2(c,:,:))); ylim([0 1.5])
    ax = gca;
    ax.XTick = [1 2];
    ax.XTickLabels = {'pre','post'}; 
    text(0,-10,title_text,'FontSize',8);
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end

layernums = unique(layer);
layerFR = zeros(length(layernums),2); %layerFR(layer,pre/post DOI) = mean FR per layer
layercellsFR = zeros(size(R,1),length(layernums),2); %layercellsFR(cell #, layer, pre/post DOI)
figure
for m = 1:2
    for n=1:length(layernums)
        layerFR(n,m) = mean(mean(R(find(layer==layernums(n)),:,m),2));
        cellrates = mean(R(find(layer==layernums(n)),:,m),2);
         darkerr (n,m) = nanstd(layerFR(:,m)/sqrt(sum(layerFR(n,m))));
      for o = 1:length(find(layer==layernums(n)));
             %darkerr (n,m) = nanstd(layerFR(:,m)/sqrt(sum(layerFR(n,m))));
             layercellsFR(o,n,m) = cellrates(o);  
      end
subplot(2,3,1)
%barweb(mean(layerFR(:,:)),mean(darkerr(:,:)));
bar(mean(layerFR(:,:))); title 'all layers'; hold on
ax = gca;
%ax.XTick = [1 2];
%ax.XTickLabels = {'pre','post'}
title ('all layers')
subplot(2,3,n+1)
bar(layerFR(n,:))%, darkerr(n,:))
%errorbar(layerFR(n,:),darkerr(n))
title(sprintf('layer %d',n+1))
ax = gca;
%ax.XTick = [1 2];
%ax.XTickLabels = {'pre','post'}
    end
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end
legend({'pre','post'})

%note:layerFR(1 = layer 2)

% bar
% high CV = more variation in spiking pattern...less bursty
%cv pre vs cv post
figure
subplot(1,2,1)
scatterhist(cv2(:,:,1),cv2(:,:,2),'Color',[0,.3,1]); xlabel 'CV2 Pre DOI'; ylabel 'CV2 Post DOI'; lsline
ylim([0.6 1.6]); axis square;

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');


for c= 1:length(preSpikes);
    cellspikes = R(c,:,:);
    normR(c,:,:) = cellspikes/max(cellspikes(:));
end


rmax = 25;
figure
subplot(2,1,1);
imagesc(R(:,:,1),[0 rmax]); title('pre'); xlabel('secs'); ylabel('cell #')
subplot(2,1,2)
imagesc(R(:,:,2),[0 rmax]); title('post'); xlabel('secs'); ylabel('cell #')

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

rmax = 10;
figure
subplot(2,1,1);
imagesc(normR(:,:,1),[0 rmax]); title('pre'); xlabel('secs'); ylabel('cell #')
subplot(2,1,2)
imagesc(normR(:,:,2),[0 rmax]); title('post'); xlabel('secs'); ylabel('cell #')

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

preCorr = corrcoef(squeeze(normR(:,:,1))');
postCorr = corrcoef(squeeze(normR(:,:,2))');

figure
subplot(1,2,1);
imagesc(preCorr,[-1 1]); colormap jet; title('Pre DOI'); xlabel('cell #'); ylabel('cell #'); axis square;
subplot(1,2,2); 
imagesc(postCorr,[-1 1]); colormap jet; title('Post DOI');xlabel('cell #'); ylabel('cell #'); 
axis square;

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

preRunV = preRunV(preRunT<dur); preRunT= preRunT(preRunT<dur);
postRunV = postRunV(postRunT<dur);postRunT= postRunT(postRunT<dur);

figure
subplot(2,1,1);
plot(preRunT,preRunV); xlim([0 dur]); xlabel('secs'); ylabel('pre speed')
subplot(2,1,2);
plot(postRunT,postRunV);xlim([0 dur]); xlabel('secs'); ylabel('post speed')

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

meanR = squeeze(mean(R,2));
figure
plot(meanR(:,2),meanR(:,1),'o');hold on; plot([0 20],[0 20]);axis equal
xlabel('post rate'); ylabel('pre rate');

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

meanRpre = mean(meanR(:,1))
meanRpost = mean(meanR(:,2))
barx = [meanRpre meanRpost]

figure
subplot(2,1,1)
barweb(meanR(:,1),meanR(:,2))

subplot(2,1,2)
bar(barx); title 'mean FR'; ylabel 'sp/sec'
set(gca,'Xticklabel',{'pre','post'});

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

allR = [squeeze(normR(:,:,1)) squeeze(normR(:,:,2))];
allV = [preRunV postRunV];
allT = [preRunT postRunT+dur];

figure 
imagesc(allR,[0 1]);
axis xy;
hold on
plot(allT/dt,5*allV/max(allV)-5 ,'g'); ylim([-5 size(allR,1)+0.5]);
plot([dur dur],[-5 size(allR,1)],'r'); colorbar;

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

[coeff score latent] = pca(allR');
figure
plot(latent(1:10)/sum(latent))

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

figure
subplot(2,1,1)
imagesc(coeff); axis square
subplot(2,1,2)

%hist of PCA loadings
hist(coeff); axis square
set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');


figure
hold on
plot(score(:,1),score(:,2))
for i=1:length(score);
plot(score(i,1),score(i,2),'.','Markersize',12,'Color',cmapVar(i,1,length(score),jet))
xlabel 'PC1'; ylabel 'PC2';
end

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

figure
for i = 1:3
    subplot(4,1,i)
    plot(score(:,i)); ylabel(sprintf('pca %d',i))
end
subplot(4,1,4)
plot(allT,allV); ylabel('speed')

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

%timeponts
figure
imagesc(corrcoef(allR)); axis square; xlabel('secs'); ylabel('secs'); colorbar; title('correlation')
colormap jet;

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

%corr units
figure
imagesc(corrcoef(allR'),[-1 1]); axis square; xlabel('cell #'); ylabel('cell #'); colorbar; title('correlation')
colormap jet;

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');
    
[f p] = uiputfile('*.pdf','pdf name');
save(fullfile(p,f),'allR', 'preSpikes', 'postSpikes', 'cv2');

ps2pdf('psfile', psfilename, 'pdffile', fullfile(p,f));
delete(psfilename);