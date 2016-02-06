clear all; close all
[f p] = uigetfile('*.mat','pre data');
load(fullfile(p,f));

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

binSpikes=10;
for c = 1:length(preSpikes);
figure
    for cond =1:2;
        if cond==1
            sp = preSpikes{c};
            preISI  = diff(sp(1:end-1));
            postISI = diff(sp(2:end))
            cv21(c) = mean(2*abs(preISI-postISI)./(preISI-postISI));
            %subplot(2,2,2)
            %bar(cv2(c),'b'); legend({'pre','post'}); axis square
           % hold on
            subplot(2,2,4)
            loglog(preISI,postISI,'b.'); xlabel 'pre ISI'; ylabel 'postISI'; axis square
            hold on 
        else
            sp = postSpikes{c};
            preISI  = diff(sp(1:end-1));
            postISI = diff(sp(2:end))
            cv22(c) = mean(2*abs(preISI-postISI)./(preISI-postISI));
            %subplot(2,2,2)
            %bar(cv2(c),'r');%legend({'post'})
            subplot(2,2,4)
            loglog(preISI,postISI,'r.');
        end
        R(c,:,cond) = hist(sp(sp<dur),dt/2:dt:dur);
        %cv2(c) = mean(2*abs(preISI-postISI)./(preISI-postISI));
    end
    subplot(2,2,1)
    plot(squeeze(R(c,:,:))); xlabel 'time(s)'; ylabel 'sp/sec'
    legend({'pre','post'})
    subplot(2,2,2)
    barcv2 = [cv21 cv22]
    bar(barcv2)
    %barcv2 = [meanRpre meanRpost]
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end


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
imagesc(preCorr,[-1 1]); colormap jet; title('pre'); xlabel('cell #'); ylabel('cell #');
subplot(1,2,2); 
imagesc(postCorr,[-1 1]); colormap jet; title('post');xlabel('cell #'); ylabel('cell #'); 

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
bar(barx)
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
plot([dur dur],[-5 size(allR,1)],'r');

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

[coeff score latent] = pca(allR');
figure
plot(latent(1:10)/sum(latent))

set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

figure
imagesc(coeff)
set(gcf, 'PaperPositionMode', 'auto');
print('-dpsc',psfilename,'-append');

figure
hold on
plot(score(:,1),score(:,2))
for i=1:length(score);
plot(score(i,1),score(i,2),'.','Markersize',12,'Color',cmapVar(i,1,length(score),jet))
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
save(fullfile(p,f),'allR', 'preSpikes', 'postSpikes');

ps2pdf('psfile', psfilename, 'pdffile', fullfile(p,f));
delete(psfilename);
