[f p] = uigetfile('*.mat','pre data');
load(fullfile(p,f));

preSpikes = blockSpike;
preRunT = tsampDark;
preRunV = vsmoothDark;

[f p] = uigetfile('*.mat','post data');
load(fullfile(p,f));

postSpikes = blockSpike;
postRunT = tsampDark;
postRunV = vsmoothDark;

dur = 600;
dt = 1;
for c = 1:length(preSpikes);
    for cond =1:2;
        if cond==1
            sp = preSpikes{c};
        else
            sp = postSpikes{c};
        end
        R(c,:,cond) = hist(sp(sp<dur),dt/2:dt:dur);
    end
    figure
    plot(squeeze(R(c,:,:)));
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

rmax = 10;
figure
subplot(2,1,1);
imagesc(normR(:,:,1),[0 rmax]); title('pre'); xlabel('secs'); ylabel('cell #')
subplot(2,1,2)
imagesc(normR(:,:,2),[0 rmax]); title('post'); xlabel('secs'); ylabel('cell #')

preRunV = preRunV(preRunT<dur); preRunT= preRunT(preRunT<dur);
postRunV = postRunV(postRunT<dur);postRunT= postRunT(postRunT<dur);

figure
subplot(2,1,1);
plot(preRunT,preRunV); xlim([0 dur]); xlabel('secs'); ylabel('pre speed')
subplot(2,1,2);
plot(postRunT,postRunV);xlim([0 dur]); xlabel('secs'); ylabel('post speed')

meanR = squeeze(mean(R,2));
figure
plot(meanR(:,1),meanR(:,2),'o');hold on; plot([0 20],[0 20]);axis equal
xlabel('pre rate'); ylabel('post rate');

allR = [squeeze(normR(:,:,1)) squeeze(normR(:,:,2))];
allV = [preRunV postRunV];
allT = [preRunT postRunT+dur];

figure 
imagesc(allR,[0 1]);
axis xy;
hold on
plot(allT/dt,5*allV/max(allV)-5 ,'g'); ylim([-5 size(allR,1)+0.5]);
plot([dur dur],[-5 size(allR,1)],'r');

[coeff score latent] = pca(allR');
figure
plot(latent(1:10)/sum(latent))

figure
imagesc(coeff)

figure
hold on
plot(score(:,1),score(:,2))
for i=1:length(score);
plot(score(i,1),score(i,2),'.','Markersize',12,'Color',cmapVar(i,1,length(score),jet))
end

figure
for i = 1:3
    subplot(4,1,i)
    plot(score(:,i)); ylabel(sprintf('pca %d',i))
end
subplot(4,1,4)
plot(allT,allV); ylabel('speed')


figure
imagesc(corrcoef(allR)); axis square; xlabel('secs'); ylabel('secs'); colorbar; title('correlation')
colormap jet;

figure
imagesc(corrcoef(allR'),[-1 1]); axis square; xlabel('cell #'); ylabel('cell #'); colorbar; title('correlation')
colormap jet;