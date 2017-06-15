function [preCorr postCorr cv2 R eigs] = prepostDOIdarkness(clustfile,afile, blocks, dt, savePDF);
blocknm = blocks{1};
spikes = getSpikes(clustfile,afile, blocknm,0);
spd = getSpeed(clustfile,afile,blocknm,0);
preSpikes = spikes.sp;
preRunT = spd.t;
preRunV = spd.v;

blocknm = blocks{2};
spikes = getSpikes(clustfile,afile, blocknm,0);
spd = getSpeed(clustfile,afile,blocknm,0);
postSpikes = spikes.sp;
postRunT = spd.t;
postRunV = spd.v;

load(afile,'layer','cells');
% 
% if savePDF
%     psfilename = 'C:\analysisPS.ps';
%     if exist(psfilename,'file')==2;delete(psfilename);end %%%
% end

dur = min(max(preRunT),max(postRunT));
%dt = 1;
%histbins = dt/2:dt:dur

clear R cv2
for c = 1:length(preSpikes);
   % figure
    for cond = 1:2;
        if cond==1
            col = 'b';
            sp = preSpikes{c};
            preISI  = diff(sp(1:end-1));
            postISI = diff(sp(2:end));
%             subplot(2,2,4)
%             loglog(preISI*10^3,postISI*10^3,'b.'); xlabel 'pre ISI'; ylabel 'postISI'; axis square
%             hold on
        else
            col = 'r';
            sp = postSpikes{c};
            preISI  = diff(sp(1:end-1));
            postISI = diff(sp(2:end));
            %subplot(2,2,4)
            %loglog(preISI*10^3,postISI*10^3,'r.'); axis([1 10^4 1 10^4])
        end
        R(c,:,cond) = (hist(sp(sp<dur),dt/2:dt:dur))/dt;
        cv2(c,:,cond) = mean(2*abs(preISI-postISI)./(preISI+postISI));
    end
    
%     subplot(2,2,1)
%    plot(squeeze(R(c,:,:))); xlabel 'time(s)'; ylabel 'sp/sec';
    
%subplot(2,2,2)
%     bar(squeeze(cv2(c,:,:))); ylim([0 1.5])
%     ax = gca;
%     ax.XTick = [1 2];
%     ax.XTickLabels = {'pre','post'};
%     title_text = sprintf('ch %d cl %d',cells(c,1),cells(c,2));
%     text(0,-10,title_text,'FontSize',8);
%     ylabel('cv2')

%     subplot(2,2,3)
%     bar(squeeze(mean(R(c,:,:),2))); ylim([0 10]);
%     ax = gca;
%     ax.XTick = [1 2];
%     ax.XTickLabels = {'pre','post'};
%     ylabel('mean sp/sec')
    
%     if savePDF
%         set(gcf, 'PaperPositionMode', 'auto');
%         print('-dpsc',psfilename,'-append');
%     end
%    close(gcf)
end


% bar
% high CV = more variation in spiking pattern...less bursty
%cv pre vs cv post
%figure
%subplot(1,2,1)
%scatterhist(cv2(:,:,1),cv2(:,:,2),'Color',[0,.3,1]); xlabel 'CV2 Pre DOI'; ylabel 'CV2 Post DOI'; lsline
%ylim([0.6 1.6]); axis square;
% 
% if savePDF
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');
% end

clear normR
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

% if savePDF
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');
% end

rmax = 10;
figure
subplot(2,1,1);
imagesc(normR(:,:,1),[0 rmax]); title('pre'); xlabel('secs'); ylabel('cell #')
subplot(2,1,2)
imagesc(normR(:,:,2),[0 rmax]); title('post'); xlabel('secs'); ylabel('cell #')

% if savePDF
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');
% end


preCorr = corrcoef(squeeze(normR(:,:,1))');
postCorr = corrcoef(squeeze(normR(:,:,2))');
% 
% 
figure
subplot(1,2,1);
imagesc(preCorr,[-1 1]); colormap jet; title('Pre DOI'); xlabel('cell #'); ylabel('cell #'); axis square;
subplot(1,2,2);
imagesc(postCorr,[-1 1]); colormap jet; title('Post DOI');xlabel('cell #'); ylabel('cell #');
axis square;

% if savePDF
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');
% end

figure
plot(preCorr(:),postCorr(:),'.')
axis equal; hold on
plot([-0.5 0.5],[ -0.5 0.5],'r')
axis([-0.5 1 -0.5 1])
xlabel('pre correlation'); ylabel('post correlation')

figure
bar([nanmean(abs(preCorr(preCorr~=1))) nanmean(abs(postCorr(postCorr~=1)))]);
hold on
errorbar([1 2],[nanmean(abs(preCorr(preCorr~=1))) nanmean(abs(postCorr(postCorr~=1)))], [nanstd(abs(preCorr(preCorr~=1))) nanstd(abs(postCorr(postCorr~=1)))]/sqrt(0.5*length(preCorr(:))),'o')
ylabel('mean abs correlation');
set(gca,'XtickLabel',{'pre','post'});

figure
bar([nanmean((preCorr(preCorr~=1))) nanmean((postCorr(postCorr~=1)))]);
hold on
errorbar([1 2],[nanmean((preCorr(preCorr~=1))) nanmean((postCorr(postCorr~=1)))], [nanstd((preCorr(preCorr~=1))) nanstd((postCorr(postCorr~=1)))]/sqrt(0.5 *length(preCorr(:))),'o')
ylabel('mean correlation');
set(gca,'XtickLabel',{'pre','post'});


figure
cbins = -1:0.05:1;
plot(cbins,hist(preCorr(:),cbins));
hold on
plot(cbins,hist(postCorr(:),cbins),'g');
xlabel('correlation');
xlim([-0.75 0.75])

preCorrClean = preCorr; preCorrClean(isnan(preCorrClean))=0;
postCorrClean = postCorr; postCorrClean(isnan(postCorrClean))=0;

clear eigs
eigs(:,1) = sort(eig(preCorrClean),'descend');
eigs(:,2) = sort(eig(postCorrClean),'descend');
% figure
% plot(eigs); hold on; plot([1 length(eigs)],[1 1],'k:')norm
% legend({'pre','post'}); ylabel('eigenvalues');

normR(isnan(normR))=0;
[coeff score prelatent] = pca(squeeze(normR(:,:,1))');
[coeff score postlatent] = pca(squeeze(normR(:,:,2))');

latents(:,1)=prelatent/sum(prelatent);
latents(:,2)=postlatent/sum(postlatent);

% figure
% plot(latents);
% ylabel('latents'); legend('pre','post');

preRunV = preRunV(preRunT<dur); preRunT= preRunT(preRunT<dur);
postRunV = postRunV(postRunT<dur);postRunT= postRunT(postRunT<dur);

%dbstop
% figure
% subplot(2,1,1);
% plot(preRunT,preRunV); xlim([0 dur]); xlabel('secs'); ylabel('pre speed')
% subplot(2,1,2);
% plot(postRunT,postRunV);xlim([0 dur]); xlabel('secs'); ylabel('post speed')

% if savePDF
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');
% end

meanR = squeeze(mean(R,2));
figure
plot(meanR(:,1),meanR(:,2),'o');hold on; plot([0 20],[0 20]);axis equal
xlabel('pre rate'); ylabel('post rate');

% if savePDF
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');
% end

% meanRpre = mean(meanR(:,1))
% meanRpost = mean(meanR(:,2))
% barx = [meanRpre meanRpost]

% figure
% subplot(2,1,1)
% barweb(meanR(:,1),meanR(:,2))
% 
% subplot(2,1,2)
% bar(barx); title 'mean FR'; ylabel 'sp/sec'
% set(gca,'Xticklabel',{'pre','post'});
% 
% if savePDF
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');
% end

%allR=[squeeze(normR(:,1:600,1)) squeeze(normR(:,1:600,2))];

allR = [squeeze(normR(:,:,1)) squeeze(normR(:,:,2))];
allV = [preRunV postRunV];
allT = [preRunT postRunT+dur];

figure
imagesc(allR,[0 1]);
axis xy;
hold on
plot(allT/dt,5*allV/max(allV)-5 ,'g');
ylim([-5 size(allR,1)+0.5]);
%dur =600
plot([dur dur],[-5 size(allR,1)],'r','Linewidth',4); %colorbar;

%dbstop
allR = [squeeze(normR(:,:,1)) squeeze(normR(:,:,2))];
% 
% vthresh=1;
% stat = find(allV(1:10:end-10)<vthresh)
% mv = find(allV(1:10:end-10)<vthresh)
% 
% test = mean(allR(:,stat))
% mean(allR(:,mv))

 
% if savePDF
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');
% end

[coeff score latent] = pca(allR');
% figure
% plot(latent/sum(latent))
% 
% if savePDF
% set(gcf, 'PaperPositionMode', 'auto');
% print('-dpsc',psfilename,'-append');
% end

figure
subplot(2,1,1)
imagesc(coeff); axis square
subplot(2,1,2)

%hist of PCA loadings
% hist(coeff); axis square
% if savePDF
% set(gcf, 'PaperPositionMode', 'auto');
% print('-dpsc',psfilename,'-append');
% end


figure
hold on
plot(score(:,1),score(:,2))
for i=1:length(score);
    plot(score(i,1),score(i,2),'.','Markersize',12,'Color',cmapVar(i,1,length(score),jet))
end
    xlabel 'PC1'; ylabel 'PC2';
    
    

% if savePDF
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');
% end

% figure
% for i = 1:3
%     subplot(4,1,i)
%     plot(score(:,i)); ylabel(sprintf('pca %d',i))
% end
% subplot(4,1,4)
% plot(allT,allV); ylabel('speed')

% if savePDF
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');
% end
% 
% %timeponts
% figure
% imagesc(corrcoef(allR)); axis square; xlabel('secs'); ylabel('secs'); colorbar; title('correlation')
% colormap jet;
% 
% if savePDF
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');
% end
% 
% %corr units
% figure
% imagesc(corrcoef(allR'),[-1 1]); axis square; xlabel('cell #'); ylabel('cell #'); colorbar; title('correlation')
% colormap jet;

% 
% if savePDF
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');
% end

% if savePDF
%     [f p] = uiputfile('*.pdf','pdf name');
%     ps2pdf('psfile', psfilename, 'pdffile', fullfile(p,f));
%     delete(psfilename);
end