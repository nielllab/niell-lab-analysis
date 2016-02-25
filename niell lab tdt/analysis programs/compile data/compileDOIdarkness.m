clear all
close all

dbstop if error

batchDOIephys; %%% load batch file

%%% select the sessions you want based on filters
%%% example
use =  find(strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark})   )

%use =  find( strcmp({files.treatment},'DOI') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}) )

sprintf('%d selected sessions',length(use))

saline=1; doi=2; lisuride=3;

savePDF=0;
redo = 1;
n=0; ncorr=0; %%% number of units loaded, ncorr= number of correlation pairs
for i = 1:3
    
    %%% extract filenames
    afile = [pathname '\' files(use(i)).dir '\' files(use(i)).analysisfile '.mat'];
    clustfile = [pathname '\' files(use(i)).dir '\' files(use(i)).clusterfile '.mat'] ;
    
    %%% get cell type based on waveform
    [inh mid] = getWaveform(clustfile,afile,0);
    nc = length(inh); cellrange = n+1:n+nc;
    inhAll(cellrange) = inh;
    
    %%% get layer info
    load(afile,'layer');
    layerAll(cellrange) = layer;
    %%% getLayers (needs histo information, but will give layers for all sites
    
    %%% session info
    sessionNum(cellrange)=i;
    for j=1:length(cellrange); expt{cellrange(j)} = files(use(i)).expt; end
    if strcmp(files(use(i)).treatment,'Saline'), treatment(cellrange)=saline, end;
    if strcmp(files(use(i)).treatment,'DOI'), treatment(cellrange)=doi, end;
    if strcmp(files(use(i)).treatment,'Lisuride'), treatment(cellrange)=lisuride, end;
    sessionTreatment(i) = treatment(cellrange(1));

for prepost = 1:2
spd = getSpeed(clustfile,afile,files(use(i)).blockWn{prepost},0);
speedHistWn(i,:,prepost) = hist(spd.v,0.5:1:100)/length(spd.v);
end


    %%% get grating responses
    for prepost = 1:2
        drift = getDrift_mv(clustfile,afile,files(use(i)).blockDrift{prepost},0);
        drift_orient(cellrange,:,:,prepost)=drift.orient_tune;
        drift_sf(cellrange,:,:,prepost) = drift.sf_tune;
        drift_spont(cellrange,:,prepost) = drift.interSpont;
    end
    
%%% get wn response
    for prepost = 1:2      
        wn = getWn_mv(clustfile,afile,files(use(i)).blockWn{prepost},0,300);
        wn_crf(cellrange,:,:,prepost)=wn.crf;
        wn_spont(cellrange,:,prepost)=wn.spont;
        wn_evoked(cellrange,:,prepost)=wn.evoked;
    end
    
    %%% lfp power 
    %%%(right now averages over all sites, should use layer info)
    for prepost=1:2
    lfpMove = getLfpMovement(clustfile,afile,files(use(i)).blockWn{prepost},0);
    LFPall(i,:,:,prepost) =squeeze(median(lfpMove.meanSpect, 1))/median(lfpMove.meanSpect(:));
    end
    
    %%% darkness / correlation analysis

    %%% need to keep track of n^2 values for correlations
    corrRange=ncorr+1:ncorr+nc^2;
    corrTreatment(corrRange)=treatment(cellrange(1));
    
    
    %%% prepost correlation for white noise
    dt = 1;
    [preCorr postCorr cv2 R eigs] = prepostDOIdarkness(clustfile,afile,files(use(i)).blockWn,dt,0);
    wnCorr(corrRange,1) = preCorr(:); wnCorr(corrRange,2)=postCorr(:);
    
    %%%% prepost correlation in darkness
    dt = 1;
    [preCorr postCorr cv2 R eigs] = prepostDOIdarkness(clustfile,afile,files(use(i)).blockDark,dt,0);
    darkCorr(corrRange,1) = preCorr(:); darkCorr(corrRange,2)=postCorr(:);
    
    cv2Dark(cellrange,:) = cv2;
    meanRdark(cellrange,:) = mean(R,2);

    %%% keep track of cell type for correlations
    corrType1 = zeros(size(preCorr)); corrType2 = corrType1;
    for j= 1:length(inh);
        corrType1(j,:)=inh(j); corrType2(:,j)=inh(j);
    end
    corrType1all(corrRange) = corrType1(:) ; corrType2all(corrRange)= corrType2(:);

  n= n+nc;
  ncorr= ncorr+nc^2;
end

titles = {'saline','doi'};
figure
for i = 1:2
subplot(1,2,i);
plot(wnCorr(corrTreatment==i,1),wnCorr(corrTreatment==i,2),'.'); hold on; axis equal
plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
xlabel('pre wn corr'); ylabel('post')
end

titles = {'saline','doi'};
figure
for i = 1:2
subplot(1,2,i);
plot(darkCorr(corrTreatment==i,1),darkCorr(corrTreatment==i,2),'.'); hold on; axis equal
plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
xlabel('pre dark corr'); ylabel('post')
end

figure
plot(drift_spont(:),wn_spont(:),'.'); hold on; plot([0 10], [0 10]); axis equal
xlabel('drift spont'),ylabel('wn spont');


figure
for i = 1:6
    subplot(2,3,i)
    plot(drift_spont(treatment==saline & layerAll ==i,1,1),drift_spont(treatment==saline& layerAll ==i,1,2),'k.');
    hold on
    plot(drift_spont(treatment==doi& layerAll ==i,1,1),drift_spont(treatment==doi& layerAll ==i,1,2),'r.');
    plot([0 10],[0 10]); axis equal
    title(sprintf('layer %d',i)); xlabel('stop drift spont'); ylabel('post');
end

figure
for i = 1:6
    subplot(2,3,i)
    plot(drift_spont(treatment==saline & layerAll ==i,2,1),drift_spont(treatment==saline& layerAll ==i,2,2),'k.');
    hold on
    plot(drift_spont(treatment==doi& layerAll ==i,2,1),drift_spont(treatment==doi& layerAll ==i,2,2),'r.');
    plot([0 10],[0 10]); axis equal
    title(sprintf('layer %d',i)); xlabel('move drift spont'); ylabel('post');
end


figure
for i = 1:6
    subplot(2,3,i)
    plot(wn_spont(treatment==saline & layerAll ==i,1,1),wn_spont(treatment==saline& layerAll ==i,1,2),'k.');
    hold on
    plot(wn_spont(treatment==doi& layerAll ==i,1,1),wn_spont(treatment==doi& layerAll ==i,1,2),'r.');
    plot([0 10],[0 10]); axis equal
    title(sprintf('layer %d',i)); xlabel('stop wn spont'); ylabel('post');
end

figure
for i = 1:6
    subplot(2,3,i)
    plot(wn_spont(treatment==saline & layerAll ==i,2,1),wn_spont(treatment==saline& layerAll ==i,2,2),'k.');
    hold on
    plot(wn_spont(treatment==doi& layerAll ==i,2,1),wn_spont(treatment==doi& layerAll ==i,2,2),'r.');
    plot([0 10],[0 10]); axis equal
    title(sprintf('layer %d',i)); xlabel('move wn spont'); ylabel('post');
end


figure
for i = 1:6
    subplot(2,3,i)
    plot(wn_evoked(treatment==saline & layerAll ==i,1,1),wn_evoked(treatment==saline& layerAll ==i,1,2),'k.');
    hold on
    plot(wn_evoked(treatment==doi& layerAll ==i,1,1),wn_evoked(treatment==doi& layerAll ==i,1,2),'r.');
    plot([0 10],[0 10]); xl = get(gca,'Xlim'); yl = get(gca,'Ylim'); axis square; axis([min(xl(1),yl(1)) max(xl(2),yl(2)) min(xl(1),yl(1)) max(xl(2),yl(2)) ])
    title(sprintf('layer %d',i)); xlabel('stop wn evoked'); ylabel('post');
end

figure
for i = 1:6
    subplot(2,3,i)
    plot(wn_evoked(treatment==saline & layerAll ==i,2,1),wn_evoked(treatment==saline& layerAll ==i,2,2),'k.');
    hold on
    plot(wn_evoked(treatment==doi& layerAll ==i,2,1),wn_evoked(treatment==doi& layerAll ==i,2,2),'r.');
    plot([0 10],[0 10]); xl = get(gca,'Xlim'); yl = get(gca,'Ylim'); axis square; axis([min(xl(1),yl(1)) max(xl(2),yl(2)) min(xl(1),yl(1)) max(xl(2),yl(2)) ])
    title(sprintf('layer %d',i)); xlabel('move wn evoked'); ylabel('post');
end

for t = 1:2
figure
    useN = find(treatment==t)
    for i = 1:length(useN)
        np = ceil(sqrt(length(useN)));
        subplot(np,np,i);
        hold on
        plot(wn_crf(useN(i),:,1,1),'Color',[0.5 0 0]);  plot(wn_crf(useN(i),:,2,1),'Color',[0 0.5 0]);
        plot(wn_crf(useN(i),:,1,2),'Color',[1 0 0]);  plot(wn_crf(useN(i),:,2,2),'Color',[0 1 0]);
        yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)])
    end
end


for t = 1:2
figure
    useN = find(treatment==t)
    for i = 1:length(useN)
        np = ceil(sqrt(length(useN)));
        subplot(np,np,i);
        hold on
        plot(drift_orient(useN(i),:,1,1),'Color',[0.5 0 0]);  plot(drift_orient(useN(i),:,2,1),'Color',[0 0.5 0]);
        plot(drift_orient(useN(i),:,1,2),'Color',[1 0 0]);  plot(drift_orient(useN(i),:,2,2),'Color',[0 1 0]);
             plot([1 12], [1 1]*drift_spont(useN(i),1,1),':','Color',[0.5 0 0]);  plot([1 12], [1 1]*drift_spont(useN(i),2,1),':','Color',[0 0.5 0]);
        plot([1 12], [1 1]*drift_spont(useN(i),1,2),':','Color',[1 0 0]);  plot([1 12], [1 1]*drift_spont(useN(i),2,2),':','Color',[0 1 0]);
yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)]); xlim([0.5 12.5])
    end
end

figure
hold on
subplot(1,2,1)
plot(mean(speedHistWn(sessionTreatment==saline,:,:),1)'); title('saline'); xlabel('speed')
subplot(1,2,2)
plot(mean(speedHistWn(sessionTreatment==doi,:,:),1)'); title('doi'); xlabel('speed')
