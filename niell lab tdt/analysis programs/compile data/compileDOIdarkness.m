clear all
close all

dbstop if error

batchDOIephys; %%% load batch file

%%% select the sessions you want based on filters
%%% example
use =  find(strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark})& ~cellfun(@isempty,{files.postdark}) )

%use =  find( strcmp({files.treatment},'DOI') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}) )
%use =  find( strcmp({files.treatment},'Saline') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}) & strcmp({files.expt},'040116'))


sprintf('%d selected sessions',length(use))

saline=1; doi=2; lisuride=3;

savePDF=0;
redo = 1;
n=0; ncorr=0; %%% number of units loaded, ncorr= number of correlation pairs
for i = 8:10 %length(use)
    
    %%% extract filenames
    afile = [pathname '\' files(use(i)).dir '\' files(use(i)).analysisfile '.mat'];
    clustfile = [pathname '\' files(use(i)).dir '\' files(use(i)).clusterfile '.mat'] ;
    
    %%% get cell type based on waveform
    [inh mid] = getWaveform(clustfile,afile,0);
    nc = length(inh); cellrange = n+1:n+nc;
    inhAll(cellrange) = inh;
    
    %%% get layer info
   % load(afile,'layer');
   % layerAll(cellrange) = layer;
    %%% getLayers (needs histo information, but will give layers for all sites)
    
    %%% getEyes  (needs camera files)d
    
    %%% session info
    sessionNum(cellrange)=i;
    for j=1:length(cellrange); expt{cellrange(j)} = files(use(i)).expt; end
    if strcmp(files(use(i)).treatment,'Saline'), treatment(cellrange)=saline, end;
    if strcmp(files(use(i)).treatment,'DOI'), treatment(cellrange)=doi, end;
    if strcmp(files(use(i)).treatment,'Lisuride'), treatment(cellrange)=lisuride, end;
    sessionTreatment(i) = treatment(cellrange(1));
    
    %%% get pre/post running speed
    for prepost = 1:2
        spd = getSpeed(clustfile,afile,files(use(i)).blockWn{prepost},0);
        speedHistWn(i,:,prepost) = hist(spd.v,0.5:1:100)/length(spd.v);
        speedTrace{i,prepost}=spd.v;
    end
    
    
    %%% get grating responses
    for prepost = 1:2
        drift = getDrift_mv(clustfile,afile,files(use(i)).blockDrift{prepost},1);
        drift_orient(cellrange,:,:,prepost)=drift.orient_tune;
        drift_sf(cellrange,:,:,prepost) = drift.sf_tune;
        drift_spont(cellrange,:,prepost) = drift.interSpont;
        drift_osi(cellrange,:,prepost) = drift.cv_osi;
        drift_F1F0(:,cellrange,prepost)= drift.F1F0;
        %drift_ot_tune(cellrange,:,prepost)=drift.orient_tune;
        %drift_sf_trial(1,:,prepost) = drift.trialSF;
        %drift_trial_psth(cellrange,:,:,prepost) = drift.trialPsth;
    end
    
    %%% get wn response
    for prepost = 1:2
        wn = getWn_mv(clustfile,afile,files(use(i)).blockWn{prepost},0,300);
        wn_crf(cellrange,:,:,prepost)=wn.crf;
        wn_spont(cellrange,:,prepost)=wn.spont;
        wn_evoked(cellrange,:,prepost)=wn.evoked;
        wn_gain(1,cellrange,prepost)=wn.gain
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

%%% plot correlation for white noise
titles = {'saline','doi'};
figure
for i = 1:2
    subplot(1,2,i);
    plot(wnCorr(corrTreatment==i,1),wnCorr(corrTreatment==i,2),'.'); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
    xlabel('pre wn corr'); ylabel('post')
end


%%% plot correlation for darkness
titles = {'saline','doi'};
figure
for i = 1:2
    subplot(1,2,i);
    plot(darkCorr(corrTreatment==i,1),darkCorr(corrTreatment==i,2),'.'); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
    xlabel('pre dark corr'); ylabel('post')
end

%%% compare spontaneous rates measured with gratings and wn
figure
plot(drift_spont(:),wn_spont(:),'.'); hold on; plot([0 10], [0 10]); axis equal
xlabel('drift spont'),ylabel('wn spont');

%%% scatter plot of drift spont
for mv = 1:2
    figure
    for i = 1:6
        subplot(2,3,i)
        plot(drift_spont(treatment==saline & layerAll ==i,mv,1),drift_spont(treatment==saline& layerAll ==i,mv,2),'k.');
        hold on
        plot(drift_spont(treatment==doi& layerAll ==i,mv,1),drift_spont(treatment==doi& layerAll ==i,mv,2),'r.');
        plot([0 10],[0 10]); axis equal
        title(sprintf('layer %d',i)); ylabel('post');
        if mv ==1 , xlabel('stop drift spont'); else  xlabel('move drift spont'); end
        
    end
end

%%% scatter plot of wn spont
for mv = 1:2
    figure
    for i = 1:6
        subplot(2,3,i)
        plot(wn_spont(treatment==saline & layerAll ==i,mv,1),wn_spont(treatment==saline& layerAll ==i,mv,2),'k.');
        hold on
        plot(wn_spont(treatment==doi& layerAll ==i,mv,1),wn_spont(treatment==doi& layerAll ==i,mv,2),'r.');
        plot([0 10],[0 10]); axis equal
        title(sprintf('layer %d',i));  ylabel('post');
        if mv ==1 , xlabel('stop wn spont'); else  xlabel('move wn spont'); end
    end
end


%%% scatter plot wn evoked
for mv = 1:2
    figure
    for i = 1:6
        subplot(2,3,i)
        plot(wn_evoked(treatment==saline & layerAll ==i,mv,1),wn_evoked(treatment==saline& layerAll ==i,mv,2),'k.');
        hold on
        plot(wn_evoked(treatment==doi& layerAll ==i,mv,1),wn_evoked(treatment==doi& layerAll ==i,mv,2),'r.');
        plot([0 10],[0 10]); xl = get(gca,'Xlim'); yl = get(gca,'Ylim'); axis square; axis([min(xl(1),yl(1)) max(xl(2),yl(2)) min(xl(1),yl(1)) max(xl(2),yl(2)) ])
        title(sprintf('layer %d',i));  ylabel('post');
        if mv ==1 , xlabel('stop wn evoked'); else  xlabel('move wn evoked'); end
    end
end

%%% plot white noise response functions for all units
for t = 1:2
    figure
    if t==1, set(gcf,'Name','saline wn CRF'), else set(gcf,'Name','doi wn CRF'),end
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

%%% plot orientation tuning curves for all units
for t = 1:2
    figure
    if t==1, set(gcf,'Name','saline drift orientation'), else set(gcf,'Name','doi drift orientation'),end
    useN = find(treatment==t)
    for i = 1:length(useN)
        np = ceil(sqrt(length(useN)));
        subplot(np,np,i);
        hold on
        plot(drift_orient(useN(i),:,1,1),'Color',[0.5 0 0]);  plot(drift_orient(useN(i),:,2,1),'Color',[0 0.5 0]); %pre sal & DOI mv & stat
        plot(drift_orient(useN(i),:,1,2),'Color',[1 0 0]);  plot(drift_orient(useN(i),:,2,2),'Color',[0 1 0]); %pre sal & DOI mv & stat
        plot([1 12], [1 1]*drift_spont(useN(i),1,1),':','Color',[0.5 0 0]);  plot([1 12], [1 1]*drift_spont(useN(i),2,1),':','Color',[0 0.5 0]);
        plot([1 12], [1 1]*drift_spont(useN(i),1,2),':','Color',[1 0 0]);  plot([1 12], [1 1]*drift_spont(useN(i),2,2),':','Color',[0 1 0]);
        yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)]); xlim([0.5 12.5])
    end
end

%plot spatial frequency tuning curves for all units
for t = 1:2
    figure
    if t==1, set(gcf,'Name','saline spatial frequency'), else set(gcf,'Name','doi spatial frequency'),end
    useN = find(treatment==t)
    for i = 1:length(useN)
        np = ceil(sqrt(length(useN)));
        subplot(np,np,i);
        hold on
        plot(drift_sf(useN(i),:,1,1),'Color',[0.5 0 0]);  plot(drift_sf(useN(i),:,2,1),'Color',[0 0.5 0]); %pre, mv=2
        plot(drift_sf(useN(i),:,1,2),'Color',[1 0 0]);plot(drift_sf(useN(i),:,2,2),'Color',[0 1 0]); %post
        plot([1 7], [1 1]*drift_spont(useN(i),1,1),':','Color',[0.5 0 0]);  plot([1 7], [1 1]*drift_spont(useN(i),2,1),':','Color',[0 0.5 0]);
        plot([1 7], [1 1]*drift_spont(useN(i),1,2),':','Color',[1 0 0]);  plot([1 7], [1 1]*drift_spont(useN(i),2,2),':','Color',[0 1 0]);
        yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)]); xlim([0.5 8.5])
    end
end

%%% plot speed histogram
figure
hold on
subplot(2,1,1)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==saline,1:25,:),1))); title('saline'); xlabel('speed')
subplot(2,1,2)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==doi,1:25,:),1))); title('doi'); xlabel('speed')
legend('pre','post')

%%%%sf pre and post%%%
%extract low and high SF responses for doi drift - mv
% this is averaging SF responses over all trials/condition -- do not use
% for decoding
% for t = 2
%     useN = find(treatment==t)
%     for i = 1:length(useN)
%         np = ceil(sqrt(length(useN)));
% lowsf_predoi(i) = (drift_sf(useN(i),2,1,1))
% hsf_predoi(i) = (drift_sf(useN(i),6,1,1))
% 
% lowsf_postdoi(i) = (drift_sf(useN(i),2,1,2))
% hsf_postdoi(i) = (drift_sf(useN(i),6,1,2))
%     end
% end
% % 
% data_pre = [lowsf_predoi; hsf_predoi]'
% 
% %prepost low sf
% figure
% plot(lowsf_predoi, '.b'); hold on; plot(lowsf_postdoi,'.r')
% 
% %prepost high sf
% figure
% plot(hsf_predoi, '.b'); hold on; plot(hsf_postdoi,'.r')
% 
% 
% data_sf_pre = [drift_sf(:,2,1,1); drift_sf(:,6,1,1)]
% %lowsf = data_sf_pre(1:106); hsf=data_sf_pre(107:end)
% sf(1:106)=1; sf(107:212) =2;
% 
% pre_full = []
% pre_full = [data_sf_pre  sf'];
% 
% figure  
% plot(data_sf_pre)
% figure
% plot(sf)

lowsf= find(drift.trialSF==2)
highsf=find(drift.trialSF==6)
%to separate high/low SF responses for classifier
%need to separate prepost treatments
% also separate treatment 
% for t=1:2
%     useN = find(treatment==t)
%     for i = 1:length(useN)
%         np = ceil(sqrt(length(useN)));
%         
%         resp_lowsf_pre= mean(drift.trialPsth(cellrange,lowsf,:,1),3) %each cell's response to lowsf, pre
%         resp_highsf_pre= mean(drift.trialPsth(useN(i),highsf,:,1),3)
%         resp_sf_pre= [resp_lowsf_pre;resp_highsf_pre]
%         
%         sf=zeros([1 length(resp_sf_pre)])
%         sf(1:length(resp_sf_pre)/2)=1
%         sf(sf~=1)=2
%         
%         data_full_sf_pre = [resp_sf_pre sf']
%     end
% end
% 
% 
% for t=1:2
%     useN = find(treatment==t)
%     for i = 1:length(useN)
%         np = ceil(sqrt(length(useN)));
%         
%         resp_lowsf_post= mean(drift_trial_psth(cellrange,lowsf,:,2),3) %each cell's response to lowsf, post
%         resp_highsf_post= mean(drift_trial_psth(cellrange,highsf,:,2),3)
%         resp_sf_post= [resp_lowsf_post;resp_highsf_post]
%         
%         sf=zeros([1 length(resp_sf_post)])
%         sf(1:length(resp_sf_post)/2)=1
%         sf(sf~=1)=2
%         data_full_sf_post = [resp_sf_post sf']
%     end
% end
%  
% figure
% subplot(1,2,1)
% plot(data_full_sf_pre);title 'Pre SF responses';
% subplot(1,2,2)
% plot(data_full_sf_post);title 'Post SF responses';

%to separate orientation responses for classifier
% orient90= find(drift.trialOrient==4)
% orient270=find(drift.trialOrient==10)
% for t=1:2
%     useN = find(treatment==t)
%     for i = 1:length(useN)
%         np = ceil(sqrt(length(useN)));
%         
%         resp_90= mean(drift.trialPsth(:,orient90,:),3);
%         resp_270= mean(drift.trialPsth(:,orient270,:),3);
%         resp_orient= [resp_90;resp_270];
%         
%         
%         ot=zeros([1 length(resp_orient)]);
%         ot(1:length(resp_orient)/2)=1;
%         ot(ot~=1)=2;
%         
%         data_full_orient = [resp_orient ot'];
%         
%     end
% end
% % 
% %%%to set up trials with orientation
% orient90= find(drift.trialOrient==4)
% orient270=find(drift.trialOrient==10)
% 
% resp_90= mean(drift.trialPsth(:,orient90,:),3)
% resp_270= mean(drift.trialPsth(:,orient270,:),3)
% resp_orient= [resp_90;resp_270] 
% ot(1:56)=1;ot(56:112)=2
% data_full_orient = [resp_orient ot']
