clear all
close all 

dbstop if error

batchDOIephys; %%% load batch file

%%% select the sessions you want based on filters
%%% example
use =  find(strcmp({files.notes},'good data')& ~cellfun(@isempty,{files.predark})& ~cellfun(@isempty,{files.postdark}) )
%use =  find( strcmp({files.treatment},'KetanserinDOI') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}) )
%use =  find( strcmp({files.treatment},'DOI') &  ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}))

%for specific experiment:
%use =  find(strcmp({files.notes},'bad data')  & ~cellfun(@isempty,{files.predark})& ~cellfun(@isempty,{files.postdark}) & strcmp({files.expt},'120915'))
sprintf('%d selected sessions',length(use))

saline=1; doi=2; ketanserin=3; ketandoi=4; lisuride=5;
 
savePDF=0;
redo = 1;
n=0; ncorr=0; %%% number of units loaded, ncorr= number of correlation pairs
for i = 1:length(use)
    close all
    %%% extract filenames
    afile = [pathname '\' files(use(i)).dir '\' files(use(i)).analysisfile '.mat'];
    clustfile = [pathname '\' files(use(i)).dir '\' files(use(i)).clusterfile '.mat'] ;
    cfile = [{[pathname '\' files(use(i)).dir '\' files(use(i)).predark_camera '.mat']}; {[pathname '\' files(use(i)).dir '\' files(use(i)).postdark_camera '.mat']}]';
 
    %%% get cell type based on waveform
    [inh mid] = getWaveform(clustfile,afile,1);
    nc = length(inh); cellrange = n+1:n+nc;
    inhAll(cellrange) = inh;
    
    
    if ~isempty(files(use(i)).prepinpFile)
        pinp =getPinpHTR2A(clustfile,afile, files(use(i)).prepinpFile,3,5,0,0,1)
        pinpAll(cellrange,1) = pinp.pinped;
%         pinpPsth(cellrange,:,1)=pinp.psth;
    else
        pinpAll(cellrange,1)=NaN;
%         pinpPsth(cellrange,:,1)=NaN;
    end
    
    if ~isempty(files(use(i)).postpinpFile)
        pinp =getPinpHTR2A(clustfile,afile, files(use(i)).postpinpFile,3,5,0,0,1)
        pinpAll(cellrange,2) = pinp.pinped;
%         pinpPsth(cellrange,:,2)=pinp.psth;
    else
        pinpAll(cellrange,2)=NaN;
%        pinpPsth(cellrange,:,2)=NaN;
    end
    
    
    %% get layer info
    clear layer
    load(afile,'layer');
    if exist('layer','var')
        layerAll(cellrange) = layer;
    else
        layerInfo = getLayer(clustfile,afile,files(use(i)).tip1,files(use(i)).tip2,files(use(i)).angle, 1); %(needs histo information, but will give layers for all sites)
        layerAll(cellrange) = layerInfo.units;
        layerSites= layerInfo.sites;
    end
    
if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2}) 
    
for prepost =1:2
eyes = getEyes(clustfile,afile,cfile{:,prepost}, files(use(i)).blockWn{prepost},1);
rad{:,prepost} = eyes.rad
t{prepost} = eyes.t
end
end

    %%% session info
    sessionNum(cellrange)=i;
    for j=1:length(cellrange); expt{cellrange(j)} = files(use(i)).expt; end
    if strcmp(files(use(i)).treatment,'Saline'), treatment(cellrange)=saline, end;
    if strcmp(files(use(i)).treatment,'DOI'), treatment(cellrange)=doi, end;
    if strcmp(files(use(i)).treatment,'Ketanserin'), treatment(cellrange)=ketanserin, end;
    if strcmp(files(use(i)).treatment,'KetanserinDOI'), treatment(cellrange)=ketandoi, end;
    if strcmp(files(use(i)).treatment,'Lisuride'), treatment(cellrange)=lisuride, end;

    sessionTreatment(i) = treatment(cellrange(1));
    
      if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2})  
    %% get pre/post running speed
    for prepost = 1:2
        spd = getSpeed(clustfile,afile,files(use(i)).blockWn{prepost},1);
        speedHistWn(i,:,prepost) = hist(spd.v,0.5:1:100)/length(spd.v);
        speedTrace{i,prepost}=spd.v;
    end
      end
      
    if ~isempty(files(use(i)).blockDrift{1}) & ~isempty(files(use(i)).blockDrift{2})
        %%% get grating responses
        for prepost = 1:2
            drift = getDrift_mv(clustfile,afile,files(use(i)).blockDrift{prepost},0);
            drift_orient(cellrange,:,:,prepost)=drift.orient_tune;
            drift_sf(cellrange,:,:,prepost) = drift.sf_tune;
            drift_spont(cellrange,:,prepost) = drift.interSpont;
            drift_osi(cellrange,:,prepost) = drift.cv_osi;
            drift_F1F0(:,cellrange,prepost)= drift.F1F0;
            drift_ot_tune(cellrange,:,:,prepost)=drift.orient_tune;
            drift_sf_trial{prepost} = drift.trialSF;
            drift_orient_trial{prepost} = drift.trialOrient;
            for c = 1:length(cellrange)
                drift_trial_psth{cellrange(c),:,prepost} = squeeze(drift.trialPsth(c,:,:));
            end
        end
    end
    
        if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2})

    %%% get wn response
    for prepost = 1:2
        wn = getWn_mv(clustfile,afile,files(use(i)).blockWn{prepost},0,300);
        wn_crf(cellrange,:,:,prepost)=wn.crf;
        wn_spont(cellrange,:,prepost)=wn.spont;
        wn_evoked(cellrange,:,prepost)=wn.evoked;
        wn_gain(1,cellrange,prepost)=wn.gain
    end
        end
    
    %%% lfp power
%     %%%(right now averages over all sites, should use layer info)
    for prepost=1:2
        lfpMove = getLfpMovement(clustfile,afile,files(use(i)).blockWn{prepost},0);
        LFPall(i,:,:,prepost) =squeeze(nanmedian(lfpMove.meanSpect, 1))/nanmedian(lfpMove.meanSpect(:));
    end
    
      for prepost=1:2
        lfpMoveDark = getLfpMovement(clustfile,afile,files(use(i)).blockDark{prepost},0);
        LFPallDark(i,:,:,prepost) =squeeze(nanmedian(lfpMove.meanSpect, 1))/nanmedian(lfpMove.meanSpect(:));
   
      end
    
    %% darkness / correlation analysis
    
    %%% need to keep track of n^2 values for correlations
    corrRange=ncorr+1:ncorr+nc^2;
    corrTreatment(corrRange)=treatment(cellrange(1));
    
    
    %%% prepost correlation for white noise
    dt = 1;
    [preCorr postCorr cv2 R eigs] =  prepostDOIdarkness(clustfile,afile,files(use(i)).blockWn,dt,0);
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

%%eyetracking

% figure % radius histogram
% h1=hist(rad{1},1:1: max(rad{1}));
% h2=hist(rad{2},1:1: max(rad{2}));
% bins1=1:1:max(rad{1})
% bins2=1:1:max(rad{2})
% plot(bins1,h1/sum(h1)); hold on;
% plot(bins2,h2/sum(h2))
% ylim([0 .8]); xlim([0 30])
% xlabel('pixels');ylabel('proportion of time')
% title ('DOI')
% % 
% % 
% figure %raw trace
% plot(rad{1});hold on;plot(rad{2});xlim([0 4000]);ylim([0 20]);


%%% plot correlation for white noise
titles = {'saline','doi','ketanserin', 'ketanserin + DOI'};
figure
for i = 1:4
    subplot(2,2,i);
    plot(wnCorr(corrTreatment==i,1),wnCorr(corrTreatment==i,2),'.'); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
    xlabel('pre wn corr'); ylabel('post')
    set(gcf,'Name','Wn Corr')
end

titles = {'saline','doi','ketanserin', 'ketanserin + DOI'};
figure
for i = 1:4
    subplot(2,2,i);
    wnCorrHist= myHist2(wnCorr(corrTreatment==i,1),wnCorr(corrTreatment==i,2),-.5:.1:1.5,-.5:.1:1.5);
    %wnCorrHist_pre= myHist2(wnCorr(corrTreatment==i,1),-.5:.1:1.5,-.5:.1:1.5)
    plot(wnCorrHist);hold on; axis square;title(titles{i})
    set(gcf,'Name','Wn Corr'); ylim([0 6000])
end


%%% plot correlation for darkness
titles = {'saline','doi','ketanserin', 'ketanserin + DOI'};
figure
for i = 1:4
    subplot(2,2,i);
    plot(darkCorr(corrTreatment==i,1),darkCorr(corrTreatment==i,2),'.'); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
    xlabel('pre dark corr'); ylabel('post')
    set(gcf,'Name','Dark Corr')
end

titles = {'saline','doi','ketanserin', 'ketanserin + DOI'};
figure
for i = 1:4
    subplot(2,2,i);
    darkCorrHist= myHist2(darkCorr(corrTreatment==i,1),darkCorr(corrTreatment==i,2),-.5:.1:1.5,-.5:.1:1.5);
    plot(darkCorrHist);hold on; axis square;title(titles{i});
    set(gcf,'Name','Dark Corr'); ylim ([0 6000])

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
        plot(drift_spont(treatment==ketanserin& layerAll ==i,mv,1),drift_spont(treatment==ketanserin& layerAll ==i,mv,2),'m.');
        plot(drift_spont(treatment==ketandoi& layerAll ==i,mv,1),drift_spont(treatment==ketandoi& layerAll ==i,mv,2),'c.');
        plot([0 10],[0 10]); axis equal
        title(sprintf('layer %d',i)); ylabel('post');
        if mv ==1 , xlabel('stop drift spont'), set(gcf, 'Name', 'prepost stationary drift spont');
        else  xlabel('move drift spont'),set(gcf, 'Name', 'prepost move drift spont'); end
          
    end
end

%LFP averaging over all sites...separate by layer

%%evoked LFP all
% figure
% for t=1:4
% subplot(2,2,t)
% set(gcf, 'Name', 'evoked LFP')
% plot(squeeze(LFPall(:,:,1,1)),'b');hold on;
% plot(squeeze(LFPall(:,:,1,2)),'r'); xlabel 'Frequency (Hz)'; ylabel 'normalized power';
% end
% 
% figure
% set(gcf, 'Name', 'darkness LFP')
% % for t=1:4
% % subplot(2,2,t)
% plot(squeeze(LFPallDark(:,:,1,1)),'b');hold on;
% plot(squeeze(LFPallDark(:,:,1,2)),'r'); xlabel 'Frequency (Hz)'; ylabel 'normalized power';
% % end

%%% scatter plot of wn spont
for mv = 1:2
    figure
    for i = 1:6
        subplot(2,3,i)
        plot(wn_spont(treatment==saline & layerAll ==i,mv,1),wn_spont(treatment==saline& layerAll ==i,mv,2),'k.');
        % ylim ([-2 50]); xlim([-2 50]);
        hold on
        plot(wn_spont(treatment==doi& layerAll ==i,mv,1),wn_spont(treatment==doi& layerAll ==i,mv,2),'r.');
        plot(wn_spont(treatment==ketanserin& layerAll ==i,mv,1),wn_spont(treatment==ketanserin& layerAll ==i,mv,2),'m.');
        plot(wn_spont(treatment==ketandoi& layerAll ==i,mv,1),wn_spont(treatment==ketandoi& layerAll ==i,mv,2),'c.');
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
        plot(wn_evoked(treatment==ketanserin& layerAll ==i,mv,1),wn_evoked(treatment==ketanserin& layerAll ==i,mv,2),'m.');
        plot(wn_evoked(treatment==ketandoi& layerAll ==i,mv,1),wn_evoked(treatment==ketandoi& layerAll ==i,mv,2),'c.');
        plot([0 10],[0 10]); xl = get(gca,'Xlim'); yl = get(gca,'Ylim'); axis square; %axis([min(xl(1),yl(1)) max(xl(2),yl(2)) min(xl(1),yl(1)) max(xl(2),yl(2)) ])
        title(sprintf('layer %d',i));  ylabel('post'); ylim([-5 50]);xlim([-5 50])
        if mv ==1 , xlabel('stop wn evoked'); else  xlabel('move wn evoked'); end
    end
end

%rows = layers col = treatment...wn evoked FR mv and stationary
titles = {'saline','doi','ketanserin', 'ketanserin + DOI'};
for mv = 1:2
    figure
    if mv ==1 , set(gcf,'Name','stop wn evoked'); else  set(gcf,'Name','move wn evoked'); end
    for c=0:1
    for t=1:4
        subplot(6,4,t)
        plot(wn_evoked(treatment==t & layerAll ==1 & inhAll==c,mv,1),wn_evoked(treatment==t& layerAll ==1& inhAll==c,mv,2),'.');axis equal;ylim([-5 20]);xlim([-5 20])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+4)
        plot(wn_evoked(treatment==t & layerAll ==2& inhAll==c,mv,1),wn_evoked(treatment==t& layerAll ==2& inhAll==c,mv,2),'.');axis equal;ylim([-5 30]);xlim([-5 30])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+8)
        plot(wn_evoked(treatment==t & layerAll ==3& inhAll==c,mv,1),wn_evoked(treatment==t& layerAll ==3& inhAll==c,mv,2),'.');axis equal;ylim([-5 40]);xlim([-5 40])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+12)
        plot(wn_evoked(treatment==t & layerAll ==4& inhAll==c,mv,1),wn_evoked(treatment==t& layerAll ==4& inhAll==c,mv,2),'.');axis equal;ylim([-5 40]);xlim([-5 40])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+16)
        plot(wn_evoked(treatment==t & layerAll ==5& inhAll==c,mv,1),wn_evoked(treatment==t& layerAll ==5& inhAll==c,mv,2),'.');axis equal;ylim([-5 30]);xlim([-5 30])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+20)
        plot(wn_evoked(treatment==t & layerAll ==6& inhAll==c,mv,1),wn_evoked(treatment==t& layerAll ==6& inhAll==c,mv,2),'.');axis equal; ylim([-5 20]);xlim([-5 20]);hold on
        plot([0 45],[0 45]); %axis([min(xl(1),yl(1)) max(xl(2),yl(2)) min(xl(1),yl(1)) max(xl(2),yl(2)) ])
        ylabel('post'); xlabel('pre');
    end 
end
end

titles = {'saline','doi','ketanserin', 'ketanserin + DOI'};
for mv = 1:2
    figure
    if mv ==1 , set(gcf,'Name','stop wn evoked'); else  set(gcf,'Name','move wn evoked'); end
    for c=0:1
    for t=1:4
        subplot(6,4,t)
        plot(wn_evoked(treatment==t & layerAll ==1 & inhAll==c,mv,1),wn_evoked(treatment==t& layerAll ==1& inhAll==c,mv,2),'.');axis equal;ylim([-5 45]);xlim([-5 45])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+4)
        plot(wn_evoked(treatment==t & layerAll ==2& inhAll==c,mv,1),wn_evoked(treatment==t& layerAll ==2& inhAll==c,mv,2),'.');axis equal;ylim([-5 45]);xlim([-5 45])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+8)
        plot(wn_evoked(treatment==t & layerAll ==3& inhAll==c,mv,1),wn_evoked(treatment==t& layerAll ==3& inhAll==c,mv,2),'.');axis equal;ylim([-5 45]);xlim([-5 45])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+12)
        plot(wn_evoked(treatment==t & layerAll ==4& inhAll==c,mv,1),wn_evoked(treatment==t& layerAll ==4& inhAll==c,mv,2),'.');axis equal;ylim([-5 45]);xlim([-5 45])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+16)
        plot(wn_evoked(treatment==t & layerAll ==5& inhAll==c,mv,1),wn_evoked(treatment==t& layerAll ==5& inhAll==c,mv,2),'.');axis equal;ylim([-5 45]);xlim([-5 45])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+20)
        plot(wn_evoked(treatment==t & layerAll ==6& inhAll==c,mv,1),wn_evoked(treatment==t& layerAll ==6& inhAll==c,mv,2),'.');axis equal; ylim([-5 45]);xlim([-5 45]);hold on
        plot([0 45],[0 45]); %axis([min(xl(1),yl(1)) max(xl(2),yl(2)) min(xl(1),yl(1)) max(xl(2),yl(2)) ])
        ylabel('post'); xlabel('pre');
    end 
end
end



%dark layers and treatments
    figure
    for c=0:1
    for t=1:4
        subplot(6,4,t)
        plot(meanRdark(treatment==t & layerAll ==1 & inhAll==c,1),meanRdark(treatment==t& layerAll ==1& inhAll==c,2),'.');axis equal;ylim([0 20]);xlim([0 20])
        hold on; plot([0 45],[0 45]);
        subplot(6,4,t+4)
        plot(meanRdark(treatment==t & layerAll ==2& inhAll==c,1),meanRdark(treatment==t& layerAll ==2& inhAll==c,2),'.');axis equal;ylim([0 20]);xlim([0 20])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+8)
        plot(meanRdark(treatment==t & layerAll ==3& inhAll==c,1),meanRdark(treatment==t& layerAll ==3& inhAll==c,2),'.');axis equal;ylim([0 40]);xlim([0 40])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+12)
        plot(meanRdark(treatment==t & layerAll ==4& inhAll==c,1),meanRdark(treatment==t& layerAll ==4& inhAll==c,2),'.');axis equal;ylim([0 45]);xlim([0 45])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+16)
        plot(meanRdark(treatment==t & layerAll ==5& inhAll==c,1),meanRdark(treatment==t& layerAll ==5& inhAll==c,2),'.');axis equal;ylim([0 45]);xlim([0 45])
        hold on; plot([0 45],[0 45])
        subplot(6,4,t+20)
        plot(meanRdark(treatment==t & layerAll ==6& inhAll==c,1),meanRdark(treatment==t& layerAll ==6& inhAll==c,2),'.');axis equal; ylim([0 20]);xlim([0 20]);hold on
        plot([0 45],[0 45]); %axis([min(xl(1),yl(1)) max(xl(2),yl(2)) min(xl(1),yl(1)) max(xl(2),yl(2)) ])
        ylabel('post'); xlabel('pre');
    end 
    end
    

titles = {'Saline','DOI','Ketanserin', 'Ketanserin + DOI'};
%scatter all units/treatment prepost
for mv=1:2
    figure
    if mv==1 set(gcf,'Name','stationary'), else set(gcf,'Name','moving');end
   for c=0:1
    for t=1:4
        subplot(2,2,t)
        plot(wn_evoked(treatment==t & inhAll ==c,mv,1),wn_evoked(treatment==t & inhAll==c ,mv,2),'.');hold on; plot([0 45],[0 45]);axis equal; ylim([-5 40]);xlim([-5 40]);
        % plot(wn_evoked(treatment==t & inhAll(useN(i)) ==1,mv,1),wn_evoked(treatment==t &inhAll(useN(i)) ==1 ,mv,2),'r.');hold on; plot([0 10],[0 10]);axis equal; ylim([-5 40])
        ylabel('post'); xlabel('pre')
        hold on ;title(titles{t});
    end
   end
end
   


figure
for t = 1:4
    for c=0:1
subplot(2,2,t)
plot(meanRdark(treatment==t & inhAll==c,1),meanRdark(treatment==t & inhAll==c,2),'.');hold on; plot([0 45],[0 45]);axis equal;ylim([0 45]); xlim([0 45]);
title(titles{t});ylabel('post'); xlabel('pre')
end
end


%%% plot white noise response functions for all units
for t = 1:4
    figure
    if t==1, set(gcf,'Name','saline wn CRF'),
    elseif t==2, set(gcf,'Name','doi wn CRF'),
    elseif t==3, set(gcf,'Name','ketanserin wn CRF')
    else set(gcf,'Name','ketanserin + doi wn CRF'),end
    
    useN = find(treatment==t)
    for i = 1:ceil(length(useN))
        np = ceil(sqrt(length(useN)));
        subplot(np,np,i);
        hold on
        plot(wn_crf(useN(i),:,1,1),'Color',[0.5 0 0]);  plot(wn_crf(useN(i),:,2,1),'Color',[0 0.5 0]);
        plot(wn_crf(useN(i),:,1,2),'Color',[1 0 0]);  plot(wn_crf(useN(i),:,2,2),'Color',[0 1 0]);
        yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)]);
        if inhAll(useN(i)) ==1 , xlabel('inh'); else  xlabel('exc');
        end
    end
end

%%% plot orientation tuning curves for all units
for t = 1:4
      figure
    if t==1, set(gcf,'Name','saline OT'),
    elseif t==2, set(gcf,'Name','doi OT'),
    elseif t==3, set(gcf,'Name','ketanserin OT')
    else set(gcf,'Name','ketanserin + doi OT'),end
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
         if inhAll(useN(i)) ==1 , xlabel('inh'); else  xlabel('exc'); 
        end 
    end
end

%plot spatial frequency tuning curves for all units
for t = 1:4
    figure
    figure
    if t==1, set(gcf,'Name','saline SF'),
    elseif t==2, set(gcf,'Name','doi SF'),
    elseif t==3, set(gcf,'Name','ketanserin SF')
    else set(gcf,'Name','ketanserin + doi SF'),end
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
        if inhAll(useN(i)) ==1 , xlabel('inh'); else  xlabel('exc');
        end
    end
end



for t = 1:4
    figure
    figure
    if t==1, set(gcf,'Name','saline ev LFP'),
    elseif t==2, set(gcf,'Name','doi ev LFP'),
    elseif t==3, set(gcf,'Name','ketanserin ev LFP')
    else set(gcf,'Name','ketanserin + doi ev LFP'),end
    useN = find(treatment==t)
    for i = 1:length(useN)
        np = ceil(sqrt(length(useN)));
        subplot(np,np,i);
        hold on
        plot(LFPall(useN(i),:,1,1),'Color',[0.5 0 0]);  plot(LFPall(:,:,2,1),'Color',[0 0.5 0]); %pre, mv=2
        plot(LFPall(i,:,1,2),'Color',[1 0 0]);plot(LFPall(i,:,2,2),'Color',[0 1 0]); %post
       %         yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)]); xlim([0.5 8.5])

    end
end


%%% plot speed histogram
figure
hold on
subplot(2,2,1)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==saline,1:25,:),1))); title('saline'); xlabel('speed')
subplot(2,2,2)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==doi,1:25,:),1))); title('doi'); xlabel('speed')
subplot(2,2,3)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==ketanserin,1:25,:),1))); title('ketanserin'); xlabel('speed')
subplot(2,2,4)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==ketandoi,1:25,:),1))); title('ketanserin + DOI'); xlabel('speed')
legend('pre','post')


%use= find(wn_evoked(:,2,1)>.3 & wn_evoked(:,2,2)>.3)
%use= find((wn_evoked(:,2,1)>.3 & wn_evoked(:,2,2)>.3) | wn_evoked(:,2,:)> 0)

% right now there's an error where units between treatments seem to be
% mixed...running all saline alone gives a different MI distribution than
% when all treatments are run together
titles = {'Saline','DOI','Ketanserin', 'Ketanserin + DOI'};
figure
for t = 1:4
    useDark = meanRdark(:,1)>.3 | meanRdark(:,2)>.3;
    useN = find(treatment==t)
    for i = 1:length(useN)
        MIdark= (meanRdark(:,2)-meanRdark(:,1))./(meanRdark(:,2)+meanRdark(:,1));
    end
    %subplot(3,2,1)
    %plot(SI_mv_wn(useFR(i) & treatment==t),'.'); ylim([0 1.5]);hold on
    %legend ('Saline', 'DOI', 'Ketanserin','Ketanserin + DOI')
    subplot(3,2,t+2)
    h= hist(MIdark(useDark(treatment==t)),-1:.1:1);
    Mbins=-1:.1:1
    bar(Mbins,h/sum(useDark(treatment==t)))
    xlim([-1.5 1.5]); ylim([0 .25]);axis xy
    xlabel('MI'); ylabel('fraction of cells');title(titles{t});
    set(gcf,'Name','MI Dark')
    %hold on
    subplot(3,2,2)
    meanMIdark(t) = nanmean(MIdark(treatment==t))
    %err(t) = nanstd(MIdark(treatment==t))/sqrt(sum(MIdark (treatment==t)));
    bar(meanMIdark)
    %barweb(meanMIdark(1:1:t)',err(1:1:t)');
    ylim([-0.2 0.2]);
end

%% mod for corr in darkness %%

titles = {'Saline','DOI','Ketanserin', 'Ketanserin + DOI'};
figure
for t = 1:4
    useN = find(treatment==t)
    for i = 1:length(useN)
        MIdarkCorr= (darkCorr(:,2)-darkCorr(:,1))./(darkCorr(:,2)+darkCorr(:,1));
    end
    subplot(2,2,t)
    h= hist(MIdarkCorr(useN(i) &treatment==t),-1:.1:1);
    Mbins=-1:.1:1
    bar(Mbins,h/sum(useN(i) & treatment==t))
    xlim([-1.5 1.5]); ylim([0 .25]);axis xy
    xlabel('Dark pairwise MI'); ylabel('fraction of cells');title(titles{t});
    set(gcf,'Name','Corr MI Dark')
end

% mod for WN corr
titles = {'Saline','DOI','Ketanserin', 'Ketanserin + DOI'};
figure
for t = 1:4
    useN = find(treatment==t)
    for i = 1:length(useN)
        MIwnCorr= (wnCorr(:,2)-wnCorr(:,1))./(wnCorr(:,2)+wnCorr(:,1));
    end
    subplot(2,2,t)
    h= hist(MIwnCorr(useN(i) &treatment==t),-1:.1:1);
    Mbins=-1:.1:1
    bar(Mbins,h/sum(useN(i) & treatment==t))
    xlim([-1.5 1.5]); ylim([0 .35]);axis xy
    xlabel('wn pairwise MI'); ylabel('fraction of cells');title(titles{t});
    set(gcf,'Name','Corr MI Wn')
end


%ALL trials (not separated by moving/stationary)
figure
for t = 1:4
    thresh = 1;
    useEv =wn_evoked(:,:,1)>0 & wn_evoked(:,:,2)>0 & (wn_evoked(:,:,1)>thresh | wn_evoked(:,:,2)>thresh);
    useN = find(treatment==t)
    for i = 1:length(useN)
            MI_wn = (wn_evoked(:,:,2)-wn_evoked(:,:,1))./(wn_evoked(:,:,2)+wn_evoked(:,:,1));
    end
        subplot(3,2,t+2)
        hEv= hist(MI_wn(useEv(treatment==t)),-1:.1:1);
        Mbins=-1:.1:1
        bar(Mbins,hEv/sum(useEv(treatment==t)))
        xlabel('MI'); ylabel('fraction of cells'); xlim([-1.5 1.5]);ylim([0 .3]);
        set(gcf,'Name','MI wn all')
        title(titles{t});
        %subplot(3,2,2)
       % meanMI_mv_wn(t) = nanmean(MI_mv_wn(treatment==t))
       % bar(meanMI_mv_wn);ylim([-1 1]);
end

% MI MOVING!
figure
for t = 1:4
    thresh = 1;
    useEv =wn_evoked(:,2,1)>0 & wn_evoked(:,2,2)>0 & (wn_evoked(:,2,1)>thresh | wn_evoked(:,2,2)>thresh);
    useN = find(treatment==t)
    for i = 1:length(useN)
            MI_mv_wn = (wn_evoked(:,2,2)-wn_evoked(:,2,1))./(wn_evoked(:,2,2)+wn_evoked(:,2,1));
    end
        
        subplot(3,2,t+2)
        hEv= hist(MI_mv_wn(useEv(treatment==t)),-1:.1:1);
        Mbins=-1:.1:1
        bar(Mbins,hEv/sum(useEv(treatment==t)))
        xlabel('MI'); ylabel('fraction of cells'); xlim([-1.5 1.5]);ylim([0 .3]);
        set(gcf,'Name','MI move wn')
        title(titles{t});
        subplot(3,2,2)
       meanMI_mv_wn(t) = nanmean(MI_mv_wn(treatment==t))
       bar(meanMI_mv_wn);ylim([-1 1]);
end

%MI wn stationary:
figure
for t = 1:4
    thresh = 1;
    useEv =wn_evoked(:,1,1)>0 & wn_evoked(:,1,2)>0 & (wn_evoked(:,1,1)>thresh | wn_evoked(:,1,2)>thresh);
    useN = find(treatment==t)
    for i = 1:length(useN)
            MI_stat_wn = (wn_evoked(:,1,2)-wn_evoked(:,1,1))./(wn_evoked(:,1,2)+wn_evoked(:,1,1));
    end
        subplot(3,2,t+2)
        hEv_stat= hist(MI_stat_wn(useEv(treatment==t)),-1:.1:1);
        Mbins=-1:.1:1
        bar(Mbins,hEv_stat/sum(useEv(treatment==t)))
        xlabel('MI'); ylabel('fraction of cells'); xlim([-1.5 1.5]);ylim([0 .3]);
        set(gcf,'Name','MI stationary wn')
        title(titles{t});
        subplot(3,2,2)
       meanMI_stat_wn(t) = nanmean(MI_stat_wn(treatment==t))
       bar(meanMI_stat_wn);ylim([-1 1]);
end

% MI wn SPONTANEOUS all 
figure
for t = 1:4
    thresh = 1;
    useEv =wn_spont(:,:,1)>0 & wn_spont(:,:,2)>0 & (wn_spont(:,:,1)>thresh | wn_spont(:,:,2)>thresh);
    useN = find(treatment==t)
    for i = 1:length(useN)
            MI_wn_spont = (wn_spont(:,:,2)-wn_spont(:,:,1))./(wn_spont(:,:,2)+wn_spont(:,:,1));
    end
        subplot(3,2,t+2)
        hEv= hist(MI_wn_spont(useEv(treatment==t)),-1:.1:1);
        Mbins=-1:.1:1
        bar(Mbins,hEv/sum(useEv(treatment==t)))
        xlabel('MI'); ylabel('fraction of cells'); xlim([-1.5 1.5]);ylim([0 .4]);
        set(gcf,'Name','MI wn spontaneous all')
        title(titles{t});
        subplot(3,2,2)
        meanMI_wn_spont(t) = nanmean(MI_wn_spont(treatment==t))
        bar(meanMI_wn_spont);ylim([-1 1]);
end

% MI wn SPONTANEOUS MOVING
figure
for t = 1:4
    thresh = 1;
    useEv =wn_spont(:,2,1)>0 & wn_spont(:,2,2)>0 & (wn_spont(:,2,1)>thresh | wn_spont(:,2,2)>thresh);
    useN = find(treatment==t)
    for i = 1:length(useN)
            MI_mv_spont = (wn_spont(:,2,2)-wn_spont(:,2,1))./(wn_spont(:,2,2)+wn_spont(:,2,1));
    end
        subplot(3,2,t+2)
        hSpont= hist(MI_mv_spont(useEv(treatment==t)),-1:.1:1);
        Mbins=-1:.1:1
        bar(Mbins,hSpont/sum(useEv(treatment==t)))
        xlabel('MI'); ylabel('fraction of cells'); xlim([-1.5 1.5]);ylim([0 .4]);
        set(gcf,'Name','MI wn spontaneous move')
        title(titles{t});
        subplot(3,2,2)
       meanMI_mv_spont(t) = nanmean(MI_mv_spont(treatment==t))
       bar(meanMI_mv_spont);ylim([-1 1]);
end

figure
for t = 1:4
    thresh = 1;
    useEv =wn_spont(:,1,1)>0 & wn_spont(:,1,2)>0 & (wn_spont(:,1,1)>thresh | wn_spont(:,1,2)>thresh);
    useN = find(treatment==t)
    for i = 1:length(useN)
            MI_stat_spont = (wn_spont(:,1,2)-wn_spont(:,1,1))./(wn_spont(:,1,2)+wn_spont(:,1,1));
    end
        subplot(3,2,t+2)
        hSpont= hist(MI_stat_spont(useEv(treatment==t)),-1:.1:1);
        Mbins=-1:.1:1
        bar(Mbins,hSpont/sum(useEv(treatment==t)))
        xlabel('MI'); ylabel('fraction of cells'); xlim([-1.5 1.5]);ylim([0 .4]);
        set(gcf,'Name','MI wn Spontaneous stationary')
        title(titles{t});
        subplot(3,2,2)
       meanMI_stat_spont(t) = nanmean(MI_stat_spont(treatment==t))
       bar(meanMI_stat_spont);ylim([-1 1]);
end



titles = {'Saline','DOI','Ketanserin', 'Ketanserin + DOI'};
figure
for t=1:4
    useN = find(treatment==t)
    for i = 1:length(useN)
        subplot(2,2,t)
        plot(cv2Dark(treatment==t,1), cv2Dark(treatment==t,2),'.');title(titles{t});
        xlabel('Pre CV2 dark');ylabel('Post CV2 dark'); ylim([0 2]);xlim([0 2]); hold on;
        plot([0 2],[0 2])
    end
end

%%%%%%%% ===================================================== %%%%%%%%%%
%%% try to decode SF, compare pre and post sessions for each treatment %%%

% for t=1:4
% useN= find(treatment==t)
% for i = 1:length(useN)
% lowSFpre = find(drift_sf(:,1,:,1))
% highSFpre = find(drift_sf(:,7,:,1))
% 
% end
% end















