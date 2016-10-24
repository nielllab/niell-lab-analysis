clear all
close all 

dbstop if error

batchDOIephys_filtered; %%% load batch file
%set(groot,'defaultFigureVisible','off') %disable figure plotting
set(groot,'defaultFigureVisible','on')

%%% select the sessions you want based on filters
%use =  find(strcmp({files.notes},'good data')& ~cellfun(@isempty,{files.predark})& ~cellfun(@isempty,{files.postdark}) )
%useSess = use;
%use =  find( strcmp({files.treatment},'5HT') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}) )

%for specific experiment:
use =  find(strcmp({files.notes},'good data')  & ~cellfun(@isempty,{files.predark})& ~cellfun(@isempty,{files.postdark}) & strcmp({files.expt},'102116'))
sprintf('%d selected sessions',length(use))

saline=1; doi=2; ht=3; ketanserin=4; ketandoi=5; mglur2=6; mglur2doi=7; lisuride=8;
 
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
    
    good = ones(1,nc);
    good(files(use(i)).badsites) =0;
    goodAll(cellrange)=good;
    
%     
%     if ~isempty(files(use(i)).prepinpFile)
%         pinp =getPinpHTR2A(clustfile,afile, files(use(i)).prepinpFile,3,5,0,0,1)
%         pinpAll(cellrange,1) = pinp.pinped;
% %         pinpPsth(cellrange,:,1)=pinp.psth;
%     else
%         pinpAll(cellrange,1)=NaN;
% %         pinpPsth(cellrange,:,1)=NaN;
%     end
%     
%     if ~isempty(files(use(i)).postpinpFile)
%         pinp =getPinpHTR2A(clustfile,afile, files(use(i)).postpinpFile,3,5,0,0,1)
%         pinpAll(cellrange,2) = pinp.pinped;
% %         pinpPsth(cellrange,:,2)=pinp.psth;
%     else
%         pinpAll(cellrange,2)=NaN;
% %        pinpPsth(cellrange,:,2)=NaN;
%     end
%     
    
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
%     
% if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2}) 
%     
% for prepost =1:2
% eyes = getEyes(clustfile,afile,cfile{:,prepost}, files(use(i)).blockWn{prepost},1);
% rad{:,prepost} = eyes.rad
% t{prepost} = eyes.t
% end
% end

    %%% session info
    sessionNum(cellrange)=i;
    for j=1:length(cellrange); expt{cellrange(j)} = files(use(i)).expt; end
    if strcmp(files(use(i)).treatment,'Saline'), treatment(cellrange)=saline, end;
    if strcmp(files(use(i)).treatment,'DOI'), treatment(cellrange)=doi, end;
    if strcmp(files(use(i)).treatment,'5HT'), treatment(cellrange)=ht, end;
    if strcmp(files(use(i)).treatment,'Ketanserin'), treatment(cellrange)=ketanserin, end;
    if strcmp(files(use(i)).treatment,'KetanserinDOI'), treatment(cellrange)=ketandoi, end;
    if strcmp(files(use(i)).treatment,'MGluR2'), treatment(cellrange)=mglur2, end;
    if strcmp(files(use(i)).treatment,'MGluR2DOI'), treatment(cellrange)=mglur2doi, end;
    if strcmp(files(use(i)).treatment,'Lisuride'), treatment(cellrange)=lisuride, end;

    sessionTreatment(i) = treatment(cellrange(1));
    
    if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2})
        % get pre/post running speed
        for prepost = 1:2
            spd = getSpeed(clustfile,afile,files(use(i)).blockWn{prepost},1);
            speedHistWn(i,:,prepost) = hist(spd.v,0.5:1:100)/length(spd.v);
            speedTrace{i,prepost}=spd.v;
        end
    else
        for prepost = 1:2
                 spd = getSpeed(clustfile,afile,files(use(i)).blockDrift{prepost},1);
            speedHistDrift(i,:,prepost) = hist(spd.v,0.5:1:100)/length(spd.v);
            speedTrace{i,prepost}=spd.v;
        end
           
    end

    figure
    for prepost = 1:2
        subplot(2,1,prepost);
        plot(speedTrace{i,prepost}); hold on
        plot([1 length(speedTrace{i,prepost})] , [ 0.5 0.5], ':'); ylabel('speed'); ylim([0 10])
        if prepost==2,  title(sprintf('%s post %s',files(use(i)).expt, files(use(i)).treatment)); end
    end
    
      if ~isempty(files(use(i)).blockDrift{1}) & ~isempty(files(use(i)).blockDrift{2})
          %%% get grating responses
          for prepost = 1:2
              drift = getDrift_mv(clustfile,afile,files(use(i)).blockDrift{prepost},1);
              drift_orient(cellrange,:,:,prepost)=drift.orient_tune;
%               drift_sf(cellrange,:,:,prepost) = drift.sf_tune;
              drift_spont(cellrange,:,prepost) = drift.interSpont;
              drift_osi(cellrange,:,prepost) = drift.cv_osi;
              drift_F1F0(:,cellrange,prepost)= drift.F1F0;
              drift_ot_tune(cellrange,:,:,prepost)=drift.orient_tune;
              drift_sf_trial{prepost} = drift.trialSF;
              drift_orient_trial{prepost} = drift.trialOrient;
              for c = 1:length(cellrange)
                  drift_trial_psth{cellrange(c),:,prepost} = squeeze(drift.trialPsth(c,:,:));
                  mnPsth(cellrange(c),prepost,:) = squeeze(mean(drift.trialPsth(c,:,:),2));
              end
          
      
              for cond = 1:72
                  for mv = 1:2
                      if mv ==1
                          tr = find(drift.trialOrient(1:end-1) == ceil(cond/6) & drift.trialSF(1:end-1)==mod(cond-1,6)+1 & drift.frameSpd<1);
                      else
                          tr = find(drift.trialOrient(1:end-1) == ceil(cond/6) & drift.trialSF(1:end-1)==mod(cond-1,6)+1 & drift.frameSpd>1);
                      end
                      drift_tcourse(cellrange,cond,mv,prepost,:) = squeeze(nanmean(drift.trialPsth(:,tr,:),2));
                      sf = mod(cond-1,6)+1; ori = ceil(cond/6);
                      drift_cond_tcourse(cellrange,mv,prepost,ori,sf,:) = squeeze(nanmean(drift.trialPsth(:,tr,:),2));
                  end
              end
              
          end
      end
    
        if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2})

%     %%% get wn response
    for prepost = 1:2
        wn = getWn_mv(clustfile,afile,files(use(i)).blockWn{prepost},1,300);
        wn_crf(cellrange,:,:,prepost)=wn.crf;
        wn_spont(cellrange,:,prepost)=wn.spont;
        wn_evoked(cellrange,:,prepost)=wn.evoked;
        wn_gain(1,cellrange,prepost)=wn.gain;
        wn_frameR(cellrange,:,:,prepost) = downsamplebin(wn.frameR,2,15)/15;
        wn_reliability(1,cellrange,prepost) = wn.reliability;

    end
        end
        
%         figure
%         subplot(2,2,1); imagesc(squeeze(wn_frameR(cellrange,:,1,1)),[ 0 50]); title('pre stop')
%         subplot(2,2,2); imagesc(squeeze(wn_frameR(cellrange,:,1,2)),[ 0 50]); title('post stop')
%         subplot(2,2,3); imagesc(squeeze(wn_frameR(cellrange,:,2,1)),[ 0 50]); title('pre move');
%         subplot(2,2,4); imagesc(squeeze(wn_frameR(cellrange,:,2,2)),[ 0 50]); title('post move')
%         set(gcf,'Name',sprintf('%s %s',files(use(i)).expt, files(use(i)).treatment));
%         drawnow
        
    %%% lfp power
%     %%%(right now averages over all sites, should use layer info)
%     for prepost=1:2
%         lfpMove = getLfpMovement(clustfile,afile,files(use(i)).blockWn{prepost},0);
%         LFPall(i,:,:,prepost) =squeeze(nanmedian(lfpMove.meanSpect, 1))/nanmedian(lfpMove.meanSpect(:));
%     end
%     
%       for prepost=1:2
%         lfpMoveDark = getLfpMovement(clustfile,afile,files(use(i)).blockDark{prepost},0);
%         LFPallDark(i,:,:,prepost) =squeeze(nanmedian(lfpMove.meanSpect, 1))/nanmedian(lfpMove.meanSpect(:));
%    
%       end
    

  getSorting(clustfile,afile,sprintf('%s %s',files(use(i)).expt,files(use(i)).treatment));
  drawnow
  

    %% darkness / correlation analysis
    
    %%% need to keep track of n^2 values for correlations
    corrRange=ncorr+1:ncorr+nc^2;
    corrTreatment(corrRange)=treatment(cellrange(1));
    
    %%% prepost correlation for white noise
    if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2})
        dt = 1;
        [preCorr postCorr cv2 R eigs] =  prepostDOIdarkness(clustfile,afile,files(use(i)).blockWn,dt,0);
        wnCorr(corrRange,1) = preCorr(:); wnCorr(corrRange,2)=postCorr(:);
    end
    if ~isempty(files(use(i)).blockDrift{1}) & ~isempty(files(use(i)).blockDrift{2})
        dt = 1;
        [preCorr postCorr cv2 R eigs] =  prepostDOIdarkness(clustfile,afile,files(use(i)).blockDrift,dt,0);
        driftCorr(corrRange,1) = preCorr(:); driftCorr(corrRange,2)=postCorr(:);
    end
    
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

%  save('compileALL_101416')

% d= squeeze(drift_tcourse(goodAll==1 & treatment==doi,:,1,:,1));  %%% stationary first timepoint
% d = downsamplebin(d,2,2,1);
% d = reshape(d,size(d,1),size(d,2)*size(d,3));
% nanratio = mean(isnan(d),2);
% figure
% plot(nanratio)
% 
% for prepost = 1:2
%     
%     tcourse = squeeze(drift_tcourse(goodAll==1 & treatment==doi,:,1,prepost,:));
%     tcourse = downsamplebin(tcourse(:,:,2:49),3,4,1);
%     tcourse = downsamplebin(tcourse,2,2,1);
%     tcourse = permute(tcourse,[1 3 2]);
%     trace = reshape(tcourse,size(tcourse,1),size(tcourse,2)*size(tcourse,3));
%     
%     figure
%     imagesc(isnan(squeeze(trace)))
%     
%     if prepost ==1
%         s = nanstd(trace,[],2);
%     end
%     figure
%     hist(s,0.25:0.5:15)
%     usespks = find(s>1 & nanratio<0.04);
%     figure
%     imagesc(isnan(squeeze(trace(usespks,:))));
%  
%     mn = nanmean(trace,2);
%     
%     normtrace = (trace - repmat(mn,[1 size(trace,2)]))./repmat(s,[1 size(trace,2)]);
%     figure
%     imagesc(normtrace(usespks,:),[ -5 5])
%     
%     if prepost==1
%         [coeff score latent] = pca(normtrace(usespks,:)');  score(:,12) = mean(normtrace(usespks,:),1);
%         %[coeff score latent] = pca(normtrace(usespks,:)');
%         
%         figure
%         plot(latent(1:20)/sum(latent));
%         
%     end
%        data = normtrace(usespks,:)';
%        score = data*coeff;
%     
%     figure
%     for i = 1:5;
%         subplot(5,1,i);
%         plot(score(:,i));
%     end
%     
%     condtuning = reshape(score',size(score,2),size(tcourse,2),size(tcourse,3));
%     tuning = reshape(condtuning,size(condtuning,1),size(condtuning,2),3,12);
%     
%     figure
%     for i = 1:12
%         subplot(3,4,i)
%         d=circshift(squeeze(nanmean(tuning(i,:,:,:),4)),2); %%% average over either sf or orientation
%         plot(d);  range = max(abs(d(:)))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 size(d,1)+0.5])
%     end
%     
%     figure
%     for i = 1:12
%         subplot(3,4,i)
%         d=circshift(squeeze(mean(tuning(i,:,:,:),3)),2);
%         plot(d);  range = max(abs(d(:)))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 size(d,1)+0.5])
%     end
%     
%     figure
%     for i = 1:12
%         subplot(3,4,i)
%         d=squeeze(mean(mean(tuning(i,2:7,:,:),3),2))- squeeze(mean(mean(tuning(i,10:12,:,:),3),2));
%         plot(d); range = max(abs(d))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 12.5])
%     end
%     
%     figure
%     for i = 1:12
%         subplot(3,4,i)
%         d=squeeze(mean(tuning(i,2:7,:,:),2))- squeeze(mean(tuning(i,10:12,:,:),2));
%         plot(d'); range = max(abs(d(:)))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 12.5])
%     end
%     
%     figure
%     for i = 1:12
%         subplot(3,4,i)
%         d=squeeze(mean(mean(tuning(i,1,:,:),3),2))- squeeze(mean(mean(tuning(i,10:12,:,:),3),2));
%         plot(d); range = max(abs(d))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 12.5])
%     end
%     
% end
% 
% tcourse = squeeze(drift_tcourse(:,:,1,:,:));
% tcourse = downsamplebin(tcourse(:,:,:,2:49),4,4,1);
% tcourse = downsamplebin(tcourse,2,2,1);
% tcourse = permute(tcourse,[1 4 2 3]);
% trace = reshape(tcourse,size(tcourse,1),size(tcourse,2)*size(tcourse,3)*size(tcourse,4));
% 
% figure
% imagesc(isnan(squeeze(trace)))
% 
% if prepost ==1
%     s = nanstd(trace,[],2);
% end
% figure
% hist(s,0.25:0.5:15)
% usespks = find(s>1 & nanratio<0.04);
% figure
% imagesc(isnan(squeeze(trace(usespks,:))));
% 
% mn = nanmean(trace,2);
% 
% normtrace = (trace - repmat(mn,[1 size(trace,2)]))./repmat(s,[1 size(trace,2)]);
% figure
% imagesc(normtrace(usespks,:),[ -5 5])
% 
% [coeff scoreAll latent] = pca(normtrace(usespks,:)');
% %[coeff score latent] = pca(normtrace(usespks,:)');
% 
% 
% figure
% for i = 1:5;
%     subplot(5,1,i);
%     plot(scoreAll(:,i));
% end
% 
% for prepost = 1:2
%     if prepost==1
%         score = scoreAll(1:size(scoreAll,1)/2,:);
%     else
%         score = scoreAll((size(scoreAll,1)/2+1):end,:);
%     end
%     
%     condtuning = reshape(score',size(score,2),size(tcourse,2),size(tcourse,3));
%     tuning = reshape(condtuning,size(condtuning,1),size(condtuning,2),3,12);
%     
%     figure
%     for i = 1:12
%         subplot(3,4,i)
%         d=circshift(squeeze(nanmean(tuning(i,:,:,:),4)),2);
%         plot(d);  range = max(abs(d(:)))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 size(d,1)+0.5])
%     end
%     
%     figure
%     for i = 1:12
%         subplot(3,4,i)
%         d=circshift(squeeze(mean(tuning(i,:,:,:),3)),2);
%         plot(d);  range = max(abs(d(:)))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 size(d,1)+0.5])
%     end
%     
%     figure
%     for i = 1:12
%         subplot(3,4,i)
%         d=squeeze(mean(mean(tuning(i,2:7,:,:),3),2))- squeeze(mean(mean(tuning(i,10:12,:,:),3),2));
%         plot(d); range = max(abs(d))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 12.5])
%     end
%     
%     figure
%     for i = 1:12
%         subplot(3,4,i)
%         d=squeeze(mean(tuning(i,2:7,:,:),2))- squeeze(mean(tuning(i,10:12,:,:),2));
%         plot(d'); range = max(abs(d(:)))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 12.5])
%     end
%     
%     figure
%     for i = 1:12
%         subplot(3,4,i)
%         d=squeeze(mean(mean(tuning(i,1,:,:),3),2))- squeeze(mean(mean(tuning(i,10:12,:,:),3),2));
%         plot(d); range = max(abs(d))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 12.5])
%     end
% end


amp_wn = squeeze(max(wn_crf,[],2));
useResp = amp_wn(:,1,1)>1 & amp_wn(:,1,2)>1;
data_wn = goodAll==1 & useResp' & ~inhAll

%ketanserin & ketan doi  more stable than saline????

for t=1:3
    figure
   if t == 1
set(gcf,'Name', 'mean SALINE CRF');
   elseif t==2
       set(gcf,'Name', 'mean DOI CRF');
   else
       set(gcf,'Name', 'mean 5HT CRF');
   end
for prepost=1:2
    for l =1:6
        subplot(2,3,l)
        mn_wn = squeeze(mean(wn_crf(data_wn & layerAll==l & treatment==t,:,prepost))); %stationary pre
        plot(mn_wn);hold on
    end
end
end

useN= data_wn & treatment==doi
MI_crf = (wn_crf(useN,:,1,2)-wn_crf(useN,:,1,1))./(wn_crf(useN,:,1,2)+wn_crf(useN,:,1,1));

CRFsuppress = find(MI_crf <-.4)
CRFup = find(MI_crf >.4)

h= hist(MI_crf,-1:.1:1);
figure
subplot(2,2,1)
Mbins=-1:.1:1
bar(Mbins,h/sum(useN))

h= hist(MI_crf(CRFsuppress),-1:.1:1)
subplot(2,2,3)
Mbins=-1:.1:1
bar(Mbins,h/sum(useN));ylim([0 1]);xlim([-1.5 1.5]); title('suppressed');

h= hist(MI_crf(CRFup),-1:.1:1);
subplot(2,2,4)
Mbins=-1:.1:1
bar(Mbins,h/sum(useN));ylim([0 1]);xlim([-1.5 1.5]);title('facilitated');

figure
set(gcf, 'Name','facilitated units CRF')
for c=1:ceil(length(CRFup))
    %np = ceil(sqrt(length(CRFup)));
    %subplot(np,np,c);
    subplot(10,10,c);
    plot(wn_crf(CRFup(c),:,1,1),'Color',[0.5 0 0]); hold on;  plot(wn_crf(CRFup(c),:,2,1),'Color',[0 0.5 0]);
    plot(wn_crf(CRFup(c),:,1,2),'Color',[1 0 0]);  plot(wn_crf(CRFup(c),:,2,2),'Color',[0 1 0]);
    yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)]);
end

figure
set(gcf, 'Name','suppressed units CRF')
for c=1:ceil(length(CRFsuppress))
    np = ceil(sqrt(length(CRFsuppress)));
    subplot(np,np,c);
    % plot(wn_crf(useN(c),:,prepost))
    plot(wn_crf(CRFsuppress(c),:,1,1),'Color',[0.5 0 0]); hold on;  plot(wn_crf(CRFsuppress(c),:,2,1),'Color',[0 0.5 0]);
    plot(wn_crf(CRFsuppress(c),:,1,2),'Color',[1 0 0]);  plot(wn_crf(CRFsuppress(c),:,2,2),'Color',[0 1 0]);
    yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)]);
end

% figure
% plot(wn_reliability(1,data & treatment==doi,1)); hold on;
% plot(wn_reliability(1,data & treatment==doi,2))

amp = squeeze(max(drift_orient,[],2));
useResp = amp(:,1,1)>2 | amp(:,1,2)>2;
data = goodAll & useResp' & ~inhAll %which cells to include
dataInh =goodAll & useResp' & inhAll 
useN=data & treatment==ht

% redo this section with drift psth then look at OSI prepost
%miDrift= (drift_cond_tcourse(useN,1,2,:,:,:)-drift_cond_tcourse(useN,1,1,:,:,:))./(drift_cond_tcourse(useN,1,2,:,:,:)+drift_cond_tcourse(useN,1,1,:,:,:));
%miDrift= (drift_tcourse(useN,:,1,2,:)-drift_tcourse(useN,:,1,1,:)./(drift_tcourse(useN,:,1,2,:))+drift_tcourse(useN,:,1,1,:));
miDrift= (mnPsth(useN,2,:)-mnPsth(useN,1,:))./((mnPsth(useN,2,:)+mnPsth(useN,1,:)));
% miDrift= (drift_trial_psth(useN,:,2)-drift_trial_psth(useN,:,1))./(drift_trial_psth(useN,:,2)+drift_trial_psth(useN,:,1));

h= hist(miDrift,-1:.1:1);
figure
subplot(2,2,1)
Mbins=-1:.1:1
bar(Mbins,h/sum(useN));ylim([0 .35]);

driftSuppress = miDrift <-.4;
driftUp = miDrift >.4;

h= hist(miDrift(driftSuppress),-1:.1:1)
subplot(2,2,3)
Mbins=-1:.1:1
bar(Mbins,h/sum(useN));ylim([0 4]);xlim([-1.5 1.5]); title('suppressed');

h= hist(miDrift(driftUp),-1:.1:1);
subplot(2,2,4)
Mbins=-1:.1:1
bar(Mbins,h/sum(useN));ylim([0 4]);xlim([-1.5 1.5]);title('facilitated');

% figure
% for c=1:length(useN)
%     if driftSuppress ==1
%         osiHist=hist(drift_osi(useN(c),:,1))
%         bar(Mbins,osiHist/sum(useN));hold on
%     else
%         display 'not suppressed';
%     end
% end

figure
subplot(1,2,1)
imagesc(mnPsth(data,1)); %%% show full dataset (n x t)
subplot(1,2,2)
imagesc(mnPsth(data,2))

%drift_psth= cellfun(@mean, drift_trial_psth(data,1,:),'UniformOutput', false)
dt = 0.05;
figure
set(gcf,'Name', 'prepost response to 6 SFs')
for sf = 1:6
    subplot(2,3,sf)
    mn_pre=(squeeze(mean(drift_cond_tcourse(data & treatment==ht,1,1,3,sf,:)))); % cellrange,mv,prepost,ori,sf,:
    mn_pre = mn_pre - repmat(mn_pre(1,:),[50 1]);
    mn_pre = circshift(mn_pre,10)
    mn_post =(squeeze(mean(drift_cond_tcourse(data & treatment==ht,1,2,3,sf,:))));
    mn_post = mn_post - repmat(mn_post(1,:),[50 1]);
    mn_post = circshift(mn_post,10)
    plot((1:length(mn_pre)-5)*dt -dt/2,mn_pre(1:45,:));xlabel('secs'); ylabel('sp/sec');
    hold on;
    plot((1:length(mn_post)-5)*dt -dt/2,mn_post(1:45,:));
end


dt = 0.05;
figure
set(gcf,'Name', 'prepost response to 6 orientations')
for ori = 1:12
    subplot(3,4,ori)
    mn_pre=(squeeze(mean(drift_cond_tcourse(data & treatment==ht,1,1,ori,3,:)))); % cellrange,mv,prepost,ori,sf,:
    mn_pre = mn_pre - repmat(mn_pre(1,:),[50 1]);
    mn_pre = circshift(mn_pre,10)
    mn_post =(squeeze(mean(drift_cond_tcourse(data & treatment==ht,1,2,ori,3,:))));
    mn_post = mn_post - repmat(mn_post(1,:),[50 1]);
    mn_post = circshift(mn_post,10)
    plot((1:length(mn_pre)-5)*dt -dt/2,mn_pre(1:45,:));xlabel('secs'); ylabel('sp/sec');
    hold on;
    plot((1:length(mn_post)-5)*dt -dt/2,mn_post(1:45,:));
end

figure
plot(drift_ot_tune(data,:,1,1),'b'); hold on;
plot(drift_ot_tune(data,:,1,2),'r');

data = reshape(data,size(data,1),size(data,2))
%for t=1:7
%dist = pdist(squeeze(mnPsth(data,1,:)),'correlation');  %%% sort based on correlation coefficient...want to do this for just one orientation or SF
dist = pdist(squeeze(drift_cond_tcourse(data & treatment==doi,1,1,3,3,:)),'correlation'); 
display('doing cluster')
tic, Z = linkage(dist,'ward'); toc
figure
subplot(3,4,[1 5 9 ])
display('doing dendrogram')
[h t perm] = dendrogram(Z,0,'Orientation','Left','ColorThreshold' ,5);
axis off
subplot(3,4,[2 3 4 6 7 8 10 11 12 ]);
imagesc(drift_cond_tcourse(perm,:)); axis xy   %%% show sorted data
%end

[Y e] = mdscale(dist,1);
[y sortind] = sort(Y);
figure
imagesc(mnPsth(sortind,:),[0 1])


%%% plot mean timecourse across all drift stim and cells, by layer / cell-type
%%% this should be updated to select active cells and separate
%%% moving/stationary. Also maybe choose optimal stim for each cell?


%data2 = reshape(data,size(data,1),size(data,2),size(data,3))
titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
dt = 0.05;
for i=1:4
    figure 
    for t=1:7
        if i==1
            set(gcf,'Name','grating lyr5');
            mn = squeeze(mean(mnPsth(data & layerAll==5 & treatment==t,:,:),1))';

        elseif i==2
            set(gcf,'Name','grating lyr2-4');
            mn = squeeze(mean(mnPsth(data & layerAll<5 & treatment==t,:,:),1))'; 
        elseif i==3
            set(gcf,'Name','grating inh');
            mn = squeeze(mean(mnPsth(dataInh & treatment==t,:,:),1))';
%         elseif i==4
%             set(gcf,'Name','grating lyr6');
%             mn = squeeze(mean(mnPsth(data & layerAll==6 & treatment==t,:,:),1'));
        end
        mn = mn - repmat(mn(1,:),[50 1]);
        mn = circshift(mn,10)
        subplot(2,4,t)
        plot((1:length(mn)-5)*dt -dt/2,mn(1:45,:)); title(titles{t}); xlabel('secs'); ylabel('sp/sec');
        %ylim([0 max(mn(:))+1])
    end
end

% 
% amp = squeeze(max(drift_orient,[],2));
% useResp = amp(:,:,1)>2 | amp(:,:,2)>2;
% figure
% useN= find(goodAll==1 & useResp' & layerAll==2 & ~inhAll)

% for t=1:8
%     figure 
%     set (gcf,'Name','Lyr 2 psth/unit')
% for c=1:length(useN)
%     subplot(11,11,c)
%     plot(squeeze(mnPsth(useN(c) & treatment==t,1,:)));hold on;
%     plot(squeeze(mnPsth(useN(c) &treatment==t,2,:)),'r')
% end
% end

% % amp = squeeze(max(drift_orient,[],2));
% % useResp = amp(:,:,1)>2 | amp(:,:,2)>2;
figure
useN= find(data & layerAll==3 |layerAll==2 &treatment==doi)
for c=1:length(useN)
    set (gcf,'Name','Lyr 2/3 psth/unit')
    subplot(8,10,c)
    plot(squeeze(mnPsth(useN(c),1,:)));hold on;
    plot(squeeze(mnPsth(useN(c),2,:)),'r')
end


figure
useN= find(data & layerAll==4 & treatment==doi)
for c=1:length(useN)
    set (gcf,'Name','Lyr 4 psth/unit')
    subplot(5,10,c)
    plot(squeeze(mnPsth(useN(c),1,:)));hold on;
    plot(squeeze(mnPsth(useN(c),2,:)),'r')
end


figure
for c=1:length(useN)
    set(gcf,'Name','hist of change prepost lyr4 exc psth')
    MImnPsth= (mnPsth(useN,2,:)-mnPsth(useN,1,:))./(mnPsth(useN,2,:)+mnPsth(useN,1,:));
    hist(MImnPsth)
    %subplot(4,2,t)
    h= hist(MImnPsth(useN),-1:.1:1);
    Mbins=-1:.1:1
    bar(Mbins,h/sum(useN))
end

figure
useN= find(data & layerAll==5 & treatment==doi)
for c=1:length(useN)
    set (gcf,'Name','Lyr 5 psth/unit')
    subplot(5,10,c)
    plot(squeeze(mnPsth(useN(c),1,:)));hold on;
    plot(squeeze(mnPsth(useN(c),2,:)),'r')
end


figure
for c=1:length(useN)
    set(gcf,'Name','hist of change prepost lyr5 exc psth')
    MImnPsth= (mnPsth(useN,2,:)-mnPsth(useN,1,:))./(mnPsth(useN,2,:)+mnPsth(useN,1,:));
    hist(MImnPsth)
    %subplot(4,2,t)
    h= hist(MImnPsth(useN),-1:.1:1);
    Mbins=-1:.1:1
    bar(Mbins,h/sum(useN))
end


figure
useInh= find(dataInh)
for c=1:length(useInh)
    set (gcf,'Name','inh psth/unit')
    subplot(5,6,c)
    plot(squeeze(mnPsth(useInh(c),1,:)));axis square;hold on;
    plot(squeeze(mnPsth(useInh(c),2,:)),'r')
end

figure
useInh= find(dataInh & layerAll==5)
for c=1:length(useInh)
    set (gcf,'Name','Lyr 5 inh psth/unit')
    subplot(3,5,c)
    plot(squeeze(mnPsth(useInh(c),1,:)));axis square;hold on;
    plot(squeeze(mnPsth(useInh(c),2,:)),'r')
end

figure
useInh= find(dataInh & layerAll==5)
for c=1:length(useInh)
    set(gcf,'Name','hist of change prepost lyr5 inh psth')
    MImnPsth= (mnPsth(useInh,2,:)-mnPsth(useInh,1,:))./(mnPsth(useInh,2,:)+mnPsth(useInh,1,:));
    hist(MImnPsth)
    %subplot(4,2,t)
    h= hist(MImnPsth(useInh),-1:.1:1);
    Mbins=-1:.1:1
    bar(Mbins,h/sum(useInh))
end

titles = {'Saline','DOI','5HT','Ketanserin', 'Ketanserin + DOI', 'MGluR2','MGlur2 + DOI', 'Lisuride'};
figure
for t=1:8
        subplot(2,4,t)
    plot(drift_osi(data & treatment==t,1,1),drift_osi(data & treatment==t,1,2),'.')
    hold on; axis square
    plot([-1.5 1.5], [-1.5 1.5]);
    title(titles{t})
end



% meanOSI = [mean(drift_osi(data,1,1))  mean(drift_osi(data,1,2))];
% figure
% hb = bar(meanOSI)
% plot(hb)

amp = squeeze(max(drift_orient,[],2));
useResp = amp(:,1,1)>2 %non-selective cells

figure
for l=1:6
osi(2,:) = squeeze(nanmean(drift_osi(data & layerAll==l,1,2),1))
osierr(2,:) = squeeze(nanstd(drift_osi(data & layerAll==l,1,2),1))/sqrt(sum(drift_osi(useN & layerAll==l,1,2)));
osi(1,:) = squeeze(nanmean(drift_osi(data & layerAll==l,1,1),1))
osierr(1,:) = squeeze(nanstd(drift_osi(data & layerAll==l,1,1),1))/sqrt(sum(drift_osi(useN & layerAll==l,1,1)));
subplot(2,3,l)
% barweb(osi,osierr)
bar(osi);hold on;
end

figure
subplot(2,2,1); imagesc(squeeze(wn_frameR(find(data),:,1,1)),[ 0 50]); title('pre stop')
subplot(2,2,2); imagesc(squeeze(wn_frameR(find(data),:,1,2)),[ 0 50]); title('post stop')
subplot(2,2,3); imagesc(squeeze(wn_frameR(find(data),:,2,1)),[ 0 50]); title('pre move');
subplot(2,2,4); imagesc(squeeze(wn_frameR(find(data),:,2,2)),[ 0 50]); title('post move')
set(gcf,'Name',sprintf('%s %s',files(use(i)).expt, files(use(i)).treatment));
drawnow

      
figure
subplot(2,2,1); plot(squeeze(nanmean(wn_frameR(find(data),:,1,1),1))); title('pre stop'); ylim([0 5])
subplot(2,2,2); plot(squeeze(nanmean(wn_frameR(find(data),:,1,2),1))); title('post stop');ylim([0 5])
subplot(2,2,3); plot(squeeze(nanmean(wn_frameR(find(data),:,2,1),1))); title('pre move'); ylim([0 5])
subplot(2,2,4); plot(squeeze(nanmean(wn_frameR(find(data),:,2,2),1))); title('post move'); ylim([0 5])
set(gcf,'Name',sprintf('%s %s',files(use(i)).expt, files(use(i)).treatment));
drawnow

  clear cycR
for f = 1:20;
    cycR(:,f,:,:) = nanmean(wn_frameR(:,f:20:end,:,:),2);
end

figure
subplot(2,2,1); imagesc(squeeze(cycR(find(data),:,1,1)),[ 0 50]); title('pre stop')
subplot(2,2,2); imagesc(squeeze(cycR(find(data),:,1,2)),[ 0 50]); title('post stop')
subplot(2,2,3); imagesc(squeeze(cycR(find(data),:,2,1)),[ 0 50]); title('pre move');
subplot(2,2,4); imagesc(squeeze(cycR(find(data),:,2,2)),[ 0 50]); title('post move')

figure
subplot(2,2,1); plot(squeeze(mean(cycR(find(data),:,1,1),1))); title('pre stop'); ylim([1.5 4])
subplot(2,2,2); plot(squeeze(mean(cycR(find(data),:,1,2),1))); title('post stop'); ylim([1.5 4])
subplot(2,2,3); plot(squeeze(mean(cycR(find(data),:,2,1),1))); title('pre move'); ylim([1.5 4])
subplot(2,2,4); plot(squeeze(mean(cycR(find(data),:,2,2),1))); title('post move'); ylim([1.5 4])

spont = squeeze(mean(cycR(:,[1 2 19 20],:,:),2));
evoked = squeeze(mean(cycR(:,9:11,:,:),2)) - spont;

titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
for t=1:8
   figure
   set(gcf,'Name',titles{t});
for lyr = 2:6
    for i=1:2
        if i==1
            use = find(goodAll & layerAll==lyr & ~inhAll & treatment==t); symb = 'bo';
        else
            use = find(goodAll & layerAll==lyr  & inhAll & treatment==t); symb = 'ro';
        end
        subplot(5,4,1 + 4*(lyr-2));
        plot(spont(use,1,1),spont(use,1,2),symb); hold on; plot([0 50],[0 50]);  axis square; axis([0 30 0 30])
        if lyr==2,  title({'spont stop','layer 2'}), end; if lyr>2, title(sprintf('layer %d',lyr)); end
        subplot(5,4,2+ 4*(lyr-2));
        plot(spont(use,2,1),spont(use,2,2),symb);  hold on; plot([0 50],[0 50]); axis square; axis([0 30 0 30])
        if lyr==2,  title({'spont move',''}), end
        subplot(5,4,3+ 4*(lyr-2));
        plot(evoked(use,1,1),evoked(use,1,2),symb); hold on; plot([-30 30],[-30 30]);  axis square; axis([-10 20 -10 20])
        if lyr==2,  title({'evoked stop',''}), end
        subplot(5,4,4+ 4*(lyr-2))
        plot(evoked(use,2,1),evoked(use,2,2),symb); hold on; plot([-30 30],[-30 30]);  axis square; axis([-10 20 -10 20])
        if lyr==2,  title({'evoked move',''}), end
    end
end
end


for lyr = 2:6
    figure
    if lyr<6
        use = find(goodAll & layerAll==lyr & ~inhAll); 
        set(gcf,'Name',sprintf('layer %d',lyr));symb = 'bo';
    elseif lyr==6
          use = find(goodAll & layerAll==lyr & ~inhAll);
          set(gcf,'Name',sprintf('layer %d',lyr));symb = 'bo'; 
    else
        use = find(goodAll & inhAll);
        set(gcf,'Name',sprintf('inh',lyr));symb = 'ro';
    end
    subplot(2,2,1);
    plot(spont(use,1,1),spont(use,2,1),symb); title('spont pre stop vs move'); hold on; plot([0 50],[0 50]);  axis square; axis([0 30 0 30])
    subplot(2,2,3);
    plot(spont(use,1,2),spont(use,2,2),symb); title('spont post stop vs move'); hold on; plot([0 50],[0 50]); axis square; axis([0 30 0 30])
    subplot(2,2,2);
    plot(evoked(use,1,1),evoked(use,2,1),symb); title('evoked pre stop vs move'); hold on; plot([-30 30],[-30 30]);  axis square; axis([-10 20 -10 20])
    subplot(2,2,4)
    plot(evoked(use,1,2),evoked(use,2,2),symb); title('evoked post stop vs move'); hold on; plot([-30 30],[-30 30]);  axis square; axis([-10 20 -10 20])
end

titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
figure
for t=1:8
set(gcf,'Name', 'spont pre vs post stop')
subplot(2,4,t)
plot(spont(find(goodAll & treatment==t),1,1),spont(find(goodAll& treatment==t),1,2),'o'); title(titles{t}); hold on; plot([0 50],[0 50]);
plot(spont(find(goodAll&inhAll& treatment==t),1,1),spont(find(goodAll & inhAll& treatment==t),1,2),'ro'); hold on; plot([0 50],[0 50]); axis([0 15 0 15])
end

figure
for t=1:8
set(gcf,'Name', 'spont pre vs post stop')
subplot(2,4,t)
plot(spont(find(goodAll& treatment==t),2,1),spont(find(goodAll& treatment==t),2,2),'o'); title(titles{t}); hold on; plot([0 50],[0 50]);
plot(spont(find(goodAll&inhAll& treatment==t),2,1),spont(find(goodAll & inhAll& treatment==t),2,2),'ro'); hold on; plot([0 50],[0 50]); axis([0 15 0 15])
end

figure
for t=1:8
set(gcf,'Name', 'evoked pre vs post stop')
subplot(2,4,t)
plot(evoked(find(goodAll& treatment==t),1,1),evoked(find(goodAll& treatment==t),1,2),'o'); title(titles{t}); hold on; plot([0 50],[0 50]);
plot(evoked(find(goodAll&inhAll& treatment==t),1,1),evoked(find(goodAll & inhAll& treatment==t),1,2),'ro'); hold on; plot([0 50],[0 50]); axis([-10 10 -10 10])
end

figure
for t=1:8
set(gcf,'Name', 'evoked pre vs post move')
subplot(2,4,t)
plot(evoked(find(goodAll& treatment==t),2,1),evoked(find(goodAll& treatment==t),2,2),'o'); title(titles{t}); hold on; plot([0 50],[0 50]);
plot(evoked(find(goodAll&inhAll& treatment==t),2,1),evoked(find(goodAll & inhAll& treatment==t),2,2),'ro'); hold on; plot([0 50],[0 50]); axis([-10 10 -10 10])
end


figure
for t=1:8
set(gcf,'Name', 'spont stop vs move pre')
subplot(2,4,t)
plot(spont(find(goodAll & treatment==t),1,1),spont(find(goodAll &treatment==t),2,1),'o'); title(titles{t}); hold on; plot([0 50],[0 50]);
plot(spont(find(goodAll&inhAll &treatment==t),1,1),spont(find(goodAll & inhAll &treatment==t),2,1),'ro'); hold on; plot([0 50],[0 50]); axis([0 50 0 50])
end

figure
for t=1:8
set(gcf,'Name', 'evoked stop vs move pre')
subplot(2,4,t)
plot(evoked(find(goodAll & treatment==t),1,1),evoked(find(goodAll& treatment==t),2,1),'o'); title(titles{t}); hold on; plot([0 50],[0 50]);
plot(evoked(find(goodAll&inhAll& treatment==t),1,1),evoked(find(goodAll & inhAll& treatment==t),2,1),'ro'); hold on; plot([0 50],[0 50]); axis([0 50 0 50 ])
end

figure
for t=1:8
set(gcf,'Name', 'spont stop vs move post')    
subplot(2,4,t)
plot(spont(find(goodAll& treatment==t),1,2),spont(find(goodAll& treatment==t),2,2),'o'); title(titles{t}); hold on; plot([0 50],[0 50]);
plot(spont(find(goodAll&inhAll& treatment==t),1,2),spont(find(goodAll & inhAll& treatment==t),2,2),'ro'); hold on; plot([0 50],[0 50]); axis([0 50 0 50])
end

figure
for t=1:8
set(gcf,'Name', 'evoked stop vs move post')
subplot(2,4,t)
plot(evoked(find(goodAll& treatment==t),1,2),evoked(find(goodAll& treatment==t),2,2),'o'); title(titles{t}); hold on; plot([0 50],[0 50]);
plot(evoked(find(goodAll&inhAll& treatment==t),1,2),evoked(find(goodAll & inhAll& treatment==t),2,2),'ro');hold on; plot([0 50],[0 50]); axis([0 50 0 50])
end

figure
for t=1:8
subplot(2,4,t)
set(gcf,'Name', 'latent')
pre = [cycR(find(data & treatment==t),:,1,1) cycR(find(data & treatment==t),:,2,1)];
post = [cycR(find(data & treatment==t),:,1,2) cycR(find(data & treatment==t),:,2,2)];
[coeff score latent] = pca(pre','algorithm','als','variableweight','variance');
plot(latent(1:10)/sum(latent));title(titles{t});

end

% how much does each cell contribute to component
figure 
for t=1:8
subplot(2,4,t)
set(gcf,'Name', 'coeff')
pre = [cycR(find(goodAll & ~inhAll & treatment==t),:,1,1) cycR(find(goodAll & ~inhAll& treatment==t),:,2,1)];
post = [cycR(find(goodAll & ~inhAll& treatment==t),:,1,2) cycR(find(goodAll & ~inhAll& treatment==t),:,2,2)];
[coeff score latent] = pca(pre','algorithm','als','variableweight','variance');
imagesc(coeff(:,1:10));title(titles{t});
end

for i = 1:10; %pre stationary then moving
    subplot(10,1,i);
    plot(score(:,i))
end

%============
%first two components plotted...stationary then moving
figure
plot(score(:,1)); hold on; plot(score(:,2),'r')
figure
plot(score(:,1)); hold on; plot(score(:,3),'r')

preS = pre'*coeff;
postS = post'*coeff;

%using 8 components...can change to however many you want to look at
figure
for i = 1:8
    subplot(8,1,i)
    plot(preS(:,i)); hold on; plot(postS(:,i),'r');
end


figure
plot(preS(1:20,1),preS(1:20,2),'b'); hold on;plot(preS(21:40,1),preS(21:40,2),'b');
plot(postS(1:20,1),postS(1:20,2),'r'); plot(postS(21:40,1),postS(21:40,2),'r');


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
titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI','MGluR2', 'MGluR2 + DOI','Lisuride'};
figure
for i = 1:8
    subplot(2,4,i);
    plot(wnCorr(corrTreatment==i,1),wnCorr(corrTreatment==i,2),'.'); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
    xlabel('pre wn corr'); ylabel('post')
    set(gcf,'Name','Wn Corr')
end

% titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI', 'MGlur2', 'MGlur2 + DOI', 'Lisuride'};
% figure
% for i = 1:7
%     subplot(2,4,i);
%     wnCorrHist= myHist2(wnCorr(corrTreatment==i,1),wnCorr(corrTreatment==i,2),-.5:.1:1.5,-.5:.1:1.5);
%     %wnCorrHist_pre= myHist2(wnCorr(corrTreatment==i,1),-.5:.1:1.5,-.5:.1:1.5)
%     plot(wnCorrHist);hold on; axis square;title(titles{i})
%     set(gcf,'Name','Wn Corr'); ylim([0 6000])
% end


%%% plot correlation for darkness
titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI','MGlur2','MGlur2 + DOI', 'Lisuride'};
figure
for i = 1:8
    subplot(2,4,i);
    plot(darkCorr(corrTreatment==i,1),darkCorr(corrTreatment==i,2),'.'); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
    xlabel('pre dark corr'); ylabel('post')
    set(gcf,'Name','Dark Corr')
end

% titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI','MGluR2','MGlur2 + DOI','Lisuride'};
% figure
% for i = 1:8
%     subplot(2,4,i);
%     darkCorrHist= myHist2(darkCorr(corrTreatment==i,1),darkCorr(corrTreatment==i,2),-.5:.1:1.5,-.5:.1:1.5);
%     plot(darkCorrHist);hold on; axis square;title(titles{i});
%     set(gcf,'Name','Dark Corr'); ylim ([0 6000])
% 
% end

%%% compare spontaneous rates measured with gratings and wn
figure
plot(drift_spont(:),wn_spont(:),'.'); hold on; plot([0 10], [0 10]); axis equal
xlabel('drift spont'),ylabel('wn spont');


%%% scatter plot of drift spont
for mv = 1:2
    figure 
    for i = 1:6
        subplot(2,3,i)
        plot(drift_spont(goodAll & treatment==saline & layerAll ==i,mv,1),drift_spont(goodAll &treatment==saline& layerAll ==i,mv,2),'k.');
        hold on
        plot(drift_spont(goodAll & treatment==doi& layerAll ==i,mv,1),drift_spont(goodAll &treatment==doi& layerAll ==i,mv,2),'r.');
        plot(drift_spont(goodAll & treatment==ht& layerAll ==i,mv,1),drift_spont(goodAll &treatment==ht& layerAll ==i,mv,2),'w.');
        plot(drift_spont(goodAll &treatment==ketanserin& layerAll ==i,mv,1),drift_spont(goodAll &treatment==ketanserin& layerAll ==i,mv,2),'m.');
        plot(drift_spont(goodAll &treatment==ketandoi& layerAll ==i,mv,1),drift_spont(goodAll &treatment==ketandoi& layerAll ==i,mv,2),'c.');
        plot(drift_spont(goodAll &treatment==mglur2& layerAll ==i,mv,1),drift_spont(goodAll &treatment==mglur2& layerAll ==i,mv,2),'g.');
        plot(drift_spont(goodAll &treatment==mglur2doi& layerAll ==i,mv,1),drift_spont(goodAll &treatment==mglur2doi& layerAll ==i,mv,2),'y.');
         plot(drift_spont(goodAll &treatment==lisuride& layerAll ==i,mv,1),drift_spont(goodAll &treatment==lisuride& layerAll ==i,mv,2),'b.');
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
        plot(wn_spont(goodAll &treatment==saline & layerAll ==i,mv,1),wn_spont(goodAll &treatment==saline& layerAll ==i,mv,2),'k.');
        % ylim ([-2 50]); xlim([-2 50]);
        hold on
        plot(wn_spont(goodAll &treatment==doi& layerAll ==i,mv,1),wn_spont(goodAll &treatment==doi& layerAll ==i,mv,2),'r.');
        plot(wn_spont(goodAll &treatment==ht& layerAll ==i,mv,1),wn_spont(goodAll &treatment==ht& layerAll ==i,mv,2),'w.');
        plot(wn_spont(goodAll &treatment==ketanserin& layerAll ==i,mv,1),wn_spont(goodAll &treatment==ketanserin& layerAll ==i,mv,2),'m.');
        plot(wn_spont(goodAll &treatment==ketandoi& layerAll ==i,mv,1),wn_spont(goodAll &treatment==ketandoi& layerAll ==i,mv,2),'c.');
        plot(wn_spont(goodAll &treatment==mglur2& layerAll ==i,mv,1),wn_spont(goodAll &treatment==mglur2& layerAll ==i,mv,2),'g.');
        plot(wn_spont(goodAll &treatment==lisuride& layerAll ==i,mv,1),wn_spont(goodAll &treatment==lisuride& layerAll ==i,mv,2),'b.');
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
        plot(wn_evoked(goodAll &treatment==saline & layerAll ==i,mv,1),wn_evoked(goodAll &treatment==saline& layerAll ==i,mv,2),'k.');
        hold on
        plot(wn_evoked(goodAll &treatment==doi& layerAll ==i,mv,1),wn_evoked(goodAll &treatment==doi& layerAll ==i,mv,2),'r.');
        plot(wn_evoked(goodAll &treatment==ht& layerAll ==i,mv,1),wn_evoked(goodAll &treatment==ht& layerAll ==i,mv,2),'w.');
        plot(wn_evoked(goodAll &treatment==ketanserin& layerAll ==i,mv,1),wn_evoked(goodAll &treatment==ketanserin& layerAll ==i,mv,2),'m.');
        plot(wn_evoked(goodAll &treatment==ketandoi& layerAll ==i,mv,1),wn_evoked(goodAll &treatment==ketandoi& layerAll ==i,mv,2),'c.');
        plot(wn_evoked(goodAll &treatment==mglur2& layerAll ==i,mv,1),wn_evoked(goodAll &treatment==mglur2& layerAll ==i,mv,2),'g.');
        plot(wn_evoked(goodAll &treatment==lisuride& layerAll ==i,mv,1),wn_evoked(goodAll &treatment==lisuride& layerAll ==i,mv,2),'b.');
        plot([0 10],[0 10]); xl = get(gca,'Xlim'); yl = get(gca,'Ylim'); axis square; %axis([min(xl(1),yl(1)) max(xl(2),yl(2)) min(xl(1),yl(1)) max(xl(2),yl(2)) ])
        title(sprintf('layer %d',i));  ylabel('post'); ylim([-5 50]);xlim([-5 50])
        if mv ==1 , xlabel('stop wn evoked'); else  xlabel('move wn evoked'); end
    end
end

%rows = layers col = treatment...wn evoked FR mv and stationary
titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
for mv = 1:2
    figure
    if mv ==1 , set(gcf,'Name','stop wn evoked'); else  set(gcf,'Name','move wn evoked'); end
    for c=0:1
    for t=1:8
        subplot(7,7,t)
        plot(wn_evoked(goodAll &treatment==t & layerAll ==1 & inhAll==c,mv,1),wn_evoked(goodAll &treatment==t& layerAll ==1& inhAll==c,mv,2),'.');axis equal;ylim([-5 20]);xlim([-5 20])
        hold on; plot([0 45],[0 45])
        subplot(7,7,t+7)
        plot(wn_evoked(goodAll &treatment==t & layerAll ==2& inhAll==c,mv,1),wn_evoked(goodAll &treatment==t& layerAll ==2& inhAll==c,mv,2),'.');axis equal;ylim([-5 30]);xlim([-5 30])
        hold on; plot([0 45],[0 45])
        subplot(7,7,t+14)
        plot(wn_evoked(goodAll &treatment==t & layerAll ==3& inhAll==c,mv,1),wn_evoked(goodAll &treatment==t& layerAll ==3& inhAll==c,mv,2),'.');axis equal;ylim([-5 40]);xlim([-5 40])
        hold on; plot([0 45],[0 45])
        subplot(7,7,t+21)
        plot(wn_evoked(goodAll &treatment==t & layerAll ==4& inhAll==c,mv,1),wn_evoked(goodAll &treatment==t& layerAll ==4& inhAll==c,mv,2),'.');axis equal;ylim([-5 40]);xlim([-5 40])
        hold on; plot([0 45],[0 45])
        subplot(7,7,t+28)
        plot(wn_evoked(goodAll &treatment==t & layerAll ==5& inhAll==c,mv,1),wn_evoked(goodAll &treatment==t& layerAll ==5& inhAll==c,mv,2),'.');axis equal;ylim([-5 30]);xlim([-5 30])
        hold on; plot([0 45],[0 45])
        subplot(7,7,t+32)
        plot(wn_evoked(goodAll &treatment==t & layerAll ==6& inhAll==c,mv,1),wn_evoked(goodAll &treatment==t& layerAll ==6& inhAll==c,mv,2),'.');axis equal; ylim([-5 20]);xlim([-5 20]);hold on
        plot([0 45],[0 45]); %axis([min(xl(1),yl(1)) max(xl(2),yl(2)) min(xl(1),yl(1)) max(xl(2),yl(2)) ])
        ylabel('post'); xlabel('pre');
    end 
end
end

titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
for mv = 1:2
    figure
    if mv ==1 , set(gcf,'Name','stop wn evoked'); else  set(gcf,'Name','move wn evoked'); end
    for c=0:1
    for t=1:8
        subplot(6,7,t)
        plot(wn_evoked(goodAll &treatment==t & layerAll ==1 & inhAll==c,mv,1),wn_evoked(goodAll &treatment==t& layerAll ==1& inhAll==c,mv,2),'.');axis equal;ylim([-5 45]);xlim([-5 45])
        hold on; plot([0 45],[0 45])
        subplot(6,7,t+6)
        plot(wn_evoked(goodAll &treatment==t & layerAll ==2& inhAll==c,mv,1),wn_evoked(goodAll &treatment==t& layerAll ==2& inhAll==c,mv,2),'.');axis equal;ylim([-5 45]);xlim([-5 45])
        hold on; plot([0 45],[0 45])
        subplot(6,7,t+12)
        plot(wn_evoked(goodAll &treatment==t & layerAll ==3& inhAll==c,mv,1),wn_evoked(goodAll &treatment==t& layerAll ==3& inhAll==c,mv,2),'.');axis equal;ylim([-5 45]);xlim([-5 45])
        hold on; plot([0 45],[0 45])
        subplot(6,7,t+18)
        plot(wn_evoked(goodAll &treatment==t & layerAll ==4& inhAll==c,mv,1),wn_evoked(goodAll &treatment==t& layerAll ==4& inhAll==c,mv,2),'.');axis equal;ylim([-5 45]);xlim([-5 45])
        hold on; plot([0 45],[0 45])
        subplot(6,7,t+24)
        plot(wn_evoked(goodAll &treatment==t & layerAll ==5& inhAll==c,mv,1),wn_evoked(goodAll &treatment==t& layerAll ==5& inhAll==c,mv,2),'.');axis equal;ylim([-5 45]);xlim([-5 45])
        hold on; plot([0 45],[0 45])
        subplot(6,7,t+32)
        plot(wn_evoked(goodAll &treatment==t & layerAll ==6& inhAll==c,mv,1),wn_evoked(goodAll &treatment==t& layerAll ==6& inhAll==c,mv,2),'.');axis equal; ylim([-5 45]);xlim([-5 45]);hold on
        plot([0 45],[0 45]); %axis([min(xl(1),yl(1)) max(xl(2),yl(2)) min(xl(1),yl(1)) max(xl(2),yl(2)) ])
        ylabel('post'); xlabel('pre');
    end 
end
end

titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
%dark layers and treatments
    figure
    for c=0:1
    for t=1:8
        subplot(6,7,t)
        plot(meanRdark(goodAll &treatment==t & layerAll ==1 & inhAll==c,1),meanRdark(goodAll &treatment==t& layerAll ==1& inhAll==c,2),'.');axis equal;ylim([0 20]);xlim([0 20])
        hold on; plot([0 45],[0 45]);
        subplot(6,7,t+7)
        plot(meanRdark(goodAll &treatment==t & layerAll ==2& inhAll==c,1),meanRdark(goodAll &treatment==t& layerAll ==2& inhAll==c,2),'.');axis equal;ylim([0 20]);xlim([0 20])
        hold on; plot([0 45],[0 45])
        subplot(6,7,t+14)
        plot(meanRdark(goodAll &treatment==t & layerAll ==3& inhAll==c,1),meanRdark(goodAll &treatment==t& layerAll ==3& inhAll==c,2),'.');axis equal;ylim([0 40]);xlim([0 40])
        hold on; plot([0 45],[0 45])
        subplot(6,7,t+21)
        plot(meanRdark(goodAll &treatment==t & layerAll ==4& inhAll==c,1),meanRdark(goodAll &treatment==t& layerAll ==4& inhAll==c,2),'.');axis equal;ylim([0 45]);xlim([0 45])
        hold on; plot([0 45],[0 45])
        subplot(6,7,t+28)
        plot(meanRdark(goodAll &treatment==t & layerAll ==5& inhAll==c,1),meanRdark(goodAll &treatment==t& layerAll ==5& inhAll==c,2),'.');axis equal;ylim([0 45]);xlim([0 45])
        hold on; plot([0 45],[0 45])
        subplot(6,7,t+32)
        plot(meanRdark(goodAll &treatment==t & layerAll ==6& inhAll==c,1),meanRdark(goodAll &treatment==t& layerAll ==6& inhAll==c,2),'.');axis equal; ylim([0 20]);xlim([0 20]);hold on
        plot([0 45],[0 45]); %axis([min(xl(1),yl(1)) max(xl(2),yl(2)) min(xl(1),yl(1)) max(xl(2),yl(2)) ])
        ylabel('post'); xlabel('pre');
    end 
    end
    

titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};%scatter all units/treatment prepost
for mv=1:2
    figure
    if mv==1 set(gcf,'Name','stationary'), else set(gcf,'Name','moving');end
   for c=0:1
    for t=1:8
        subplot(4,2,t)
        plot(wn_evoked(goodAll &treatment==t & inhAll ==c,mv,1),wn_evoked(goodAll &treatment==t & inhAll==c ,mv,2),'.');hold on; plot([0 45],[0 45]);axis equal; ylim([-5 40]);xlim([-5 40]);
        % plot(wn_evoked(treatment==t & inhAll(useN(i)) ==1,mv,1),wn_evoked(treatment==t &inhAll(useN(i)) ==1 ,mv,2),'r.');hold on; plot([0 10],[0 10]);axis equal; ylim([-5 40])
        ylabel('post'); xlabel('pre')
        hold on ;title(titles{t});
    end
   end
end
   
%mean prepost darkness fr
figure
for t = 1:8
    for c=0:1
subplot(4,2,t)
plot(meanRdark(goodAll &treatment==t & inhAll==c,1),meanRdark(goodAll &treatment==t & inhAll==c,2),'.');hold on; plot([0 45],[0 45]);axis equal;ylim([0 45]); xlim([0 45]);
title(titles{t});ylabel('post'); xlabel('pre')
end
end


%%% plot white noise response functions for all units
for t = 1:8
    figure
    if t==1, set(gcf,'Name','saline wn CRF'),
    elseif t==2, set(gcf,'Name','doi wn CRF'),
    elseif t==3, set(gcf,'Name','ht wn CRF'),
    elseif t==4, set(gcf,'Name','ketanserin wn CRF')
    elseif t==5, set(gcf,'Name','ketanserin + doi wn CRF')
    elseif t==6, set(gcf,'Name','MGluR2 wn CRF')
    elseif t==7, set(gcf, 'Name','MGluR2+ DOI wn CRF')
    else set(gcf,'Name','Lisuride'),end
    
    useN = find(goodAll &treatment==t)
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
for t = 1:8
    figure
    if t==1, set(gcf,'Name','saline OT'),
    elseif t==2, set(gcf,'Name','doi OT'),
    elseif t==3, set(gcf,'Name','ht OT'),
    elseif t==4, set(gcf,'Name','ketanserin OT')
    elseif t==5, set(gcf,'Name','ketanserin + doi OT')
    elseif t==6, set(gcf,'Name','MGlur2 OT')
    elseif t==7, set(gcf,'Name','MGlur2+DOI OT')
    else set(gcf,'Name','Lisuride'),end
    useN = find(goodAll & treatment==t)
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
for t = 1:8
    figure
    if t==1, set(gcf,'Name','saline SF'),
    elseif t==2, set(gcf,'Name','doi SF'),
    elseif t==3, set(gcf,'Name','ht SF'),
    elseif t==4, set(gcf,'Name','ketanserin SF')
    elseif t==5, set(gcf,'Name','ketanserin + doi SF')
    elseif t==6, set(gcf,'Name','mglur2 SF')
    elseif t==7, set(gcf,'Name','mglur2 + DOI SF')
    else set(gcf,'Name','lisuride SF'),end
    useN = find(goodAll & treatment==t)
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

% amp = squeeze(max(drift_orient,[],2));
% useOS = amp(:,1,1)>2 %non-selective cells
% titles = {'Saline','DOI','5HT','Ketanserin', 'Ketanserin + DOI', 'MGluR2','MGlur2 + DOI', 'Lisuride'};
% figure
% for t=1:8
%     subplot(2,4,t)
% plot(drift_osi(goodAll==1 & useOS & treatment==t,:,1),drift_osi(goodAll==1 & useOS & treatment==t,:,2),'.')
% hold on;
% plot([-1.5 1.5], [-1.5 1.5]);
% title(titles{t})
% %h=hist(drift_osi(goodAll==1,:,:),0:.01:1)
% % figure
% % plot(h)
% end


figure
for t=1:8
    set(gcf,'Name','drift osi dist pre/post')
    subplot(2,4,t)
    h= hist(drift_osi(goodAll==1 & treatment==t,:,1),-1:.1:1);
    h2= hist(drift_osi(goodAll==1 & treatment==t,:,2),-1:.1:1);
    Mbins=-1:.1:1
    bar(Mbins,h/sum(goodAll&treatment==t));hold on
    bar(Mbins,h2/sum(goodAll & treatment==t))
end

figure
set(gcf,'Name','drift osi dist pre/post')
for t=1:8
    subplot(2,4,t)
    h= [ mean(abs(drift_osi(goodAll==1 & treatment==t & layerAll==5,1,1))) mean(abs(drift_osi(goodAll==1 & treatment==t & layerAll==5,1,2)))];
    bar(h)
    set(gca,'XTickLabel',{'pre' 'post'})
    title(titles{t})
    %     Mbins=-1:.1:1
    %     bar(Mbins,h/sum(goodAll&treatment==t))
    %     bar(Mbins,h2/sum(goodAll & treatment==t))
end


% 
% for t = 1:4
%     figure
%     figure
%     if t==1, set(gcf,'Name','saline ev LFP'),
%     elseif t==2, set(gcf,'Name','doi ev LFP'),
%     elseif t==3, set(gcf,'Name','ketanserin ev LFP')
%     else set(gcf,'Name','ketanserin + doi ev LFP'),end
%     useN = find(treatment==t)
%     for i = 1:length(useN)
%         np = ceil(sqrt(length(useN)));
%         subplot(np,np,i);
%         hold on
%         plot(LFPall(useN(i),:,1,1),'Color',[0.5 0 0]);  plot(LFPall(:,:,2,1),'Color',[0 0.5 0]); %pre, mv=2
%         plot(LFPall(i,:,1,2),'Color',[1 0 0]);plot(LFPall(i,:,2,2),'Color',[0 1 0]); %post
%        %         yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)]); xlim([0.5 8.5])
% 
%     end
% end


%%% plot speed histogram
figure
hold on
subplot(4,2,1)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==saline,1:25,:),1))); title('saline'); xlabel('speed')
subplot(4,2,2)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==doi,1:25,:),1))); title('doi'); xlabel('speed')
subplot(4,2,3)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==ht,1:25,:),1))); title('ht'); xlabel('speed')
subplot(4,2,4)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==ketanserin,1:25,:),1))); title('ketanserin'); xlabel('speed')
subplot(4,2,5)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==ketandoi,1:25,:),1))); title('ketanserin + DOI'); xlabel('speed')
subplot(4,2,6)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==mglur2,1:25,:),1))); title('mglur2'); xlabel('speed')
subplot(4,2,7)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==mglur2doi,1:25,:),1))); title('mglur2 + DOI'); xlabel('speed')
subplot(4,2,8)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==lisuride,1:25,:),1))); title('Lisuride'); xlabel('speed')
legend('pre','post')


% % right now there's an error where units between treatments seem to be
% % mixed...running all saline alone gives a different MI distribution than
% % when all treatments are run together
titles = {'Saline','DOI','5HT','Ketanserin', 'Ketanserin + DOI', 'MGluR2','MGlur2 + DOI', 'Lisuride'};
figure
for t = 1:8
    useDark = meanRdark(:,1)>.3 | meanRdark(:,2)>.3;
    useN = find(goodAll & treatment==t)
    MIdark= (meanRdark(:,2)-meanRdark(:,1))./(meanRdark(:,2)+meanRdark(:,1));
    
    subplot(5,2,t+2)
    h= hist(MIdark(useDark(treatment==t)),-1:.1:1);
    Mbins=-1:.1:1
    bar(Mbins,h/sum(useDark(treatment==t)))
    xlim([-1.5 1.5]); ylim([0 .25]);axis xy
    xlabel('MI'); ylabel('fraction of cells');title(titles{t});
    set(gcf,'Name','MI Dark')
    %hold on
    subplot(5,2,2)
    meanMIdark(t) = nanmean(MIdark(treatment==t))
    %err(t) = nanstd(MIdark(treatment==t))/sqrt(sum(MIdark (treatment==t)));
    bar(meanMIdark)
    %barweb(meanMIdark(1:1:t)',err(1:1:t)');
    ylim([-0.2 0.2]);
end
% 
% %% mod for corr in darkness %%

titles = {'Saline','DOI','5HT','Ketanserin', 'Ketanserin + DOI','MGluR2', 'MGluR2 + DOI', 'Lisuride'};
figure
for t = 1:8
    useN = find(goodAll & treatment==t)
    MIdarkCorr= (darkCorr(:,2)-darkCorr(:,1))./(darkCorr(:,2)+darkCorr(:,1));
    subplot(4,2,t)
    h= hist(MIdarkCorr(useN),-1:.1:1);
    Mbins=-1:.1:1
    bar(Mbins,h/sum(useN))
    xlim([-1.5 1.5]); ylim([0 .25]);axis xy
    xlabel('Dark pairwise MI'); ylabel('fraction of cells');title(titles{t});
    set(gcf,'Name','Corr MI Dark')
end

% % mod for WN corr
% titles = {'Saline','DOI','5HT','Ketanserin', 'Ketanserin + DOI','MGluR2', 'MGluR2 + DOI'};
% figure
% for t = 1:8
%     useN = find(treatment==t)
%     MIwnCorr= (wnCorr(:,2)-wnCorr(:,1))./(wnCorr(:,2)+wnCorr(:,1));
%     subplot(5,2,t)
%     h= hist(MIwnCorr(useN),-1:.1:1);
%     Mbins=-1:.1:1
%     bar(Mbins,h/sum(useN))
%     xlim([-1.5 1.5]); ylim([0 .35]);axis xy
%     xlabel('wn pairwise MI'); ylabel('fraction of cells');title(titles{t});
%     set(gcf,'Name','Corr MI Wn')
% end
% 
% % MI MOVING! WN evoked 
% figure
% for t = 1:8
%     thresh = 1;
%     useEv =wn_evoked(:,2,1)>0 & wn_evoked(:,2,2)>0 & (wn_evoked(:,2,1)>thresh | wn_evoked(:,2,2)>thresh) & treatment'==t;
%     useN = find(treatment==t)
%     MI_mv_wn = (wn_evoked(:,2,2)-wn_evoked(:,2,1))./(wn_evoked(:,2,2)+wn_evoked(:,2,1));
%     subplot(3,3,t+3)
%     hEv= hist(MI_mv_wn(useEv(treatment==t)),-1:.1:1);
%     Mbins=-1:.1:1
%     bar(Mbins,hEv/sum(useEv(treatment==t)))
%     xlabel('MI'); ylabel('fraction of cells'); xlim([-1.5 1.5]);ylim([0 .3]);
%     set(gcf,'Name','MI move wn')
%     title(titles{t});
%     subplot(3,2,2)
%     meanMI_mv_wn(t) = nanmean(MI_mv_wn(treatment==t))
%     bar(meanMI_mv_wn);ylim([-1 1]);
% end
% 
% %MI wn EVOKED stationary:
% figure
% for t = 1:8
%     thresh = 1;
%     useEv =wn_evoked(:,1,1)>0 & wn_evoked(:,1,2)>0 & (wn_evoked(:,1,1)>thresh | wn_evoked(:,1,2)>thresh) & treatment'==t;
%     useN = find(treatment==t)
%     MI_stat_wn = (wn_evoked(:,1,2)-wn_evoked(:,1,1))./(wn_evoked(:,1,2)+wn_evoked(:,1,1));
%     subplot(3,3,t+3)
%     hEv_stat= hist(MI_stat_wn(useEv(treatment==t)),-1:.1:1);
%     Mbins=-1:.1:1
%     bar(Mbins,hEv_stat/sum(useEv(treatment==t)))
%     xlabel('MI'); ylabel('fraction of cells'); xlim([-1.5 1.5]);ylim([0 .3]);
%     set(gcf,'Name','MI stationary wn')
%     title(titles{t});
%     subplot(3,3,2)
%     meanMI_stat_wn(t) = nanmean(MI_stat_wn(treatment==t))
%     bar(meanMI_stat_wn);ylim([-1 1]);
% end
% 
% % MI wn SPONTANEOUS stationary
% figure
% for t = 1:8
%     thresh = 1;
%     useSpont = wn_spont(:,1,1)>0 & wn_spont(:,1,2)>0 & (wn_spont(:,1,1)>thresh | wn_spont(:,1,2)>thresh) & treatment'==t;
%     MI_stat_spont = (wn_spont(:,1,2)-wn_spont(:,1,1))./(wn_spont(:,1,2)+wn_spont(:,1,1));
%     subplot(3,3,t+3)
%     hSpont= hist(MI_stat_spont(useSpont),-1:.1:1);
%     Mbins=-1:.1:1
%     bar(Mbins,hSpont/sum(useSpont))
%     xlabel('MI'); ylabel('fraction of cells'); xlim([-1.5 1.5]);ylim([0 .4]);
%     set(gcf,'Name','MI wn spontaneous stationary')
%     title(titles{t});
%     subplot(3,2,2)
%     meanMI_stat_spont(t) = nanmean(MI_stat_spont(useSpont))
%     bar(meanMI_stat_spont);ylim([-1 1]);
% end
% 
% % MI wn SPONTANEOUS MOVING %%% EDIT ALL BASED ON THIS!!
% figure
% for t = 1:8
%     thresh = 1;
%     useSpont = wn_spont(:,2,1)>0 & wn_spont(:,2,2)>0 & (wn_spont(:,2,1)>thresh | wn_spont(:,2,2)>thresh) & treatment'==t;
%     MI_mv_spont = (wn_spont(:,2,2)-wn_spont(:,2,1))./(wn_spont(:,2,2)+wn_spont(:,2,1));
%     subplot(3,3,t+3)
%     hSpont= hist(MI_mv_spont(useSpont),-1:.1:1);
%     Mbins=-1:.1:1
%     bar(Mbins,hSpont/sum(useSpont))
%     xlabel('MI'); ylabel('fraction of cells'); xlim([-1.5 1.5]);ylim([0 .4]);
%     set(gcf,'Name','MI wn spontaneous move')
%     title(titles{t});
%     subplot(3,3,2)
%     meanMI_mv_spont(t) = nanmean(MI_mv_spont(useSpont))
%     bar(meanMI_mv_spont);ylim([-1 1]);
% end

titles = {'Saline','DOI','5HT','Ketanserin', 'Ketanserin + DOI','MGluR2', 'MGluR2 + DOI', 'Lisuride'};
figure
for t=1:8
    useN = find(goodAll & treatment==t)
    for i = 1:length(useN)
        subplot(4,2,t)
        plot(cv2Dark(treatment==t,1), cv2Dark(treatment==t,2),'.');title(titles{t});
        xlabel('Pre CV2 dark');ylabel('Post CV2 dark'); ylim([0 2]);xlim([0 2]); hold on;
        plot([0 2],[0 2])
    end
end



% %%%%%%%% ===================================================== %%%%%%%%%%
% %%% try to decode SF, compare pre and post sessions for each treatment %%%
% 
% drift_trial_psth
% 
% driftSF_psth = [drift_trial_psth(
% 
% 
% for t=1:4
%     useN= find(treatment==t)
% 
%        SFresp = squeeze(drift_sf(:,:,2,1));
%       lowSFpre = squeeze(drift_sf(treatment==t,1,2,1));
%       highSFpre = squeeze(drift_sf(treatment==t,7,2,1));
%     SFresp(treatment==t & drift_sf_trial{1}
%     preSF = [lowSFpre ;highSFpre];
%     squeeze(preSF)
%   SF(1:length(lowSFpre))=1; SF(length(lowSFpre)+1:length(preSF)) =2;
% end
% 
% predata =[preSF SF']
% 
% 
% % % %=====
% for t=1:4
%     useN= find(treatment==t)
%     for i = 1:length(useN)
%         lowSFpost = squeeze(drift_sf(treatment==t,1,2,2));
%         highSFpost = squeeze(drift_sf(treatment==t,7,2,2));
%     end
%     
%     postSF= [lowSFpost ;highSFpost];
%     squeeze(postSF)
%     SFpost(1:length(lowSFpost))=1; SF(length(lowSFpost)+1:length(postSF)) =2;
% end
% % % 
% postdata =[postSF SF']
% PredictorNames = {'column_1'};
% 
% % X = table2array(varfun(@double, postdata(:,trainedClassifier_fineknn.PredictorNames)));
% % yfit = predict(trainedClassifier_fineknn, X)
% 
% yfit = predict(trainedClassifier_fineknn, postdata{:,trainedClassifier_fineknn.PredictorNames})
