clear all
close all

dbstop if error

batchDOIephys_filtered; %%% load batch file
set(groot,'defaultFigureVisible','off') %disable figure plotting
%set(groot,'defaultFigureVisible','on')

%%% select the sessions you want based on filters
use =  find(strcmp({files.notes},'good data'))%useSess = use;
%use =  find( strcmp({files.treatment},'5HT') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}) )

%for specific experiment:
%use =  find(strcmp({files.notes},'good data') & strcmp({files.expt},'022417'))
sprintf('%d selected sessions',length(use))

saline=1; doi=2; ht=3; ketanserin=4; ketandoi=5; mglur2=6; mglur2doi=7; lisuride=8;
% movieFile = 'C:\Users\Angie Michaiel\Desktop\movie files\cortex\wn_cortex_012alpha1_5hzLg30Hz.mat';
%      load(movieFile);
savePDF=0;
redo = 1;
n=0; ncorr=0; %%% number of units loaded, ncorr= number of correlation pairs
for i = 1:length(use)
    % close all
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
        layerInfo = getLayer(clustfile,afile,files(use(i)).tip1,files(use(i)).tip2,files(use(i)).angle, 0); %(needs histo information, but will give layers for all sites)
        layerAll(cellrange) = layerInfo.units;
        layerSites(i,:) = layerInfo.sites;
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
            spd = getSpeed(clustfile,afile,files(use(i)).blockWn{prepost},0);
            speedHistWn(i,:,prepost) = hist(spd.v,0.5:1:100)/length(spd.v);
            speedTrace{i,prepost}=spd.v;
        end
    else
        for prepost = 1:2
            spd = getSpeed(clustfile,afile,files(use(i)).blockDrift{prepost},0);
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
        
         for prepost=1:2
            lfpMoveDrift = getLfpMovement(clustfile,afile,files(use(i)).blockDrift{prepost},0);
            LFPallDrift(i,:,:,prepost) =squeeze(nanmedian(lfpMoveDrift.meanSpect, 1));
            if size(lfpMoveDrift.meanSpect,1)==16
                LFPallChDrift(i,:,:,:,prepost) = lfpMoveDrift.meanSpect;
                LFPfreqDrift(:,:) = lfpMoveDrift.freq;
                display('good lfp');
                lfpMoveDrift
            else
                display('lfp wrong size')
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
            wn_gain(1,cellrange,prepost)=wn.gain;
            wn_frameR(cellrange,:,:,prepost) = downsamplebin(wn.frameR,2,15)/15;
            wn_reliability(1,cellrange,prepost) = wn.reliability;
            
        end
        
        for prepost=1:2
            try
                sta = getSTA(clustfile,afile,files(use(i)).blockWn{prepost},0)
                sta_nx(cellrange,prepost) = sta.nx
                sta_ny(cellrange,prepost)=sta.ny
                sta_sigx(cellrange,prepost)=sta.sigx
                sta_sigy(cellrange,prepost)=sta.sigy
                sta_exp_var(cellrange,prepost)=sta.exp_var
                sta_all_fit(cellrange,prepost)=sta.all_fit;
                sta_all_img(cellrange,prepost)=sta.all_img;
                sta_params(cellrange,prepost)=sta.params;
            catch
            end
        end
        
        
        
        %         figure
        %         subplot(2,2,1); imagesc(squeeze(wn_frameR(cellrange,:,1,1)),[ 0 50]); title('pre stop')
        %         subplot(2,2,2); imagesc(squeeze(wn_frameR(cellrange,:,1,2)),[ 0 50]); title('post stop')
        %         subplot(2,2,3); imagesc(squeeze(wn_frameR(cellrange,:,2,1)),[ 0 50]); title('pre move');
        %         subplot(2,2,4); imagesc(squeeze(wn_frameR(cellrange,:,2,2)),[ 0 50]); title('post move')
        %         set(gcf,'Name',sprintf('%s %s',files(use(i)).expt, files(use(i)).treatment));
        %         drawnow
        
        %% lfp power
        %     %%%(right now averages over all sites, should use layer info)
        for prepost=1:2
            lfpMove = getLfpMovement(clustfile,afile,files(use(i)).blockWn{prepost},0);
            LFPall(i,:,:,prepost) =squeeze(nanmedian(lfpMove.meanSpect, 1));
            if size(lfpMove.meanSpect,1)==16
                LFPallCh(i,:,:,:,prepost) = lfpMove.meanSpect;
                LFPfreq(:,:) = lfpMove.freq;
                display('good lfp');
                lfpMove
            else
                display('lfp wrong size')
            end
        end
     end

for prepost = 1:2
    try
        getSpikes(clustfile,afile,files(use(i)).blockWn{prepost},0);
        getLFPraw(clustfile,afile,files(use(i)).blockWn{prepost},0);
    catch
        display('cant get wn')
    end
    try
        getSpikes(clustfile,afile,files(use(i)).blockDrift{prepost},0);
        getLFPraw(clustfile,afile,files(use(i)).blockDrift{prepost},0);
    catch
        display('cant get drift')
    end
    try
        getSpikes(clustfile,afile,files(use(i)).blockDark{prepost},0);
        getLFPraw(clustfile,afile,files(use(i)).blockDark{prepost},0);
    catch
        display('cant get dark')
    end
    
    
end

      if ~isempty(files(use(i)).blockDark{1}) & ~isempty(files(use(i)).blockDark{2})
      
        
        for prepost=1:2
            lfpMoveDark = getLfpMovement(clustfile,afile,files(use(i)).blockDark{prepost},0);
            LFPallDark(i,:,:,prepost) =squeeze(nanmedian(lfpMoveDark.meanSpect, 1))/nanmedian(lfpMoveDark.meanSpect(:));
            if size(lfpMoveDark.meanSpect,1)==16
                LFPallChDark(i,:,:,:,prepost) = lfpMoveDark.meanSpect;
                  display('good lfp');
              %  lfpMove
            else
                display('lfp wrong size')
            end 
       
        end
      end
   
%     getSorting(clustfile,afile,sprintf('%s %s',files(use(i)).expt,files(use(i)).treatment));
%     drawnow    
    
    %% darkness / correlation analysis
    
    %%% need to keep track of n^2 values for correlations
    corrRange=ncorr+1:ncorr+nc^2;
    corrTreatment(corrRange)=treatment(cellrange(1));
    
    %%% prepost correlation for white noise
    if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2})
        dt = 1;
        [preCorr postCorr cv2 R, eigs] =  prepostDOIdarkness(clustfile,afile,files(use(i)).blockWn,dt,0);
        wnCorr(corrRange,1) = preCorr(:); wnCorr(corrRange,2)=postCorr(:);
    end
    if ~isempty(files(use(i)).blockDrift{1}) & ~isempty(files(use(i)).blockDrift{2})
        dt = 1;
        [preCorr postCorr cv2 R eigs] =  prepostDOIdarkness(clustfile,afile,files(use(i)).blockDrift,dt,0);
        driftCorr(corrRange,1) = preCorr(:); driftCorr(corrRange,2)=postCorr(:);
    end
    
    %%%% prepost correlation in darkness
        if ~isempty(files(use(i)).blockDark{1}) & ~isempty(files(use(i)).blockDark{2})

    dt = 1;
    [preCorr postCorr cv2 R eigs] = prepostDOIdarkness(clustfile,afile,files(use(i)).blockDark,dt,0);
    darkCorr(corrRange,1) = preCorr(:); darkCorr(corrRange,2)=postCorr(:);
    
    cv2Dark(cellrange,:) = cv2;
    meanRdark(cellrange,:) = mean(R,2);
        end
    %%% keep track of cell type for correlations
    corrType1 = zeros(size(preCorr)); corrType2 = corrType1;
    for j= 1:length(inh);
        corrType1(j,:)=inh(j); corrType2(:,j)=inh(j);
    end
    corrType1all(corrRange) = corrType1(:) ; corrType2all(corrRange)= corrType2(:);
    
    n= n+nc;
    ncorr= ncorr+nc^2;
end

%save('compileALL_120316_reselected2')


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
clear max
low_wn = squeeze(min(wn_crf,[],2));
max_wn = squeeze(max(wn_crf,[],2));
amp_wn = max_wn-low_wn
useResp = amp_wn(:,1,1)>0& amp_wn(:,1,2)>0 & amp_wn(:,2,1)>0 & amp_wn(:,2,2)>0;
data_wn = goodAll==1 & useResp';
%data_wnInh = goodAll==1 & useResp' & inhAll;

figure
titles={'Saline','DOI','5HT'};
useSTA = sta_exp_var(:,1) & sta_exp_var(:,2)>.4
useSTA = useSTA & data_wn'% &(layerAll==4)'
for t=1:2 %5ht doesnt show up with exp var smaller than .4
   % mdl_nx = fitlm(sta_nx(useSTA & (treatment==t)',1),sta_nx(useSTA&(treatment==t)',2))
  %  mdl_ny= fitlm(sta_ny(useSTA & (treatment==t)',1),sta_ny(useSTA&(treatment==t)',2))
 %   rsquared_nx(t) = mdl_nx.Rsquared.Ordinary
   % rsquared_ny(t) = mdl_ny.Rsquared.Ordinary
   % subplot(2,2,t)
   subplot(1,2,1)
    plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==t)',2),'.','Markersize',10);
    hold on;
    plot(sta_nx(useSTA & (inhAll==1)'&(treatment==t)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==t)',2),'.r','Markersize',10);
    axis square;xlim([0 .7]); ylim([0 .7]); set(gca,'FontSize',18);
   % text(.02, .63, ['r^2 = ' num2str(rsquared_nx(t))],'FontSize',18)
    mpre=nanmean(mean(sta_nx(useSTA & (treatment==t)',1)))
    mpost=nanmean(mean(sta_nx(useSTA & (treatment==t)',2)))
    plot(mpre,mpost,'pg','Markersize',10);hold on
    title(titles{t}); xlabel('Pre nx');ylabel('Post nx');
    plot([0 1],[0 1]);
   % subplot(2,2,t+2)
   subplot(1,2,2)
    plot(sta_ny(useSTA & (inhAll==0)'& (treatment==t)',1),sta_ny(useSTA &(inhAll==0)'&(treatment==t)',2),'.','Markersize',10);
    hold on;
    plot(sta_ny(useSTA & (inhAll==1)'& (treatment==t)',1),sta_ny(useSTA &(inhAll==1)'&(treatment==t)',2),'.r','Markersize',10);
    axis square;xlim([0 .7]);ylim([0 .7]);set(gca,'FontSize',18);
    hold on; plot([0 1],[0 1]);
   % text(.02, .55, ['r^2 = ' num2str(rsquared_ny(t))],'FontSize',18)
    mpre=nanmean(mean(sta_ny(useSTA  & (treatment==t)',1)))
    mpost=nanmean(mean(sta_ny(useSTA & (treatment==t)',2)))
    plot(mpre,mpost,'pg','Markersize',10)
    title(titles{t},'FontSize',30); xlabel('Pre ny');ylabel('Post ny');
   % n_cells(t) = sum(useSTA &(treatment==t)')
    %text(.02 ,.65, ['n = ' num2str(n_cells(t))],'FontSize',18)
end

figure
set(gcf,'Name', 'nx prepost all layers')
subplot(3,3,1)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==1)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==1)' & (layerAll==2|3)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==1)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==1)' & (layerAll==2|3)',2),'r.','Markersize',10);
title('Layer 2/3');xlabel('pre nx');ylabel('post nx');
subplot(3,3,2)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==1)' & (layerAll==4)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==1)' & (layerAll==4)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==1)' & (layerAll==4)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==1)' & (layerAll==4)',2),'r.','Markersize',10);
title('Layer 4');xlabel('pre nx');ylabel('post nx');
subplot(3,3,3)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==1)' & (layerAll==5)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==1)' & (layerAll==5)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==1)' & (layerAll==5)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==1)' & (layerAll==5)',2),'r.','Markersize',10);
title('Layer 5');xlabel('pre nx');ylabel('post nx');
subplot(3,3,4)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==2)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==2)' & (layerAll==2|3)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==2)' & (layerAll==2|3)',2),'r.','Markersize',10);
xlabel('pre nx');ylabel('post nx');
subplot(3,3,5)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==2)' & (layerAll==4)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==2)' & (layerAll==4)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==4)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==2)' & (layerAll==4)',2),'r.','Markersize',10);
xlabel('pre nx');ylabel('post nx');
subplot(3,3,6)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==2)' & (layerAll==5)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==2)' & (layerAll==5)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==5)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==2)' & (layerAll==5)',2),'r.','Markersize',10);
xlabel('pre nx');ylabel('post nx');
subplot(3,3,8)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==2|3)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==2|3)',2),'r.','Markersize',10);
xlabel('pre nx');ylabel('post nx')
subplot(3,3,8)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==4)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==4)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==4)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==4)',2),'r.','Markersize',10);
xlabel('pre nx');ylabel('post nx');
subplot(3,3,9)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==5)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==5)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==5)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==5)',2),'r.','Markersize',10);
xlabel('pre nx');ylabel('post nx');


figure
set(gcf,'Name', 'ny prepost all layers')
subplot(3,3,1)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==1)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==1)' & (layerAll==2|3)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==1)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==1)' & (layerAll==2|3)',2),'r.','Markersize',10);
title('Layer 2/3');xlabel('pre ny');ylabel('post ny')
subplot(3,3,2)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==1)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==1)' & (layerAll==4)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==1)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==1)' & (layerAll==4)',2),'r.','Markersize',10);
title('Layer 4');xlabel('pre ny');ylabel('post ny')
subplot(3,3,3)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==1)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==1)' & (layerAll==5)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==1)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==1)' & (layerAll==5)',2),'r.','Markersize',10);
title('Layer 5');xlabel('pre ny');ylabel('post ny')
subplot(3,3,4)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==2)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==2)' & (layerAll==2|3)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==2)' & (layerAll==2|3)',2),'r.','Markersize',10);
xlabel('pre ny');ylabel('post ny')
subplot(3,3,5)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==2)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==2)' & (layerAll==4)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==2)' & (layerAll==4)',2),'r.','Markersize',10);
xlabel('pre ny');ylabel('post ny')
subplot(3,3,6)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==2)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==2)' & (layerAll==5)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==2)' & (layerAll==5)',2),'r.','Markersize',10);
xlabel('pre ny');ylabel('post ny')
subplot(3,3,8)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==2|3)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==2|3)',2),'r.','Markersize',10);
xlabel('pre ny');ylabel('post ny')
subplot(3,3,8)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==4)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==4)',2),'r.','Markersize',10);
xlabel('pre ny');ylabel('post ny')
subplot(3,3,9)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==5)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==5)',2),'r.','Markersize',10);
xlabel('pre ny');ylabel('post ny')



useSTA = sta_exp_var(:,1) & sta_exp_var(:,2)>.4
useSTA = useSTA & data_wn'% &(layerAll==4)'
figure
titles={'Saline', 'DOI'};%,'5HT'}
for t=1:2
subplot(2,3,t)
set(gcf,'Name', 'layer 2/3 prepost nx and ny')
%mdl_nxny_pre = fitlm(sta_nx(useSTA & (treatment==t)' &(layerAll==2|3)',1),sta_ny(useSTA&(treatment==t)'&(layerAll==2|3)',1))
%rsquared_nxny_pre(t) = mdl_nxny_pre.Rsquared.Ordinary
%mdl_nxny_post = fitlm(sta_nx(useSTA & (treatment==t)' &(layerAll==2|3)',2),sta_ny(useSTA&(treatment==t)'&(layerAll==2|3)',2))
%rsquared_nxny_post(t) = mdl_nxny_post.Rsquared.Ordinary
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==t)' & (layerAll==2|3)',1),'.','Markersize',10);
hold on;plot([0 1],[0 1]); axis square;xlim([0 1])
%text(.02, .55, ['r^2 = ' num2str(rsquared_nxny_pre(t))],'FontSize',18)
subplot(2,3,t+3)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)' & (layerAll==2|3)',2),sta_ny(useSTA&(inhAll==0)'&(treatment==t)' & (layerAll==2|3)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]);axis square;xlim([0 1])
title(titles{t});
%text(.02, .55, ['r^2 = ' num2str(rsquared_nxny_post(t))],'FontSize',18)
end

figure
titles={'Saline', 'DOI','5HT'}
for t=1:2
subplot(2,3,t)
set(gcf,'Name', 'all layers prepost nx and ny')
%mdl_nxny_pre = fitlm(sta_nx(useSTA & (treatment==t)',1),sta_ny(useSTA&(treatment==t)',1))
%rsquared_nxny_pre(t) = mdl_nxny_pre.Rsquared.Ordinary
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==t)' ,1),'.','Markersize',10);
hold on;plot([0 1],[0 1]); axis square; xlabel('pre nx');ylabel('pre ny');
%text(.02, .55, ['r^2 = ' num2str(rsquared_nxny_pre(t))],'FontSize',18)
plot(sta_nx(useSTA & (inhAll==1)' & (treatment==t)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==t)' ,1),'r.','Markersize',10);
subplot(2,3,t+3)
%mdl_nxny_post = fitlm(sta_nx(useSTA & (treatment==t)',1),sta_ny(useSTA&(treatment==t)',2))
%rsquared_nxny_post(t) = mdl_nxny_post.Rsquared.Ordinary
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)' ,2),sta_ny(useSTA&(inhAll==0)'&(treatment==t)' ,2),'.','Markersize',10);
hold on;plot([0 1],[0 1]);xlim([0 1]);axis square; xlabel('post nx');ylabel('post ny');
plot(sta_nx(useSTA & (inhAll==1)' & (treatment==t)' ,2),sta_ny(useSTA&(inhAll==1)'&(treatment==t)' ,2),'r.','Markersize',10);
%text(.02, .55, ['r^2 = ' num2str(rsquared_nxny_post(t))],'FontSize',18)
title(titles{t});
end

%useSTA = find((sta_exp_var(:,1) & sta_exp_var(:,2)>.6) &(treatment==doi)' &(layerAll==4)');
% for prepost=1:2
%     figure
%       if prepost==1
%         set(gcf,'Name','pre')
%     else
%          set(gcf,'Name','post');end
%     for c=1:length(useSTA)
%         subplot(10,3,c)
%         imagesc(sta_all_img{useSTA(c),prepost});axis square;
%         %imagesc(sta_all_fit{useSTA(c),prepost});axis square
%     end
% end

useSTA = find((sta_exp_var(:,1) & sta_exp_var(:,2)>.4) &(treatment==saline)')% &(layerAll==4)');
for prepost=1:2
    figure
    if prepost==1
        set(gcf,'Name','pre')
    else
        set(gcf,'Name','post');end
    for c=1:length(useSTA)
        subplot(10,3,c)
        %imagesc(sta_all_img{useSTA(c),prepost});axis square
        imagesc(sta_all_fit{useSTA(c),prepost},[-40 40]);axis square; colormap jet
    end
end



titles = {'Saline','DOI','5HT'}
figure
for t=1:3
    subplot(1,3,t)
    set(gcf,'Name','wn evoke prepost')
    plot(wn_evoked(data_wn&treatment==t&~inhAll,1,1),wn_evoked(data_wn&treatment==t&~inhAll,1,2),'.');hold on;axis square;
    xlim([-10 30]); ylim([-10 30]);
    plot([-30 30],[-30 30]);
    plot(wn_evoked(data_wn&treatment==t&inhAll==1,1,1),wn_evoked(data_wn&treatment==t&inhAll==1,1,2),'.r');
    title(titles{t})
end

%ketanserin & ketan doi  more stable than saline????
    titles = {'layer 1','layer 2','layer 3','layer 4','layer 5','layer 6'};
C = {[1 0 0],[.5 0 0]} %bright = pre
D = {[0 1 0],[0 .5 0]}
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
            mn_wn = squeeze(mean(wn_crf(data_wn & layerAll==l & treatment==t,:,1,prepost),1)); %stationary prepost
            mn_wn_mv = squeeze(mean(wn_crf(data_wn & layerAll==l & treatment==t,:,2,prepost),1)); %running prepost
            plot(mn_wn,'Color',C{prepost}) ;hold on;
            plot(mn_wn_mv,'Color', D{prepost});
            title(titles{l});

    end
end
end
useN= data_wn & treatment==saline
MI_crf = (wn_crf(useN,:,:,2)-wn_crf(useN,:,:,1))./(wn_crf(useN,:,:,2)+wn_crf(useN,:,:,1));

CRFsuppress = find(MI_crf <-.5)
CRFup = find(MI_crf >.5)

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
    %subplot(np,np,c);
    subplot(5,4,c)
    % plot(wn_crf(useN(c),:,prepost))
    plot(wn_crf(CRFsuppress(c),:,1,1),'Color',[0.5 0 0]); hold on;  plot(wn_crf(CRFsuppress(c),:,2,1),'Color',[0 0.5 0]);
    plot(wn_crf(CRFsuppress(c),:,1,2),'Color',[1 0 0]);  plot(wn_crf(CRFsuppress(c),:,2,2),'Color',[0 1 0]);
    yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)]);
end


% % figure
% % for c=1:length(data_wn)
% % subplot(5,5,c)
% % imagesc(sta_all_img{c &treatment==doi,1});axis square
% % end
% % figure
% % for c=1:length(data_wn)
% % subplot(5,5,c)
% % imagesc(sta_all_img{c,2});axis square
% % end
% %
% % figure
% % for c=1:length(data_wn)
% % subplot(5,5,c)
% % imagesc(sta_all_fit{c,1});axis square
% % end
% % figure
% % for c=1:length(data_wn)
% % subplot(5,5,c)
% % imagesc(sta_all_fit{c,2});axis square
% % end

% figure
% plot(wn_reliability(1,data & treatment==doi,1)); hold on;
% plot(wn_reliability(1,data & treatment==doi,2))
clear max amp low
%need to subract min from amplitude
peakresp = squeeze(max(drift_orient,[],2));
low = squeeze(min(drift_orient,[],2));
amp = peakresp-low;
useResp = amp(:,1,1)>1 & amp(:,1,2)>1 |amp(:,2,1)>1 & amp(:,2,2)>1;
data = goodAll & useResp' & ~inhAll ;%which cells to include
dataInh =goodAll & useResp' & inhAll ;

for t=1:3
    figure
    if t == 1
        set(gcf,'Name', 'hist saline');
    elseif t==2
        set(gcf,'Name', 'mi hist doi');
    else
        set(gcf,'Name', 'mi hist 5ht');
    end
    useN=goodAll & treatment==t
    
    %miDrift= (drift_cond_tcourse(useN,1,2,:,:,:)-drift_cond_tcourse(useN,1,1,:,:,:))./(drift_cond_tcourse(useN,1,2,:,:,:)+drift_cond_tcourse(useN,1,1,:,:,:));
    %miDrift= (drift_tcourse(useN,:,1,2,:)-drift_tcourse(useN,:,1,1,:)./(drift_tcourse(useN,:,1,2,:))+drift_tcourse(useN,:,1,1,:));
    miDrift= (mnPsth(useN &(layerAll==2|layerAll==3) ,2,:)-mnPsth(useN&(layerAll==2|layerAll==3),1,:))./((mnPsth(useN&(layerAll==2|layerAll==3),2,:)+mnPsth(useN&(layerAll==2|layerAll==3),1,:)));
    h= hist(miDrift,-1:.1:1);
    Mbins=-1:.1:1
    subplot(1,3,1)
    bar(Mbins,h/sum(useN));ylim([0 .12]);axis square
    for i=4:5
        miDrift= (mnPsth(useN &layerAll==i ,2,:)-mnPsth(useN&layerAll==i,1,:))./((mnPsth(useN&layerAll==i,2,:)+mnPsth(useN&layerAll==i,1,:)));
        % miDrift= (drift_trial_psth(useN,:,2)-drift_trial_psth(useN,:,1))./(drift_trial_psth(useN,:,2)+drift_trial_psth(useN,:,1));
        h= hist(miDrift,-1:.1:1);
        subplot(1,3,i-2)
        Mbins=-1:.1:1
        bar(Mbins,h/sum(useN));ylim([0 .12]); axis square
        
        driftSuppress = miDrift <-.4;
        driftUp = miDrift >.4;
    end
    % h= hist(miDrift(driftSuppress),-1:.1:1)
    % subplot(2,2,3)
    % bar(Mbins,h/sum(useN));ylim([0 4]);xlim([-1.5 1.5]); title('suppressed');
    %
    % h= hist(miDrift(driftUp),-1:.1:1);
    % subplot(2,2,4)
    % bar(Mbins,h/sum(useN));ylim([0 4]);xlim([-1.5 1.5]);title('facilitated');
end

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
imagesc(mnPsth(data&treatment==doi,1)); %%% show full dataset (n x t)
subplot(1,2,2)
imagesc(mnPsth(data&treatment==doi,2))

clear max c h g prefOri
useN =data & treatment==doi
%for prepost=1:2

%%% plot all unit responses
usespks = find(goodAll);
for i = 1:length(usespks);
    if mod(i,48)==1
        figure
    end
    subplot(6,8,mod(i-1,48)+1);
    hold on
    plot(drift_orient(usespks(i),:,1,1),'Color',[0.5 0 0]);
    plot(drift_orient(usespks(i),:,1,2),'Color',[1 0 0]);
    plot(drift_orient(usespks(i),:,2,1),'Color',[0 0.5 0]);
    plot(drift_orient(usespks(i),:,2,2),'Color',[0 1 0]);
    xlim([1 12]); yl = get(gca,'Ylim'); ylim([0 max(5,yl(2))]);
end

%drift_cond_tcourse(cellrange,mv,prepost,ori,sf,:)
titles = {'pre stop','post stop','pre move','post move'};
for c = 1:3
    figure
    if c==1
        select = (~inhAll & layerAll<5); set(gcf,'Name','layer 2-4');
    elseif c==2
        select = (~inhAll & layerAll==5); set(gcf,'Name','layer 5');
    elseif c==3
        select = (inhAll); set(gcf,'Name','inh');
    end
    
    
    for prepost=1:2
        clear ori_data
        for ori = 1:6
            ori_data(:,:,ori,:) = squeeze(nanmean(nanmean(drift_cond_tcourse(:,:,prepost,[ori ori+6],:,:),5),4));
        end
        for n= 1:size(ori_data,1);
            [y peak] = max(squeeze(mean(mean(ori_data(n,:,:,5:30),4),2)));
            ori_data(n,:,:,:) = circshift(ori_data(n,:,:,:),-peak,3);
        end
        
        ori_data =downsamplebin(ori_data,4,2,1);
        ori_data = circshift(ori_data,5,4);
        
        for mv = 1:2
            subplot(2,2,prepost+2*(mv-1));
            plot(squeeze(nanmean(ori_data(find(goodAll & select),mv,:,:),1))');
            title(titles{prepost+2*(mv-1)}); if c<=2; ylim([0 6]); else ylim([0 12]); end
        end
    end
end



mv = 1
clear ori_data
for ori = 1:6
    ori_data(:,:,ori,:) = squeeze(nanmean(nanmean(drift_cond_tcourse(:,mv,:,[ori ori+6],:,:),5),4));
end
for n= 1:size(ori_data,1);
    [y peak] = max(squeeze(mean(mean(ori_data(n,:,:,5:30),4),2)));
    ori_data(n,:,:,:) = circshift(ori_data(n,:,:,:),-peak,3);
end

ori_data =downsamplebin(ori_data,4,2,1);
ori_data = circshift(ori_data,5,4);

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'prepost', 'orientation', 'Condition-independent', 'treatment/stim Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

d = ori_data(:,:);
clean = ~isnan(mean(d,2));
data = ori_data(find(clean & goodAll'),:,:,:);
data = permute(data,[1 3 2 4]);

d = permute(data,[1 4 2 3]);
d = d(:,:);

celltype = layerAll; celltype(inhAll) = 7;
[y order] = sort(celltype(find(clean & goodAll')));

resp = d(:,:);
for n= 1:size(resp,1);
    m = max(resp(n,:)); m= max(5,m); resp(n,:) = resp(n,:)/m;  %%% normalized
    data(n,:,:,:) = data(n,:,:,:)/m;
end






%%

for c= 1:length(useN)
    for prepost =1:2
        [resp oind] = max(drift_orient(useN,:,1,prepost)');
        [g h]=(max(resp(:)));
        prefOri(prepost,c)=squeeze(oind(h));
    end
end

%
% size(prefOri')
%
% prefOri_pre = prefOri(1,:)'
% prefOri_post = prefOri(2,:)'
%
% % figure
% % scatterhist(prefOri_pre,prefOri_post,12)
%
% figure
% subplot(1,2,1)
% hist(prefOri_pre);ylim([0 180]);title('hist pref orientation pre');
% subplot(1,2,2)
% hist(prefOri_post);ylim([0 180]);title('post');
%
% prefOri = [prefOri_pre'; prefOri_post']'


%drift_psth= cellfun(@mean, drift_trial_psth(data,1,:),'UniformOutput', false)
% data = goodAll & useResp' & ~inhAll
% dt = 0.05;
% titles = {'.01','.02','.04','.08','.16','.32'};
% figure
% set(gcf,'Name', 'prepost response to 6 SFs')
% for sf = 1:6
%     subplot(2,3,sf)
%     mn_pre=(squeeze(mean(drift_cond_tcourse(data & treatment==doi,1,1,3,sf,:),1))); % cellrange,mv,prepost,ori,sf,:
%     mn_pre = mn_pre - repmat(mn_pre(1,:),[50 1]);
%     mn_pre = circshift(mn_pre,10)
%     mn_post =(squeeze(mean(drift_cond_tcourse(data & treatment==doi,1,2,3,sf,:),1)));
%     mn_post = mn_post - repmat(mn_post(1,:),[50 1]);
%     mn_post = circshift(mn_post,10)
%     plot((1:length(mn_pre)-5)*dt -dt/2,mn_pre(1:45,:));xlabel('secs'); ylabel('sp/sec');
%     title(titles{sf});
%     hold on;
%     plot((1:length(mn_post)-5)*dt -dt/2,mn_post(1:45,:));
% end

%titles = {};
dt = 0.05;
figure
set(gcf,'Name', 'prepost response to 6 orientations')
for ori=1:6
    subplot(2,3,ori)
    mn_pre=(squeeze(mean(ori_data(data &treatment==ht,1,ori,:),1))); % cellrange,mv,prepost,ori,sf,:
    mn_pre = mn_pre - repmat(mn_pre(1,:),[25 1]);
    mn_pre = circshift(mn_pre,10)
    mn_post=(squeeze(mean(ori_data(data &treatment==ht,2,ori,:),1)));
    mn_post = mn_post - repmat(mn_post(1,:),[25 1]);
    mn_post = circshift(mn_post,10)
    plot((1:length(mn_pre)-2.5)*dt -dt/2,mn_pre(1:22.5,:));xlabel('secs'); ylabel('sp/sec');
    % title(titles{ori});
    hold on;
    plot((1:length(mn_post)-2.5)*dt -dt/2,mn_post(1:22.5,:)); hold on
end

% titles = {'0','30','60','90','120','150','180','210','240','270','300','330'};
% dt = 0.05;
% figure
% set(gcf,'Name', 'prepost response to 6 orientations')
% for ori = 1:12
%     subplot(3,4,ori)
%     mn_pre=(squeeze(mean(drift_cond_tcourse(data & treatment==ht,1,1,ori,3,:),1))); % cellrange,mv,prepost,ori,sf,:
%     mn_pre = mn_pre - repmat(mn_pre(1,:),[50 1]);
%     mn_pre = circshift(mn_pre,10)
%     mn_post =(squeeze(mean(drift_cond_tcourse(data & treatment==ht,1,2,ori,3,:),1)));
%     mn_post = mn_post - repmat(mn_post(1,:),[50 1]);
%     mn_post = circshift(mn_post,10)
%     plot((1:length(mn_pre)-5)*dt -dt/2,mn_pre(1:45,:));xlabel('secs'); ylabel('sp/sec');
%     title(titles{ori});
%     hold on;
%     plot((1:length(mn_post)-5)*dt -dt/2,mn_post(1:45,:));
% end



% ori=1:12
% for x=1:12
% test = drift_cond_tcourse(data&treatment==doi,1,:,ori==prefOri,:,:);
%
% end



% figure
% plot(drift_ot_tune(data,:,1,1),'b'); hold on;
% plot(drift_ot_tune(data,:,1,2),'r');

% data = reshape(data,size(data,1),size(data,2))
% 
% dist = pdist(squeeze(wn_evoked(data_wn &treatment==ht,1,:)),'correlation');  %%% sort based on correlation coefficient...want to do this for just one orientation or SF
% %dist = pdist(squeeze(drift_cond_tcourse(data & treatment==ht,1,1,3,3,:)),'correlation');
% display('doing cluster')
% tic, Z = linkage(dist,'ward'); toc
% figure
% subplot(3,4,[1 5 9 ])
% display('doing dendrogram')
% [h t perm] = dendrogram(Z,0,'Orientation','Left','ColorThreshold' ,5);
% axis off
% subplot(3,4,[2 3 4 6 7 8 10 11 12 ]);
% imagesc(wn_evoked(perm,:)); axis xy   %%% show sorted data


% [Y e] = mdscale(dist,1);
% [y sortind] = sort(Y);
% figure
% imagesc(mnPsth(sortind,:),[0 1])


%%% plot mean timecourse across all drift stim and cells, by layer / cell-type
%%% this should be updated to select active cells and separate
%%% moving/stationary. Also maybe choose optimal stim for each cell?


% useResp = amp(:,1,1)>2& amp(:,1,2)>2 | amp(:,2,1)>2 & amp(:,2,2)>2;
% useResp= amp(:,1,1)>0
% data = goodAll & useResp' & ~inhAll ;%which cells to include
% dataInh =goodAll & useResp' & inhAll ;

%data2 = reshape(data,size(data,1),size(data,2),size(data,3))
titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
dt = 0.05;
figure
for i= 1:5
    figure
    for t=1:3
        if i==1
            set(gcf,'Name','grating lyr5');
            mn = squeeze(mean(mnPsth(data & layerAll==5 & treatment==t,:,:),1))';
            ylim([-.5 3]); xlim([0 2.5]);
        elseif i==2
            set(gcf,'Name','grating lyr 2/3');
            mn = squeeze(mean(mnPsth(data & (layerAll==2 |layerAll==3) & treatment==t,:,:),1))';
            ylim([-.5 5]);xlim([0 2.5]);
        elseif i ==3
            set(gcf,'Name','grating lyr 4');
            mn = squeeze(mean(mnPsth(data & layerAll==4 & treatment==t,:,:),1))';
            ylim([-.5 4]);xlim([0 2.5]);
        elseif i==4
            set(gcf,'Name','grating inh');
            mn = squeeze(mean(mnPsth(dataInh & treatment==t,:,:),1))';
            ylim([-.5 12]);xlim([0 2.5]);
        elseif i==5
            set(gcf,'Name','grating lyr6');
            mn = squeeze(mean(mnPsth(data & layerAll==6 & treatment==t,:,:),1))';
        end
        mn = mn - repmat(mn(1,:),[50 1]);
        mn = circshift(mn,10);
        subplot(1,3,t)
        plot((1:length(mn)-5)*dt -dt/2,mn(1:45,:),'LineWidth',2);axis square; set(gca,'fontsize', 18);
        %title(titles{t}); xlabel('time (s)');
        %ylabel('spikes/sec');
        %ylim([0 max(mn(:))+1])
    end
end

%need to make scatter of transient & sustained responses/cell/treatment pre and post
for i= 1:4
    figure
    for t=1:3
        if i==1
            set(gcf,'Name','grating lyr5');
            trans = mean(mnPsth(data & layerAll==5 & treatment==t,:,1:10));
            % ylim([0 1]); xlim([0 2.5]);
        elseif i==2
            set(gcf,'Name','grating lyr 2/3');
            trans = mean(mnPsth(data & (layerAll==2| layerAll==3) & treatment==t,:,1:10))
            %ylim([-.5 5]);xlim([0 2.5]);
        elseif i ==3
            set(gcf,'Name','grating lyr 4');
            trans = mean(mnPsth(data & (layerAll==4) & treatment==t,:,1:10))
            % ylim([-.5 4]);xlim([0 2.5]);
        else i==4
            set(gcf,'Name','grating inh');
            trans = mean(mnPsth(dataInh & treatment==t,:,1:10))
            %ylim([-.5 12]);xlim([0 2.5]);
            %         mn = mn - repmat(mn(1,:),[50 1]);
            %         mn = circshift(mn,10);
            subplot(1,3,t)
            plot(trans(i,1,:),trans(i,2,:),'.');axis square
        end
    end
end


% titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
% dt = 0.05;
% for i= 1
%     figure
%     for t=1:3
%         if i==1
%             set(gcf,'Name','grating lyr5');
%             mn = squeeze(nanmean(ori_data(data & layerAll==5 & treatment==t,:,:,:),1));
%             ylim([-.5 3]); xlim([0 2.5]);
%         elseif i==2
%             set(gcf,'Name','grating lyr 2/3');
%             mn = squeeze(mean(ori_data(data & (layerAll==2 |layerAll==3) & treatment==t,:,:,:),1))';
%             ylim([-.5 5]);xlim([0 2.5]);
%         elseif i ==3
%             set(gcf,'Name','grating lyr 4');
%             mn = squeeze(mean(ori_data(data & layerAll==4 & treatment==t,:,:,:),1))';
%             ylim([-.5 4]);xlim([0 2.5]);
%         elseif i==4
%             set(gcf,'Name','grating inh');
%             mn = squeeze(mean(ori_data(dataInh & treatment==t,:,:,:),1))';
%             ylim([-.5 12]);xlim([0 2.5]);
%         elseif i==5
%             set(gcf,'Name','grating lyr6');
%             mn = squeeze(mean(ori_data(data & layerAll==6 & treatment==t,:,:,:),1))';
%         end
%        %mn = mn - repmat(mn(1,:,:),[25 1]);
%        mn = circshift(mn,10);
%        % subplot(1,3,t)
%        plot((1:length(mn)-5)*dt -dt/2,mn(:,1:20),'LineWidth',2);%axis square; set(gca,'fontsize', 18);
%        % title(titles{t}); xlabel('time (s)');
%        % ylabel('spikes/sec');
%         %ylim([0 max(mn(:))+1])
%     end
% end




titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
figure
for t=1:3
    mnpre = squeeze(mean(mnPsth(data & layerAll==5 & treatment==t,1,:),1))';
    mnpost = squeeze(mean(mnPsth(data & layerAll==5 &treatment==t,2,:),1))';
    %mnpre = mnpre - repmat(mnpre,[1 50]);
    mnpre = circshift(mnpre,10);
    %mnpost = mnpost - repmat(mnpost(1,:),[50 1]);
    scalefact(t)=mean(mnpre)./mean(mnpost)
    mnpost_scale = scalefact(t)*(circshift(mnpost,10));
    mnpost=circshift(mnpost,10);
    subplot(1,3,t)
    plot((1:length(mnpre)-5)*dt -dt/2,mnpre(:,1:45));axis square; title(titles{t}); xlabel('secs'); ylabel('sp/sec'); hold on;
    plot((1:length(mnpost)-5)*dt -dt/2,mnpost(:,1:45),'r')
    plot((1:length(mnpost_scale)-5)*dt -dt/2,mnpost_scale(:,1:45),'g')
    text(.02, .63, ['scale factor = ' num2str(scalefact(t))],'FontSize',12)

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


figure
useN= find(data & (layerAll==2|3) &treatment==doi)
length(useN)
for c=1:length(useN)
    set (gcf,'Name','Lyr 2/3 psth/unit')
    subplot(11,11,c)
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
useN= find(data & layerAll==5 & treatment==ht)
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
useInh= find(dataInh &treatment==saline)
for c=1:length(useInh)
    set (gcf,'Name','inh psth/unit')
    subplot(10,10,c)
    plot(squeeze(mnPsth(useInh(c),1,:)));axis square;hold on;
    plot(squeeze(mnPsth(useInh(c),2,:)),'r')
end

figure
useInh= find(dataInh & layerAll==5 & treatment==ht)
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

for t =1:3
    preOsi = drift_osi(data &treatment==t,1,1)
    postOsi = drift_osi(data &treatment==t,1,2)
    preOsi(isnan(preOsi)) = [];
    postOsi(isnan(postOsi)) = [];
    preOsi=abs(preOsi)
    postOsi=abs(postOsi)
%     mdl = fitlm(preOsi(1:length(postOsi)),postOsi)
%     rsquared(t) = mdl.Rsquared.Ordinary
end

titles = {'Saline','DOI','5HT'};% ,'Ketanserin', 'Ketanserin + DOI', 'MGluR2','MGlur2 + DOI', 'Lisuride'};
figure
for t=1:3
    subplot(1,3,t)
    plot(abs(drift_osi(treatment==t,1,1)),abs(drift_osi(treatment==t,1,2)),'.')
    set(gca,'FontSize',20)
    hold on; axis square; xlim([0 1.1])
    plot([-1.5 1.5], [-1.25 1.25]);
    title(titles{t})
    %text(-.25, .85, ['r^2 = ' num2str(rsquared(t))])
end


% meanOSI = [mean(drift_osi(data,1,1))  mean(drift_osi(data,1,2))];
% figure
% hb = bar(meanOSI)
% plot(hb)

figure
for l=1:6
    osi(2,:) = squeeze(nanmean(drift_osi(data & layerAll==l,1,2),1))
    osierr(2,:) = squeeze(nanstd(drift_osi(data & layerAll==l,1,2),1))/sqrt(sum(drift_osi(data& layerAll==l,1,2)));
    osi(1,:) = squeeze(nanmean(drift_osi(data & layerAll==l,1,1),1))
    osierr(1,:) = squeeze(nanstd(drift_osi(data & layerAll==l,1,1),1))/sqrt(sum(drift_osi(data & layerAll==l,1,1)));
    subplot(2,3,l)
    %barweb(osi,osierr)
    bar(osi);hold on;
end

figure
subplot(2,2,1); imagesc(squeeze(wn_frameR(find(data_wn),:,1,1)),[ 0 50]); title('pre stop')
subplot(2,2,2); imagesc(squeeze(wn_frameR(find(data_wn),:,1,2)),[ 0 50]); title('post stop')
subplot(2,2,3); imagesc(squeeze(wn_frameR(find(data_wn),:,2,1)),[ 0 50]); title('pre move');
subplot(2,2,4); imagesc(squeeze(wn_frameR(find(data_wn),:,2,2)),[ 0 50]); title('post move')
set(gcf,'Name',sprintf('%s %s',files(use(i)).expt, files(use(i)).treatment));
drawnow


figure
subplot(2,2,1); plot(squeeze(nanmean(wn_frameR(find(data_wn),:,1,1),1))); title('pre stop'); ylim([0 5])
subplot(2,2,2); plot(squeeze(nanmean(wn_frameR(find(data_wn),:,1,2),1))); title('post stop');ylim([0 5])
subplot(2,2,3); plot(squeeze(nanmean(wn_frameR(find(data_wn),:,2,1),1))); title('pre move'); ylim([0 5])
subplot(2,2,4); plot(squeeze(nanmean(wn_frameR(find(data_wn),:,2,2),1))); title('post move'); ylim([0 5])
set(gcf,'Name',sprintf('%s %s',files(use(i)).expt, files(use(i)).treatment));
drawnow

clear cycR
for f = 1:20;
    cycR(:,f,:,:) = nanmean(wn_frameR(:,f:20:end,:,:),2);
end

figure
subplot(2,2,1); imagesc(squeeze(cycR(find(data_wn),:,1,1)),[ 0 50]); title('pre stop')
subplot(2,2,2); imagesc(squeeze(cycR(find(data_wn),:,1,2)),[ 0 50]); title('post stop')
subplot(2,2,3); imagesc(squeeze(cycR(find(data_wn),:,2,1)),[ 0 50]); title('pre move');
subplot(2,2,4); imagesc(squeeze(cycR(find(data_wn),:,2,2)),[ 0 50]); title('post move')

figure
subplot(2,2,1); plot(squeeze(mean(cycR(find(data_wn),:,1,1),1))); title('pre stop'); ylim([1.5 4])
subplot(2,2,2); plot(squeeze(mean(cycR(find(data_wn),:,1,2),1))); title('post stop'); ylim([1.5 4])
subplot(2,2,3); plot(squeeze(mean(cycR(find(data_wn),:,2,1),1))); title('pre move'); ylim([1.5 4])
subplot(2,2,4); plot(squeeze(mean(cycR(find(data_wn),:,2,2),1))); title('post move'); ylim([1.5 4])

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
    plot(spont(find(data & treatment==t),1,1),spont(find(data& treatment==t),1,2),'.'); title(titles{t}); hold on; axis square; plot([0 50],[0 50]);
    plot(spont(find(dataInh&inhAll& treatment==t),1,1),spont(find(dataInh & inhAll& treatment==t),1,2),'.r'); hold on; plot([0 50],[0 50]); axis([0 10 0 10])
end

figure
for t=1:8
    set(gcf,'Name', 'spont pre vs post mv')
    subplot(2,4,t)
    plot(spont(find(data& treatment==t),2,1),spont(find(data& treatment==t),2,2),'.'); title(titles{t}); hold on; plot([0 50],[0 50]); axis square
    plot(spont(find(dataInh&inhAll& treatment==t),2,1),spont(find(dataInh & inhAll& treatment==t),2,2),'.r'); hold on; plot([0 50],[0 50]); axis([0 15 0 15])
end

figure
for t=1:8
    set(gcf,'Name', 'evoked pre vs post stop')
    subplot(2,4,t)
    plot(evoked(find(data& treatment==t),1,1),evoked(find(data& treatment==t),1,2),'.'); title(titles{t}); hold on; plot([0 50],[0 50]);
    plot(evoked(find(dataInh&inhAll& treatment==t),1,1),evoked(find(dataInh & inhAll& treatment==t),1,2),'.r'); hold on; plot([0 50],[0 50]); axis([-10 10 -10 10])
end

figure
for t=1:8
    set(gcf,'Name', 'evoked pre vs post move')
    subplot(2,4,t)
    plot(evoked(find(goodAll& treatment==t),2,1),evoked(find(goodAll& treatment==t),2,2),'.'); title(titles{t}); hold on; plot([0 50],[0 50]);
    plot(evoked(find(goodAll&inhAll& treatment==t),2,1),evoked(find(goodAll & inhAll& treatment==t),2,2),'.r'); hold on; plot([0 50],[0 50]); axis([-10 10 -10 10])
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
    plot(spont(find(goodAll& treatment==t),1,2),spont(find(goodAll& treatment==t),2,2),'.'); title(titles{t}); hold on; plot([0 50],[0 50]);
    plot(spont(find(goodAll&inhAll& treatment==t),1,2),spont(find(goodAll & inhAll& treatment==t),2,2),'.r'); hold on; plot([0 50],[0 50]); axis([0 50 0 50])
end

figure
for t=1:8
    set(gcf,'Name', 'evoked stop vs move post')
    subplot(2,4,t)
    plot(evoked(find(goodAll& treatment==t),1,2),evoked(find(goodAll& treatment==t),2,2),'.'); title(titles{t}); hold on; plot([0 50],[0 50]);
    plot(evoked(find(goodAll&inhAll& treatment==t),1,2),evoked(find(goodAll & inhAll& treatment==t),2,2),'.r');hold on; plot([0 50],[0 50]); axis([0 50 0 50])
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
    pre = [cycR(find(data_wn & ~inhAll & treatment==t),:,1,1) cycR(find(data_wn & ~inhAll& treatment==t),:,2,1)];
    post = [cycR(find(data_wn & ~inhAll& treatment==t),:,1,2) cycR(find(data_wn & ~inhAll& treatment==t),:,2,2)];
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
figure
for t=1:3
    preWnCorr = wnCorr(corrTreatment==t,1)
    postWnCorr=wnCorr(corrTreatment==t,2)
    plot(preWnCorr, postWnCorr,'.'); hold on; axis equal;plot(lsline);
    mdl = fitlm(preWnCorr,postWnCorr)
    rsquared(t) = mdl.Rsquared.Adjusted
    
end

titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI','MGluR2', 'MGluR2 + DOI','Lisuride'};
figure
for i = 1:3
    subplot(1,3,i);
    plot(wnCorr(corrTreatment==i,1),wnCorr(corrTreatment==i,2),'.'); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
    xlabel('pre wn corr'); ylabel('post')
    set(gcf,'Name','Wn Corr')
    text(-.25, .85, ['r^2 = ' num2str(rsquared(i))])
    
end

titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI','MGluR2', 'MGluR2 + DOI','Lisuride'};
figure
for i = 1:3
    h= histc(wnCorr(corrTreatment==i,1),-1:.1:1);
    h2=histc(wnCorr(corrTreatment==i,2),-1:.1:1);
    subplot(3,1,i);
    plot(h,'-b')
    hold on;
    plot(h2,'-r')
    %xlim([-1 1]);
    title(titles{i});
    set(gcf,'Name','Wn Corr Hist')
end

clear rsquared
for t=1:3
    preDriftCorr = driftCorr(corrTreatment==t,1)
    postDriftCorr=driftCorr(corrTreatment==t,2)
    
%     mdl = fitlm(preDriftCorr,postDriftCorr)
%     rsquared(t) = mdl.Rsquared.Adjusted
end

titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI','MGluR2', 'MGluR2 + DOI','Lisuride'};
figure
for i = 1:3
    subplot(1,3,i);
    plot(driftCorr(corrTreatment==i,1),driftCorr(corrTreatment==i,2),'.'); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
    xlabel('pre drift corr'); ylabel('post')
  %  text(-.25, .85, ['r^2 = ' num2str(rsquared(i))])
    set(gcf,'Name','Drift Corr')
end



titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI', 'MGlur2', 'MGlur2 + DOI', 'Lisuride'};
figure
for i = 1:7
    subplot(2,4,i);
    wnCorrHist= myHist2(wnCorr(corrTreatment==i,1),wnCorr(corrTreatment==i,2),-.5:.1:1.5,-.5:.1:1.5);
    %wnCorrHist_pre= myHist2(wnCorr(corrTreatment==i,1),-.5:.1:1.5,-.5:.1:1.5)
    plot(wnCorrHist);hold on; axis square;title(titles{i})
    set(gcf,'Name','Wn Corr'); %ylim([0 6000])
end


%%% plot correlation for darkness
clear rsquared
for t=1:3
    preDarkCorr = darkCorr(corrTreatment==t,1)
    postDarkCorr=darkCorr(corrTreatment==t,2)
    
%     mdl = fitlm(preDarkCorr,postDarkCorr)
%     rsquared(t) = mdl.Rsquared.Adjusted
end
titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI','MGlur2','MGlur2 + DOI', 'Lisuride'};
figure
for i = 1:3
    subplot(1,3,i);
    plot(darkCorr(corrTreatment==i,1),darkCorr(corrTreatment==i,2),'.'); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
    xlabel('pre dark corr'); ylabel('post')
%     text(-.25, .85, ['r^2 = ' num2str(rsquared(i))])
    set(gcf,'Name','Dark Corr')
end

preDarkCorr = darkCorr(corrTreatment==3,1)
postDarkCorr=darkCorr(corrTreatment==3,2)

mdl = fitlm(preDarkCorr,postDarkCorr)
mdl.Rsquared.Adjusted

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

titles ={'Saline', 'DOI','5HT'};
for mv =1:2
    figure
    if mv==1
        set(gcf,'Name','stationary')
    else
        set(gcf,'Name','mv')
    end
    for t=1:3
        subplot(1,3,t)
        plot(drift_spont(goodAll==1 &treatment==t &inhAll==0,mv,1),drift_spont(goodAll==1&treatment==t&inhAll==0,mv,2),'.');
        ylim([0 20]); xlim([0 20]); xlabel('Pre spikes/sec');ylabel('Post spikes/sec');
        axis square
        hold on;
        plot([0 20], [0 20])
        plot(drift_spont(goodAll==1 & inhAll==1 &treatment==t,mv,1),drift_spont(goodAll==1& inhAll==1&treatment==t,mv,2),'r.');
        title(titles{t})
    end
end
clear mpre mpost

titles ={'Saline', 'DOI','5HT'};
for mv =1:2
    figure
    if mv==1
        set(gcf,'Name','stat')
    else
        set(gcf,'Name','mv')
    end
    for t=1:3
        subplot(1,3,t)
        plot(mean(drift_orient(goodAll==1 & treatment==t,:,mv,1),2),mean(drift_orient(goodAll==1&treatment==t,:,mv,2),2),'.','Markersize',10);
        hold on;
        set(gca,'FontSize',18);
        axis square; xlim([0 12]);ylim([0 12]); %title(titles{t}); xlabel('Pre spikes/sec');ylabel('Post spikes/sec');
        
        mpre=nanmean(mean(drift_orient(goodAll==1 & treatment==t,:,mv,1),2))
        mpost=nanmean(mean(drift_orient(goodAll==1 & treatment==t,:,mv,2),2))
        
        allPre = drift_orient(goodAll==1 & treatment==t,mv,1,1)
        allPost = drift_orient(goodAll==1 & treatment==t,mv,1,2)
        
%         mdl_drift = fitlm(allPre,allPost)
%         rsquared_drift(t) = mdl_drift.Rsquared.Ordinary
%         
        plot(mpre,mpost,'+k','Markersize',10,'Linewidth',2)
        plot([0 35], [0 35])
        plot(mean(drift_orient(goodAll==1 & inhAll==1  &treatment==t,:,mv,1),2),mean(drift_orient(goodAll==1 & inhAll==1 &treatment==t,:,mv,2),2),'r.','Markersize',10);
        n_cells(t) = sum(goodAll==1 &treatment==t)
        text(.5, 11.5, ['n = ' num2str(n_cells(t))],'FontSize',18)
%         text(.5, 10.75, ['r^2 = ' num2str(rsquared_drift(t))],'FontSize',18)
    end
end

for t=1:3
    figure
    if t == 1
        set(gcf,'Name', 'hist saline');
    elseif t==2
        set(gcf,'Name', 'mi hist doi');
    else
        set(gcf,'Name', 'mi hist 5ht');
    end
    useN=goodAll & treatment==t
    miDrift= (mean(drift_orient(useN &(layerAll==2|layerAll==3),:,1,2),2)-mean(drift_orient(useN&(layerAll==2|layerAll==3),:,1,1),2))./(mean(drift_orient(useN&(layerAll==2|layerAll==3),:,1,2),2)+mean(drift_orient(useN&(layerAll==2|layerAll==3),:,1,1),2));
    h= hist(miDrift,-1:.2:1);
    Mbins=-1:.2:1
    subplot(1,3,1)
    bar(Mbins,h/sum(useN &(layerAll==2|layerAll==3)),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .4]); xlim([-1.5 1.5]);axis square
    for i=4:5
        miDrift= (mean(drift_orient(useN &(layerAll==i),:,1,2),2)-mean(drift_orient(useN&(layerAll==i),:,1,1),2))./(mean(drift_orient(useN&(layerAll==i),:,1,2),2)+mean(drift_orient(useN&(layerAll==i),:,1,1),2));
        h= hist(miDrift,-1:.2:1);
        subplot(1,3,i-2)
        Mbins=-1:.2:1
        bar(Mbins,h/sum(useN&layerAll==i),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .4]);xlim([-1.5 1.5]); axis square
    end
end

%%% scatter plot of drift spont
for mv = 1:2
    figure
    for i = 1:6
        subplot(2,3,i)
        plot(drift_spont(data & treatment==saline & layerAll ==i,mv,1),drift_spont(data &treatment==saline& layerAll ==i,mv,2),'k.');
        hold on
        plot(drift_spont(data & treatment==doi& layerAll ==i,mv,1),drift_spont(data &treatment==doi& layerAll ==i,mv,2),'r.');
        plot(drift_spont(data & treatment==ht& layerAll ==i,mv,1),drift_spont(data &treatment==ht& layerAll ==i,mv,2),'w.');
        plot(drift_spont(data &treatment==ketanserin& layerAll ==i,mv,1),drift_spont(data &treatment==ketanserin& layerAll ==i,mv,2),'m.');
        plot(drift_spont(data &treatment==ketandoi& layerAll ==i,mv,1),drift_spont(data &treatment==ketandoi& layerAll ==i,mv,2),'c.');
        plot(drift_spont(data &treatment==mglur2& layerAll ==i,mv,1),drift_spont(data &treatment==mglur2& layerAll ==i,mv,2),'g.');
        plot(drift_spont(data &treatment==mglur2doi& layerAll ==i,mv,1),drift_spont(data &treatment==mglur2doi& layerAll ==i,mv,2),'y.');
        plot(drift_spont(data &treatment==lisuride& layerAll ==i,mv,1),drift_spont(data &treatment==lisuride& layerAll ==i,mv,2),'b.');
        plot([0 10],[0 10]); axis equal
        title(sprintf('layer %d',i)); ylabel('post');
        if mv ==1 , xlabel('stop drift spont'), set(gcf, 'Name', 'prepost stationary drift spont');
        else  xlabel('move drift spont'),set(gcf, 'Name', 'prepost move drift spont'); end
        
    end
end

% 
% figure
% for prepost=1:2
% subplot(1,2,prepost)
% imagesc(squeeze(find(LFPallCh(sessionTreatment==2,16,:,1,prepost))))
% end

salineSessions = (sessionTreatment==1)
doiSessions = (sessionTreatment==2)
htSessions = sessionTreatment==3

% doi4= layerTet(doiSessions,:)==4
% saline4=layerTet(salineSessions,:)==4

% figure
% for s=1:length(salineSessions)
% subplot(5,5,s)
% imagesc(squeeze(LFPallCh(salineSessions==s,LFPallCh,:,1,2)))
% end

%%% is averaging across treatments for each layer

 layerTet = layerSites(:,2:4:64)
% figure
% imagesc(squeeze(LFPallChDark(1,:,:,1,1)))

use =find(doiSessions==1)
for prepost=1:2
    figure
    if prepost==1, set(gcf,'Name','DOI LFP PRE');
    else set(gcf,'Name','DOI LFP POST'); end
    for s =1:sum(doiSessions)
        subplot(4,4,s)
        imagesc(squeeze(LFPallCh(use(s),:,:,1,prepost)))
    end
end

use =find(salineSessions==1)
for prepost=1:2
    figure
    if prepost==1, set(gcf,'Name','Saline LFP PRE');
    else set(gcf,'Name','Saline LFP POST'); end
    for s =1:sum(salineSessions)
        subplot(4,4,s)
        imagesc(squeeze(LFPallCh(use(s),:,:,1,prepost)))
    end
end

use = find(htSessions==1)
for prepost=1:2
    figure
    if prepost==1, set(gcf,'Name','5HT LFP PRE');
    else set(gcf,'Name','5HT LFP POST'); end
    for s =1:sum(htSessions)
        subplot(4,4,s)
        imagesc(squeeze(LFPallCh(use(s),:,:,1,prepost)))
    end
end

clear s use
use = find(salineSessions==1)
for mv = 1:2
figure
if mv==1, set(gcf,'Name', 'moving');
else set(gcf,'Name', 'stationary')
end
    for prepost =1:2
        for s =  1:length(use)
            subplot(4,5,s)
            plot(squeeze(LFPfreq(:,prepost)),squeeze(nanmean(LFPallCh(use(s), layerTet(s,:)==4,:,mv,prepost))));
            hold on
        end
    end
end

squeeze(nanmean(LFPallCh(sessionTreatment==1',layerTet(sessionTreatment==1,:)==4,:,1,1)))

%%evoked LFP all
figure
for t=1:4
subplot(2,2,t)
set(gcf, 'Name', 'evoked LFP')
plot(squeeze(LFPall(salineSessions==1,:,1,1)),'b');hold on;
plot(squeeze(LFPall(salineSessions==1,:,1,2)),'r'); xlabel 'Frequency (Hz)'; ylabel 'normalized power';
end
%
layerTet = layerSites(:,2:4:64)
use = layerTet(:,:)==4

figure
set(gcf, 'Name', 'all ch darkness LFP')
% for t=1:4
% subplot(2,2,t)
plot(squeeze(LFPallDark(:,:,1,1)),'b');hold on;
plot(squeeze(LFPallDark(:,:,1,2)),'r'); xlabel 'Frequency (Hz)'; ylabel 'normalized power';
% end

C = {[1 0 0],[.5 0 0]}; %red = mv=1, light =pre, dark =post
D = {[0 1 0],[0 .5 0]}; %green = mv=2
figure
for prepost =1:2
for tet=1:16
subplot(4,4,tet)
plot(squeeze(LFPallChDark(9,tet,:,1,prepost)),'Color',C{prepost});hold on; xlim([0 100]);
plot(squeeze(LFPallChDark(9,tet,:,2,prepost)),'Color',D{prepost}); xlim([0 100]);
xlabel 'Frequency (Hz)'; ylabel 'normalized power'; %mv
title(sprintf('layer %d',layerTet(tet)));
end
end

figure
subplot(1,2,1)
imagesc(squeeze(LFPallChDark(1,layerTet==3,:,1,1)));
subplot(1,2,2)
imagesc(squeeze(LFPallChDark(1,layerTet==3,:,1,2))); xlabel 'Frequency (Hz)'; ylabel 'normalized power';


%%% scatter plot of wn spont
for mv = 1:2
    figure
    for i = 1:6
        subplot(2,3,i)
        plot(wn_spont(data_wn &treatment==saline & layerAll ==i,mv,1),wn_spont(data_wn &treatment==saline& layerAll ==i,mv,2),'k.');
        % ylim ([-2 50]); xlim([-2 50]);
        hold on
        plot(wn_spont(data_wn &treatment==doi& layerAll ==i,mv,1),wn_spont(data_wn &treatment==doi& layerAll ==i,mv,2),'r.');
        plot(wn_spont(data_wn &treatment==ht& layerAll ==i,mv,1),wn_spont(data_wn &treatment==ht& layerAll ==i,mv,2),'w.');
        plot(wn_spont(data_wn &treatment==ketanserin& layerAll ==i,mv,1),wn_spont(data_wn &treatment==ketanserin& layerAll ==i,mv,2),'m.');
        plot(wn_spont(data_wn &treatment==ketandoi& layerAll ==i,mv,1),wn_spont(data_wn &treatment==ketandoi& layerAll ==i,mv,2),'c.');
        plot(wn_spont(data_wn &treatment==mglur2& layerAll ==i,mv,1),wn_spont(data_wn &treatment==mglur2& layerAll ==i,mv,2),'g.');
        plot(wn_spont(data_wn &treatment==lisuride& layerAll ==i,mv,1),wn_spont(data_wn &treatment==lisuride& layerAll ==i,mv,2),'b.');
        plot([0 10],[0 10]); axis equal
        title(sprintf('layer %d',i));  ylabel('post');
        if mv ==1 , xlabel('stop wn spont'); else  xlabel('move wn spont'); end
    end
end

figure
for l=1:5
    mean_fr(2,:) =mean(nanmean(drift_orient(data&treatment==doi & layerAll==l,:,2)))
    mean_fr(1,:) =mean(nanmean(drift_orient(data&treatment==doi & layerAll==l,:,1)))
    subplot(2,3,l)
    bar(mean_fr)
end


figure
for l=1:5
    osi(2,:) = squeeze(nanmean(drift_osi(data & treatment==doi &layerAll==l,1,2),1))
    %osierr(2,:) = squeeze(nanstd(drift_osi(data & treatment==doi & layerAll==l,1,2),1))/sqrt(sum(drift_osi(useN & layerAll==l,1,2)));
    osi(1,:) = squeeze(nanmean(drift_osi(data& treatment==doi  & layerAll==l,1,1),1))
    %osierr(1,:) = squeeze(nanstd(drift_osi(data &  treatment==doi &layerAll==l,1,1),1))/sqrt(sum(drift_osi(useN & layerAll==l,1,1)));
    subplot(2,3,l)
    % barweb(osi,osierr)
    bar(osi);hold on;
end

useOsi_pre = drift_osi(:,1,1)>0
useOsi_post = drift_osi(:,1,2)>0
useOsi=(useOsi_pre & useOsi_post)'
%useResp'
data = goodAll & useResp' &  useOsi==1 ;%which cells to include

figure
for t=1:3
    osi(1,:) = nanmean(abs(drift_osi(data&treatment==t,1,1)))
    osierr(1,:) = squeeze(nanstd(drift_osi(data & treatment==t,1,1)/sqrt(sum(drift_osi(data & treatment==t,1,1),1))));
    osi(2,:) = nanmean(abs(drift_osi(data&treatment==t,1,2)))
    osierr(2,:) = squeeze(nanstd(drift_osi(data & treatment==t,1,2)/sqrt(sum(drift_osi(data & treatment==t,1,2),1))));
    subplot(1,3,t)
    barweb(osi,osierr); ylim([0 1]); axis square
end

figure
for t=1:3
    osi(1,:) = nanmean(abs(drift_osi(data&treatment==t&layerAll==2|3,1,1)))
    osierr(1,:) = squeeze(nanstd(drift_osi(data & treatment==t&layerAll==2|3,1,1)/sqrt(sum(drift_osi(data & treatment==t&layerAll==2|3,1,1),1))));
    osi(2,:) = nanmean(abs(drift_osi(data&treatment==t&layerAll==5,1,2)))
    osierr(2,:) = squeeze(nanstd(drift_osi(data & treatment==t&layerAll==2|3,1,2)/sqrt(sum(drift_osi(data & treatment==t&layerAll==2|3,1,2),1))));
    subplot(1,3,t)
    barweb(osi,osierr); ylim([0 1]); axis square
    %[p h]=ranksum(drift_osi(data&treatment==t&layerAll==5,1,1),drift_osi(data&treatment==t&layerAll==5,1,2))
    % text(.9, 1, ['p = ' num2str(p)],'FontSize',18)
end



data = goodAll & useResp' & ~inhAll ;%which cells to include

%%% scatter plot wn evoked
for mv = 1:2
    figure
    for i = 1:6
        subplot(2,3,i)
        plot(wn_evoked(data &treatment==saline & layerAll ==i,mv,1),wn_evoked(data &treatment==saline& layerAll ==i,mv,2),'k.');
        hold on
        plot(wn_evoked(data &treatment==doi& layerAll ==i,mv,1),wn_evoked(data &treatment==doi& layerAll ==i,mv,2),'r.');
        plot(wn_evoked(data &treatment==ht& layerAll ==i,mv,1),wn_evoked(data &treatment==ht& layerAll ==i,mv,2),'w.');
        plot(wn_evoked(data &treatment==ketanserin& layerAll ==i,mv,1),wn_evoked(data &treatment==ketanserin& layerAll ==i,mv,2),'m.');
        plot(wn_evoked(data &treatment==ketandoi& layerAll ==i,mv,1),wn_evoked(data &treatment==ketandoi& layerAll ==i,mv,2),'c.');
        plot(wn_evoked(data &treatment==mglur2& layerAll ==i,mv,1),wn_evoked(data &treatment==mglur2& layerAll ==i,mv,2),'g.');
        plot(wn_evoked(data &treatment==lisuride& layerAll ==i,mv,1),wn_evoked(data &treatment==lisuride& layerAll ==i,mv,2),'b.');
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

titles={'Saline', 'DOI', '5HT'}
data=(goodAll==1)
figure
for t=1:3
    subplot(1,3,t)
    plot(meanRdark(~inhAll &goodAll & treatment==t,1), meanRdark(~inhAll& goodAll& treatment==t,2),'.','Markersize',12);hold on; axis square;
    set(gca,'FontSize',18)
    plot(meanRdark(inhAll &dataInh & treatment==t,1), meanRdark(inhAll& dataInh& treatment==t,2),'r.','Markersize',12);hold on; axis square;
    mdl_dark= fitlm(meanRdark(~inhAll &data & treatment==t,1), meanRdark(~inhAll& data& treatment==t,2),2)
    rsquared_dark(t) = mdl_dark.Rsquared.Ordinary
    xlim([0 5]);ylim([0 5]); %xlabel('Pre spikes/sec');ylabel('Post spikes/sec');
    % mpre=nanmean(mean(meanRdark(goodAll==1 & treatment==t,1),2))
    % mpost=nanmean(mean(meanRdark(goodAll==1 & treatment==t,2),2))
    % plot(mpre,mpost,'+k','Markersize',12,'Linewidth',2)
    plot([0 35],[0 35],'Linewidth',2);
    n_cells(t) = sum(goodAll==1 &treatment==t)
%     text(.25, 14, ['n = ' num2str(n_cells(t))],'FontSize',18)
%     text(.25, 13, ['r^2 = ' num2str(rsquared_dark(t))],'FontSize',18)
    %title(titles{t});
end


titles={'Saline', 'DOI', '5HT'}
figure
for t=1:3
    subplot(3,2,t)
    hist(meanRdark(treatment==t,1),1:2:40);xlim([0 35]);
    subplot(3,2,t+2)
    hist(meanRdark(treatment==t,2),1:2:40);xlim([0 35]);
    %     xlim([0 35]);ylim([0 35]); xlabel('Pre FR');ylabel('Post FR');
    %     plot([0 35],[0 35]);
    %plot(meanRdark(inhAll==1 &treatment==t,1), meanRdark(inhAll==1&treatment==t,2),'r.');hold on;
    title(titles{t});
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
            plot(wn_evoked(data &treatment==t & inhAll ==c,mv,1),wn_evoked(data &treatment==t & inhAll==c ,mv,2),'.');hold on; plot([0 45],[0 45]);axis equal; ylim([-5 40]);xlim([-5 40]);
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
        plot(meanRdark(data &treatment==t & inhAll==c,1),meanRdark(data &treatment==t & inhAll==c,2),'.');hold on; plot([0 45],[0 45]);axis equal;ylim([0 45]); xlim([0 45]);
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
        plot(drift_orient(useN(i),:,1,1),'Color',[0.5 0 0]);  plot(drift_orient(useN(i),:,2,1),'Color',[0 0.5 0]); %pre sal & DOI mv & stat dark =pre
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
for t = 1:3
    useDark = meanRdark(:,1)>1 | meanRdark(:,2)>1;
    useN = find(data & treatment==t)
    MIdark= (meanRdark(:,2)-meanRdark(:,1))./(meanRdark(:,2)+meanRdark(:,1));
    
    subplot(2,3,t+2)
    h= hist(MIdark(useDark(treatment==t)),-1:.1:1);
    Mbins=-1:.1:1
    bar(Mbins,h/sum(useDark(treatment==t)))
    xlim([-1.5 1.5]); ylim([0 .25]);axis xy
    xlabel('MI'); ylabel('fraction of cells');title(titles{t});
    set(gcf,'Name','MI Dark')
    %hold on
    %     subplot(2,3,2)
    %     meanMIdark(t) = nanmean(MIdark(treatment==t))
    %     %err(t) = nanstd(MIdark(treatment==t))/sqrt(sum(MIdark (treatment==t)));
    %     bar(meanMIdark)
    %     %barweb(meanMIdark(1:1:t)',err(1:1:t)');
    ylim([-0.2 0.2]);
end
%
% %% mod for corr in darkness %%

titles = {'Saline','DOI','5HT','Ketanserin', 'Ketanserin + DOI','MGluR2', 'MGluR2 + DOI', 'Lisuride'};
figure
for t = 1:8
    useN = find(data & treatment==t)
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
    useN = find(data & treatment==t)
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
