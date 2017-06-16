clear all
close all
dbstop if error

batchEphys_Ian; %%% load batch file
% set(groot,'defaultFigureVisible','off') %disable figure plotting
set(groot,'defaultFigureVisible','on')

%%% select the sessions you want based on filters
% use =  find(strcmp({files.notes},'good data'))%useSess = use;
% use =  find( strcmp({files.treatment},'5HT') & strcmp({files.notes},'good data'))

%for specific experiment:
use =  find(strcmp({files.notes},'good data') & strcmp({files.expt},'061417'))
sprintf('%d selected sessions',length(use))

saline=1; doi=2; ht=3; ketanserin=4; ketandoi=5; mglur2=6; mglur2doi=7; lisuride=8;
savePDF=0;
redo = 1;
n=0; ncorr=0; %%% number of units loaded, ncorr= number of correlation pairs


for i = 1:length(use)
    % close all
    %%% extract filenames
    afile = [pathname '\' files(use(i)).dir '\' files(use(i)).analysisfile '.mat'];
    clustfile = [pathname '\' files(use(i)).dir '\' files(use(i)).clusterfile '.mat'] ;
%     cfile = [{[pathname '\' files(use(i)).dir '\' files(use(i)).predark_camera '.mat']}; {[pathname '\' files(use(i)).dir '\' files(use(i)).postdark_camera '.mat']}]';
    
    %%% get cell type based on waveform
    [inh mid] = getWaveform(clustfile,afile,1);
    nc = length(inh); cellrange = n+1:n+nc;
    inhAll(cellrange) = inh;
    
    good = ones(1,nc);
    good(files(use(i)).badsites) =0;
    goodAll(cellrange)=good;
    

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
      %  try
     hasDrift(cellrange)=1;
        for prepost = 1:2
            drift = getDrift_mv(clustfile,afile,files(use(i)).blockDrift{prepost},0);
            drift_orient(cellrange,:,:,prepost)=drift.orient_tune;
            drift_sf(cellrange,:,:,prepost) = drift.sf_tune;
            drift_spont(cellrange,:,prepost) = drift.interSpont;
            drift_fl_spont(cellrange,:,prepost) = drift.spont;
            drift_osi(cellrange,:,prepost) = drift.cv_osi;
            drift_F1F0(:,cellrange,prepost)= drift.F1F0;
            drift_ot_tune(cellrange,:,:,prepost)=drift.orient_tune;
            drift_sf_trial{prepost} = drift.trialSF;
            drift_orient_trial{prepost} = drift.trialOrient;
            drift_pref_theta(cellrange,:,prepost) = drift.pref_theta;
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
%             
        end
   
         for prepost=1:2
            lfpMoveDrift = getLfpMovement(clustfile,afile,files(use(i)).blockDrift{prepost},0);
            LFPallDrift(i,:,:,prepost) =squeeze(nanmedian(lfpMoveDrift.meanSpect,1));
            if size(lfpMoveDrift.meanSpect,1)==16
                LFPallChDrift(i,:,:,:,prepost) = lfpMoveDrift.meanSpect;
                LFPfreqDrift(:,:) = lfpMoveDrift.freq;
                display('good lfp');
            else
                display('lfp wrong size')
            end
         end  
        %catch end
    else
        hasDrift(cellrange)=0;
        
        drift = NaN;
        drift_orient(cellrange,12,1:2,prepost)=NaN;
        drift_sf(cellrange,1:7,1:2,prepost) = NaN;
        drift_spont(cellrange,1:2,prepost) = NaN;
        drift_osi(cellrange,1:2,prepost) = NaN;
        drift_F1F0(1,cellrange,prepost)= NaN;
        drift_ot_tune(cellrange,12,1:2,prepost)=NaN;
        drift_sf_trial{prepost} = NaN;
        drift_orient_trial{prepost} = NaN;
        drift_pref_theta(cellrange,1:2,prepost) = NaN;
        for c = 1:length(cellrange)
            drift_trial_psth{cellrange(c),1,prepost} = NaN;
            mnPsth(cellrange(c),prepost,50) = NaN;
        end

   for cond = 1:72
                for mv = 1:2
%                     if mv ==1
%                         tr = find(drift.trialOrient(1:end-1) == ceil(cond/6) & drift.trialSF(1:end-1)==mod(cond-1,6)+1 & drift.frameSpd<1);
%                     else
%                         tr = find(drift.trialOrient(1:end-1) == ceil(cond/6) & drift.trialSF(1:end-1)==mod(cond-1,6)+1 & drift.frameSpd>1);
%                     end
                    drift_tcourse(cellrange,1:72,1:2,prepost,:) = NaN;
                    sf = mod(cond-1,6)+1; ori = ceil(cond/6);
                    drift_cond_tcourse(cellrange,mv,prepost,ori,sf,:) = NaN;
                end
            end
          

lfpMoveDrift = NaN;
LFPallDrift(i,200,1:2,prepost) =NaN;
LFPallChDrift(i,16,200,1:2,prepost) = NaN;
LFPfreqDrift(1,200) = NaN;

    end

    if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2})
        hasWn(cellrange)=1;
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
    else
        hasWn(cellrange)=0;
        %    wn = NaN;
        wn_crf(cellrange,1:20,1:2,prepost)= NaN;
        wn_spont(cellrange,1:2,prepost)= NaN;
        wn_evoked(cellrange,1:2,prepost)= NaN;
        wn_gain(1,cellrange,prepost)= NaN;
        wn_frameR(cellrange,:,1:2,prepost) = NaN;
        wn_reliability(1,cellrange,prepost) = NaN;
    end
    
%     if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2})
%         hasWn(cellrange)=1;
%         
%         for prepost=1:2
%             
%             sta = getSTA(clustfile,afile,files(use(i)).blockWn{prepost},0)
%                 sta_nx(cellrange,prepost) = sta.nx
%                 sta_ny(cellrange,prepost)=sta.ny
%                 sta_sigx(cellrange,prepost)=sta.sigx
%                 sta_sigy(cellrange,prepost)=sta.sigy
%                 sta_exp_var(cellrange,prepost)=sta.exp_var
%                 sta_all_fit(cellrange,prepost)=sta.all_fit;
%                 sta_all_img(cellrange,prepost)=sta.all_img;
%                  sta_params(cellrange,prepost)=sta.params;
% 
% 
%         end
%     else
%         hasWn(cellrange)=0;
        
% sta = NaN;
% sta_nx(cellrange,prepost) =  NaN;
% sta_ny(cellrange,prepost)= NaN;
% sta_sigx(cellrange,prepost)= NaN;
% sta_sigy(cellrange,prepost)= NaN;
% sta_exp_var(cellrange,prepost)= NaN;
% % % sta_all_fit(cellrange,prepost) = NaN;
% % % sta_all_img(cellrange,prepost)= NaN;
% %sta_params(cellrange,prepost) = [];
% 
% 
% for c = 1:length(cellrange)
%     sta_all_fit{cellrange(c),prepost} = NaN;
%     sta_all_img{cellrange(c),prepost} = NaN;
% %     sta_params(cellrange(c),prepost)= [];
% end
%    end

     if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2})

        %% lfp power
        %     %%%(right now averages over all sites, should use layer info)
     try
        for prepost=1:2
            lfpMove = getLfpMovement(clustfile,afile,files(use(i)).blockWn{prepost},0);
            LFPall(i,:,:,prepost) =squeeze(nanmedian(lfpMove.meanSpect, 0));
            if size(lfpMove.meanSpect,1)==16
                LFPallCh(i,:,:,:,prepost) = lfpMove.meanSpect;
                LFPfreq(:,:) = lfpMove.freq;
                display('good lfp');
                lfpMove
            else
                display('lfp wrong size')
         end
         end  
        catch end
     else
                LFPallCh(i,200,1:2,:,prepost) = NaN;
                LFPfreq(1,200) = NaN;
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
                LFPallDark(i,200,1:2,prepost) = NaN;

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
        cv2Wn(cellrange,:) = cv2;
    else wnCorr(corrRange,1) = NaN; 
         wnCorr(corrRange,2)=NaN;
         cv2Wn(cellrange,:) = NaN;
    end
    if ~isempty(files(use(i)).blockDrift{1}) & ~isempty(files(use(i)).blockDrift{2})
        dt = 1;
        [preCorr postCorr cv2 R eigs] =  prepostDOIdarkness(clustfile,afile,files(use(i)).blockDrift,dt,0);
        driftCorr(corrRange,1) = preCorr(:); driftCorr(corrRange,2)=postCorr(:);
        cv2Drift(cellrange,:) = cv2;
    else driftCorr(corrRange,1) = NaN; 
         driftCorr(corrRange,2)=NaN;
         cv2Drift(cellrange,:) = NaN;
    end
    
    %%%% prepost correlation in darkness
    if ~isempty(files(use(i)).blockDark{1}) & ~isempty(files(use(i)).blockDark{2})
        
        dt = 1;
        [preCorr postCorr cv2 R eigs] = prepostDOIdarkness(clustfile,afile,files(use(i)).blockDark,dt,0);
        darkCorr(corrRange,1) = preCorr(:); darkCorr(corrRange,2)=postCorr(:);
        cv2Dark(cellrange,:) = cv2;
        meanRdark(cellrange,:) = mean(R,2);
    else
        cv2Dark(cellrange,:) = NaN;
        meanRdark(cellrange,:) = NaN;
        darkCorr(corrRange,1) = NaN;
        darkCorr(corrRange,2)= NaN;
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


% %%% plot speed histogram
% titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
% figure 
% for t =1:8
% subplot(2,4,t)
% plot(0.5:1:15,squeeze(mean(speedHistWn(sessionTreatment==t,1:15,:),1))); 
% title(titles{t}); xlabel('speed');axis square; 
% end
% legend('pre','post')

figure 
for t =1:3
subplot(1,3,t)
plot(0.5:1:25,squeeze(mean(speedHistDrift(sessionTreatment==t,1:25,:),1))); 
title(titles{t}); xlabel('speed');axis square; 
end
legend('pre','post')


%%%filter data before plotting%%%
clear max_wn
% low_wn = squeeze(min(wn_crf,[],2));
max_wn = wn_evoked-wn_spont% squeeze(max(wn_evoked,[],2))-(squeeze(wn_spont))
amp_wn = max_wn  %max_wn-low_wn
useResp = amp_wn(:,1,1)>3 | amp_wn(:,1,2)>3 %| amp_wn(:,2,1)>3 | amp_wn(:,2,2)>3;
data_wn = goodAll==1 & useResp';

clear max amp low
peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont
%peakresp = squeeze(mean(drift_orient,2))-squeeze(drift_spont);%subtract driftspont
% low = squeeze(min(drift_orient,[],2));
amp = peakresp; 
spontAmp = squeeze(max(drift_spont,[],2));
useResp = amp(:,1,1)>2 | amp(:,1,2)>2 %|amp(:,2,1)>2 | amp(:,2,2)>2;
data = goodAll & useResp' & ~inhAll & hasDrift==1 ;%which cells to include
dataInh =goodAll & useResp' & inhAll & hasDrift==1 ;
dataAll =goodAll & useResp' & hasDrift==1 ;

%evoked FR from drift
titles ={'Saline', 'DOI','5HT'};
for mv =1:2
    figure
    if mv==1
        set(gcf,'Name','mean stat drift FR')
    else
        set(gcf,'Name','mean mv drift FR')
    end
    for t=1:3
        subplot(1,3,t)
        plot(peakresp(~inhAll & goodAll==1 & treatment==t & hasDrift==1,mv,1),peakresp(~inhAll&goodAll==1&treatment==t& hasDrift==1,mv,2),'.','Markersize',10);
        hold on;
        set(gca,'FontSize',18);
        axis square;
        % xlim([0 18]);ylim([0 18]); %title(titles{t}); xlabel('Pre spikes/sec');ylabel('Post spikes/sec');
        ylim ([-10 35]);xlim([-10 35]);
        mpre=nanmean(peakresp(goodAll==1 & treatment==t& hasDrift==1,mv,1))
        mpost=nanmean(peakresp(goodAll==1 & treatment==t& hasDrift==1,mv,2))
        %         allPre = drift_orient(goodAll==1 & treatment==t& hasDrift==1,mv,1,1)
        %         allPost = drift_orient(goodAll==1 & treatment==t& hasDrift==1,mv,1,2)
        allPre = peakresp(goodAll==1 & treatment==t& hasDrift==1,mv,1)
        allPost = peakresp(goodAll==1 & treatment==t& hasDrift==1,mv,2) 
%         c = corrcoef(allPre,allPost)
%         c(t) = c(2,1);
%         mdl_drift = fitlm(allPre,allPost)
%         rsquared_drift(t) = mdl_drift.Rsquared.Ordinary
        plot(mpre,mpost,'+k','Markersize',18,'Linewidth',4)
        plot([0 35], [0 35])
        plot(peakresp(goodAll==1 & inhAll==1  &treatment==t& hasDrift==1,mv,1),peakresp(goodAll==1 & inhAll==1 & hasDrift==1&treatment==t,mv,2),'r.','Markersize',10);
        n_cells(t) = sum(goodAll==1 &treatment==t& hasDrift==1)
        text(2, 32, ['n = ' num2str(n_cells(t))],'FontSize',20)
%         text(2, 30, ['r^2 = ' num2str(rsquared_drift(t))],'FontSize',20)
%        text(2, 30, ['c.c. = ' num2str(c(t))],'FontSize',20)

    end
end


%Drift-evoked Proportion of cells suppressed or facilitated pre/post
figure
for mv = 1:2
    for t=1:3
        preR = mean(drift_orient(goodAll==1&treatment==t& hasDrift==1,:,mv,1),2);
        postR = mean(drift_orient(goodAll==1  &treatment==t& hasDrift==1,:,mv,2),2);
        deltaR = postR-preR;
        stderror = nanstd(deltaR) / sqrt(length(deltaR));
        fracSupp(t) = sum(deltaR<-2)/length(deltaR);
        stderrorSupp(t) = nanstd(fracSupp) / sqrt(length(fracSupp));
        fracInc(t) = sum(deltaR>2)/length(deltaR);
        stderrorInc(t) = nanstd(fracInc) / sqrt(length(fracInc));
    end
    subplot(1,2,mv)
    set(gcf,'Name','mean drift fraction, evoked')
    if mv==1, title ('stationary'), else title('moving'); end
    %else set(gcf,'Name','mean mv drift fractions, evoked')
    %end
    hold on; errorbar(fracSupp,stderrorSupp, '-o','LineWidth',2);
    errorbar(fracInc,stderrorInc,'-^','MarkerEdgeColor','r','Color','r','LineWidth', 2);
    legend('Suppressed','Facilitated');
    axis square; ylabel('Proportion of cells','FontSize',18); ylim([0 0.2]);
    Labels = {'Saline','DOI', '5HT'};set(gca, 'XTick', 1:3, 'XTickLabel', Labels,'FontSize',18);
    if mv==1, title ('stationary'), else title('moving'); end
end

figure
for mv = 1:2
    for t=1:3
        preR = mean(drift_spont(goodAll==1&treatment==t& hasDrift==1,mv,1),2);
        postR = mean(drift_spont(goodAll==1  &treatment==t& hasDrift==1,mv,2),2);
        deltaR = postR-preR;
        fracSupp(t) = sum(deltaR<-2)/length(deltaR);
        fracInc(t) = sum(deltaR>2)/length(deltaR);
        stderrorSupp(t) = nanstd(fracSupp) / sqrt(length(fracSupp));
        fracInc(t) = sum(deltaR>2)/length(deltaR);
        stderrorInc(t) = nanstd(fracInc) / sqrt(length(fracInc));
    end
    
    subplot(1,2,mv)
    set(gcf,'Name','mean drift fractions, spont')
%     if mv==1
%         set(gcf,'Name','mean stat drift fractions, spont')
%     else set(gcf,'Name','mean mv drift fractions, spont')
%     end
    errorbar(fracSupp, stderrorSupp,'-o','LineWidth', 2); hold on;
    errorbar(fracInc,stderrorInc,'-^','MarkerEdgeColor','r','Color','r','LineWidth', 2); % plot(1-fracInc-fracSupp);
    legend('Suppressed','Facilitated');
    axis square; ylabel('Proportion of cells','FontSize',18); ylim([0 0.15]);
    Labels = {'Saline','DOI', '5HT'};set(gca, 'XTick', 1:3, 'XTickLabel', Labels,'FontSize',18);
   if mv==1, title ('stationary'), else title('moving'); end

end

for t=1:3
        preRdark= meanRdark(goodAll==1 &treatment==t,1);
        postRdark = meanRdark(goodAll==1&treatment==t,2);
        deltaRdark = postRdark-preRdark;
        fracSuppDark(t) = sum(deltaRdark<-2)/length(deltaRdark);
        fracIncDark(t) = sum(deltaRdark>2)/length(deltaRdark);
        stderrorSuppDark(t) = nanstd(fracSuppDark) / sqrt(length(fracSuppDark));
        stderrorIncDark(t) = nanstd(fracIncDark) / sqrt(length(fracIncDark));
    end
    figure
    errorbar(fracSuppDark, stderrorSuppDark,'-o','LineWidth', 2); hold on; 
    errorbar(fracIncDark,stderrorIncDark,'-^','MarkerEdgeColor','r','Color','r','LineWidth', 2); % plot(1-fracInc-fracSupp);
    legend('Suppressed','Facilitated'); title('Spontanous(dark)');
    axis square; ylabel('Proportion of cells','FontSize',18); ylim([0 0.15]);
    Labels = {'Saline','DOI', '5HT'};set(gca, 'XTick', 1:3, 'XTickLabel', Labels,'FontSize',18);



%%MI for evoked drift, stationary %%
layerz={'l1','l2','l3','layer 4','layer 5','L6'};
for t=1:3
    figure
    if t == 1
        set(gcf,'Name', 'hist saline');
    elseif t==2
        set(gcf,'Name', 'mi hist doi');
    else
        set(gcf,'Name', 'mi hist 5ht');
    end
    useN=goodAll==1 & treatment==t & hasDrift==1
    miDrift= (mean(drift_orient(useN &(layerAll==2|layerAll==3),:,1,2),2)-mean(drift_orient(useN&(layerAll==2|layerAll==3),:,1,1),2))./...
        (mean(drift_orient(useN&(layerAll==2|layerAll==3),:,1,2),2)+mean(drift_orient(useN&(layerAll==2|layerAll==3),:,1,1),2));
    h= hist(miDrift,-1:.2:1);
    Mbins=-1:.2:1
    subplot(1,3,1)
    bar(Mbins,h/sum(useN &(layerAll==2|layerAll==3)),'FaceColor',[0 .5 .5],'Linewidth',2); xlim([-1.5 1.5]);axis square; title('layer 2/3')%ylim([0 .4])
    for i=4:5
        miDrift= (mean(drift_orient(useN &(layerAll==i),:,1,2),2)-mean(drift_orient(useN&(layerAll==i),:,1,1),2))./(mean(drift_orient(useN&(layerAll==i),:,1,2),2)...
            +mean(drift_orient(useN&(layerAll==i),:,1,1),2));
        h= hist(miDrift,-1:.2:1);
        subplot(1,3,i-2)
        Mbins=-1:.2:1
        bar(Mbins,h/sum(useN&layerAll==i),'FaceColor',[0 .5 .5],'Linewidth',2);xlim([-1.5 1.5]); axis square %ylim([0 .5])
        title(layerz{i});
    end
end

%%MI for evoked drift, move%%
layerz={'l1','l2','l3','layer 4','layer 5'};
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
    miDrift= (mean(drift_orient(useN &(layerAll==2|layerAll==3),:,2,2),2)-mean(drift_orient(useN&(layerAll==2|layerAll==3),:,2,1),2))./(mean(drift_orient(useN&(layerAll==2|layerAll==3),:,2,2),2)+mean(drift_orient(useN&(layerAll==2|layerAll==3),:,2,1),2));
    h= hist(miDrift,-1:.2:1);
    Mbins=-1:.2:1
    subplot(1,3,1)
    bar(Mbins,h/sum(useN &(layerAll==2|layerAll==3)),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .4]); xlim([-1.5 1.5]);axis square; title('layer 2/3')
    for i=4:5
        miDrift= (mean(drift_orient(useN &(layerAll==i),:,2,2),2)-mean(drift_orient(useN&(layerAll==i),:,2,1),2))./(mean(drift_orient(useN&(layerAll==i),:,2,2),2)+mean(drift_orient(useN&(layerAll==i),:,2,1),2));
        h= hist(miDrift,-1:.2:1);
        subplot(1,3,i-2)
        Mbins=-1:.2:1
        bar(Mbins,h/sum(useN&layerAll==i),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .4]);xlim([-1.5 1.5]); axis square
        title(layerz{i});
        
    end
end

%spont from drift
titles ={'Saline', 'DOI','5HT'};
for mv =1:2
    figure
    if mv==1
        set(gcf,'Name','mean stationary')
    else
        set(gcf,'Name','mv')
    end
    for t=1:3
        subplot(1,3,t)
        plot(drift_spont(goodAll==1 &treatment==t &inhAll==0 & hasDrift==1,mv,1),drift_spont(goodAll==1&treatment==t&inhAll==0& hasDrift==1,mv,2),'.','Markersize',10);
        ylim([0 35]); xlim([0 35]);
        xlabel('Pre spikes/sec');ylabel('Post spikes/sec');
        axis square
        hold on;
        plot([0 35], [0 35])
        plot(drift_spont(goodAll==1 & inhAll==1 &treatment==t& hasDrift==1,mv,1),...
            drift_spont(goodAll==1& inhAll==1&treatment==t& hasDrift==1,mv,2),'r.','Markersize',10);
%        mdl_dSpont= fitlm(drift_spont(goodAll&treatment==t,mv,1),drift_spont(goodAll&treatment==t,mv,2),2)
%        rsquared_driftSpont(t) = mdl_dSpont.Rsquared.Ordinary
        title(titles{t})
        n_cells(t) = sum(goodAll==1 &treatment==t& hasDrift==1)
%         text(2, 32, ['n = ' num2str(n_cells(t))],'FontSize',20)
       % text(2, 30, ['n = ' num2str(rsquared_driftSpont(t))],'FontSize',20)

    end
end

%MI for drift spont mv
layerz={'l1','l2','l3','layer 4','layer 5'};
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
    miDrift= (mean(drift_spont(useN &(layerAll==2|layerAll==3),2,2),2)-mean(drift_spont(useN&(layerAll==2|layerAll==3),2,1),2))./(mean(drift_spont(useN&(layerAll==2|layerAll==3),2,2),2)+mean(drift_spont(useN&(layerAll==2|layerAll==3),2,1),2));
    h= hist(miDrift,-1:.2:1);
    Mbins=-1:.2:1
    subplot(1,3,1)
    bar(Mbins,h/sum(useN &(layerAll==2|layerAll==3)),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .4]); xlim([-1.5 1.5]);axis square; title('layer 2/3')
    for i=4:5
        miDrift= (mean(drift_spont(useN &(layerAll==i),2,2),2)-mean(drift_spont(useN&(layerAll==i),2,1),2))./(mean(drift_spont(useN&(layerAll==i),2,2),2)+mean(drift_spont(useN&(layerAll==i),2,1),2));
        h= hist(miDrift,-1:.2:1);
        subplot(1,3,i-2)
        Mbins=-1:.2:1
        bar(Mbins,h/sum(useN&layerAll==i),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .4]);xlim([-1.5 1.5]); axis square
        title(layerz{i});
    end
end

% useDrift = (abs((amp(:,1,1)-amp(:,1,2)))>2 | abs((amp(:,2,1)-amp(:,2,2))>2))';
% useSpont = (abs(spontAmp(:,1)-spontAmp(:,2))>2)';
% % useN=goodAll==1 & hasDrift==1 &useDrift==1; 
% useC = goodAll==1 & hasDrift==1 & useDrift==0;
% useSpont = goodAll==1 & hasDrift==1 & useSpont ==1&(useResp==1)';;
% useSpontC = goodAll==1 & hasDrift==1 & useSpont ==0&(useResp==1)';;
% % allC = goodAll==1 & hasDrift==1;
% nullmiDrift= (mean(drift_orient(useC,:,1,2),2)-mean(drift_orient(useC,:,1,1),2))./(mean(drift_orient(useC,:,1,2),2)+mean(drift_orient(useC,:,1,1),2));
% %miDrift= (mean(drift_orient(useN,:,1,2),2)-mean(drift_orient(useN,:,1,1),2))./(mean(drift_orient(useN,:,1,2),2)+mean(drift_orient(useN,:,1,1),2));
% figure
% subplot(2,2,1)
% for t=1:3
%     useN=goodAll==1 & hasDrift==1 &treatment==t&useDrift==1&(useResp==1)';;
%     allC = goodAll==1 & hasDrift==1 & treatment==t&(useResp==1)';
%     miDrift= (mean(drift_orient(useN,:,1,2),2)-mean(drift_orient(useN,:,1,1),2))./(mean(drift_orient(useN,:,1,2),2)+mean(drift_orient(useN,:,1,1),2));
%     suppressedEv = miDrift<0
%     fracSuppressedEv(t) = sum(suppressedEv)/sum(allC)
%     bar(fracSuppressedEv);ylabel('proportion of cells suppressed');
%     title ('suppression of evoked rate');
%     Labels = {'Saline', 'DOI', '5HT'};
%     set(gca, 'XTick', 1:3, 'XTickLabel', Labels,'FontSize',12);
% end
% subplot(2,2,2)
% for t=1:3
%     useN=goodAll==1 & hasDrift==1 &treatment==t &useSpont==1&(useResp==1)';
%     allC = goodAll==1 & hasDrift==1 & treatment==t&(useResp==1)';
%     miDrift= (mean(drift_spont(useN,1,2),2)-mean(drift_spont(useN,1,1),2))./(mean(drift_spont(useN,1,2),2)+mean(drift_spont(useN,1,1),2));
%     suppressedSpont = miDrift<0
%     fracSuppressedSpont(t) = sum(suppressedSpont)/sum(allC)
%     bar(fracSuppressedSpont);ylabel('proportion of cells suppressed');
%     title ('suppression of spontaneous rate');
%     Labels = {'Saline', 'DOI', '5HT'};
%     set(gca, 'XTick', 1:3, 'XTickLabel', Labels,'FontSize',12);
% end
% subplot(2,2,3)
% for t=1:3
%     useN=goodAll==1 & hasDrift==1 &treatment==t &useDrift==1&(useResp==1)';
%     allC = goodAll==1 & hasDrift==1 & treatment==t&(useResp==1)';
%     miDrift= (mean(drift_orient(useN,:,1,2),2)-mean(drift_orient(useN,:,1,1),2))./(mean(drift_orient(useN,:,1,2),2)+mean(drift_orient(useN,:,1,1),2));
%     incEv = miDrift>0
%     fracIncEv(t) = sum(incEv)/sum(allC)
%     bar(fracIncEv);ylabel('proportion of cells facilitated');
%     title ('facilitation of evoked rate');
%     Labels = {'Saline', 'DOI', '5HT'};
%     set(gca, 'XTick', 1:3, 'XTickLabel', Labels,'FontSize',12);
% end
% subplot(2,2,4)
% clear fracInc
% for t=1:3
%     useN=goodAll==1 & hasDrift==1 &treatment==t &useSpont==1 & (useResp==1)';
%     allC = goodAll==1 & hasDrift==1 & treatment==t&(useResp==1)';
%     miDrift= (mean(drift_spont(useN,1,2),2)-mean(drift_spont(useN,1,1),2))./(mean(drift_spont(useN,1,2),2)+mean(drift_spont(useN,1,1),2));
%     incSpont = miDrift>0
%     fracIncSpont(t) = sum(incSpont)/(sum(allC))
%     bar(fracIncSpont);ylabel('proportion of cells facilitated');
%     title ('facilitation of spontaneous rate');
%     Labels = {'Saline', 'DOI', '5HT'};
%     set(gca, 'XTick', 1:3, 'XTickLabel', Labels,'FontSize',12);
% end
% 
% for t=1:3
% useC = goodAll==1 & hasDrift==1 & treatment==t & useDrift==0&(useResp==1)';;
% useSpontC = goodAll==1 & hasDrift==1 & useSpont ==0 & treatment==t&(useResp==1)';;
% allC = goodAll==1 & hasDrift==1 & treatment==t&(useResp==1)';;
% % nullmiDriftEv= (mean(drift_orient(useC,:,1,2),2)-mean(drift_orient(useC,:,1,1),2))./(mean(drift_orient(useC,:,1,2),2)+mean(drift_orient(useC,:,1,1),2));
% % nullmiDriftSpont= (mean(drift_spont(useSpontC,1,2),2)-mean(drift_spont(useSpontC,1,1),2))./(mean(drift_spont(useSpontC,1,2),2)+mean(drift_spont(useSpontC,1,1),2));
% nochangeEv(t) = (sum(useC)/sum(allC))
% nochangeSpont(t) = (sum(useSpontC)/sum(allC))
% end
% figure
% subplot(1,2,1)
% plot(fracSuppressedEv,'LineWidth', 2);hold on; plot(fracIncEv,'r','LineWidth', 2); plot(nochangeEv,'g','LineWidth', 2);
% axis square; ylabel('proportion of cells');ylim([0 .6]);
% Labels = {'Saline', 'DOI', '5HT'};set(gca, 'XTick', 1:3, 'XTickLabel', Labels,'FontSize',14);
% legend('Suppressed', 'Facilitated')%, 'No change');
% title ('Evoked rate');
% 
% subplot(1,2,2)
% plot(fracSuppressedSpont,'LineWidth', 2); ylim([0 .2]);hold on;
% plot(fracIncSpont,'r','LineWidth', 2);%plot(nochangeSpont,'g','LineWidth', 2); 
% title ('Spontaneous rate'); axis square;ylabel('proportion of cells')
% set(gca, 'XTick', 1:3, 'XTickLabel', Labels,'FontSize',12);


clear cycR
for f = 1:20;
    cycR(:,f,:,:) = nanmean(wn_frameR(:,f:20:end,:,:),2);
end
spont = squeeze(mean(cycR(:,[1 2 19 20],:,:),2));
evoked = squeeze(mean(cycR(:,9:11,:,:),2)) - spont;

clear max_wn
max_wn = wn_evoked-wn_spont% squeeze(max(wn_evoked,[],2))-(squeeze(wn_spont))
amp_wn = max_wn  %max_wn-low_wn
useResp = amp_wn(:,1,1)>3 | amp_wn(:,1,2)>3 %| amp_wn(:,2,1)>3 | amp_wn(:,2,2)>3;
data_wn = goodAll==1 & useResp';
%wn evoked mv and stationary
titles = {'Saline','DOI','5HT'}
for mv=1:2
    figure
if mv==1, set(gcf,'Name', 'wn evoked stationary');
else set(gcf,'Name', 'wn evoked mv');end

for t=1:3
    subplot(1,3,t)
    plot(amp_wn(goodAll==1&treatment==t&~inhAll & hasWn==1,mv,1),...
        amp_wn(goodAll==1&treatment==t&~inhAll& hasWn==1,mv,2),'.', 'Markersize',10);hold on;axis square;
    xlim([-10 30]); ylim([-10 30]);
    plot([-30 30],[-30 30]);
    plot(amp_wn(goodAll==1&treatment==t&inhAll==1& hasWn==1,mv,1),...
        amp_wn(goodAll==1&treatment==t&inhAll==1& hasWn==1,mv,2),'.r','Markersize',10);
    pre= amp_wn(goodAll==1&treatment==t& hasWn==1,mv,1)
    post=amp_wn(goodAll==1&treatment==t& hasWn==1,mv,2)
    title(titles{t})
    c = corrcoef(pre,post)
    c = c(2,1);
%     mdl_wnEv= fitlm(amp_wn(goodAll&treatment==t,mv,1),amp_wn(goodAll&treatment==t,mv,2),2)
%     rsquared_wnEv(t) = mdl_wnEv.Rsquared.Ordinary
    n_cells(t) = sum(goodAll==1 &treatment==t& hasWn==1)
    text(-7, 28, ['n = ' num2str(n_cells(t))],'FontSize',20)
%     text(-7, 26, ['r^2 = ' num2str(rsquared_wnEv(t))],'FontSize',20)
    text(-7, 26, ['c.c = ' num2str(c)],'FontSize',20)

end
end

titles = {'Saline','DOI','5HT'}
for mv=1:2
   figure
if mv==1, set(gcf,'Name', 'wn spont stationary');
else set(gcf,'Name', 'wn spont mv');end
for t=1:3
    subplot(1,3,t)
    plot(wn_spont(goodAll==1&treatment==t&~inhAll& hasWn==1,mv,1),wn_spont(goodAll==1&treatment==t&~inhAll& hasWn==1,mv,2),'.', 'Markersize', 10);
    hold on;axis square;
    %     xlim([min(wn_spont(data_wn,mv,1))-.5 max(wn_spont(data_wn,mv,1))+.5]); ylim([min(wn_spont(data_wn,mv,1))-.5 max(wn_spont(data_wn,mv,1))+.5]);
    plot([-30 30],[-30 30]);
    plot(wn_spont(goodAll==1&treatment==t&inhAll==1& hasWn==1,mv,1),wn_spont(goodAll==1&treatment==t&inhAll==1& hasWn==1,mv,2),'.r', 'Markersize', 10);
    xlim([0 25]); ylim([0 25]);
    title(titles{t})
%     mdl_wnSpont= fitlm(wn_spont(goodAll&treatment==t,mv,1),wn_spont(goodAll&treatment==t,mv,2),2)
%     rsquared_wnSpont(t) = mdl_wnSpont.Rsquared.Ordinary
    pre= wn_spont(goodAll==1&treatment==t& hasWn==1,mv,1)
    post=wn_spont(goodAll==1&treatment==t& hasWn==1,mv,2)
    c = corrcoef(pre,post)
    c = c(2,1);
    n_cells(t) = sum(goodAll==1 &treatment==t& hasWn==1)
    text(2, 23, ['n = ' num2str(n_cells(t))],'FontSize',20)
    text(2, 21, ['c.c = ' num2str(c)],'FontSize',20)

%     text(2, 21, ['r^2 = ' num2str(rsquared_wnSpont(t))],'FontSize',20)
end
end

%evoked calculcated from cyc avg
clear n_cells
for mv = 1:2
    figure
    if mv==1, set(gcf,'Name', 'evoked pre vs post stationary')
    else set(gcf,'Name', 'evoked pre vs post mv'), end
    for t=1:3
        subplot(1,3,t)
        plot(evoked(find(goodAll& treatment==t),mv,1),evoked(find(goodAll& treatment==t),mv,2),'.','MarkerSize',10);
        title(titles{t});
        hold on; plot([-50 50],[-50 50]); axis xy
        %xlabel('pre'); ylabel('post');% xlim([min(evoked(data_wn,mv,1))-.5 max(evoked(data_wn,mv,1))+.5]); ylim([min(evoked(data_wn,mv,1))-.5 max(evoked(data_wn,mv,1))+.5]);
        plot(evoked(find(goodAll==1&inhAll& treatment==t),mv,1),evoked(find(goodAll==1&inhAll& treatment==t),mv,2),'.r','MarkerSize',10); axis square;
        xlim([-10 30]);ylim([-10 30])
        %         mdl_ev= fitlm(evoked(goodAll&treatment==t,mv,1),evoked(goodAll&treatment==t,mv,2),2)
        %         rs_ev(t) = mdl_ev.Rsquared.Ordinary
        %         text(-8, 28, ['r^2 = ' num2str(rs_ev(t),'%.2f')],'FontSize',20)
        pre= evoked(goodAll==1&treatment==t& hasWn==1,mv,1)
        post=evoked(goodAll==1&treatment==t& hasWn==1,mv,2)
        c = corrcoef(pre,post)
        c = c(2,1);
        text(-8, 23, ['c.c = ' num2str(c)],'FontSize',20)
        n_cells(t) = sum(goodAll==1 &treatment==t& hasWn==1)
        text(-8, 25, ['n = ' num2str(n_cells(t))],'FontSize',20)
        set(gca,'FontSize',20)
        text(-8, 23, ['c.c = ' num2str(c)],'FontSize',20)

    end
end

%spont calculcated from cyc avg
titles = {'Saline','DOI','5HT'};
for mv = 1:2
    figure
    if mv==1, set(gcf,'Name', 'wn spont stationary')
    else set(gcf,'Name', 'wn spont mv'), end
    for t=1:3
        subplot(1,3,t)
        plot(spont(goodAll==1 & treatment==t& hasWn==1 & ~inhAll,mv,1),spont(goodAll==1& treatment==t& hasWn==1& ~inhAll,mv,2),'.', 'Markersize',10); 
        title(titles{t}); hold on; axis square; plot([0 50],[0 50]);
        xlabel('pre'); ylabel('post');
        plot(spont(goodAll==1& inhAll & treatment==t& hasWn==1,mv,1),spont(goodAll==1& inhAll& treatment==t& hasWn==1,mv,2),'.r', 'Markersize',10);
        hold on; plot([0 50],[0 50]);
        axis([0 35 0 35]);
        %         mdl_spont= fitlm(spont(goodAll&treatment==t& hasWn==1,mv,1),spont(goodAll&treatment==t& hasWn==1,mv,2),2)
        %         rs_spont(t) = mdl_spont.Rsquared.Ordinary
        %         text(2, 33, ['r^2 = ' num2str(rs_spont(t),'%.2f')],'FontSize',20)
        pre= spont(goodAll==1&treatment==t& hasWn==1,mv,1)
        post=spont(goodAll==1&treatment==t& hasWn==1,mv,2)
        c = corrcoef(pre,post)
        c = c(2,1);
        text(2, 28, ['c.c = ' num2str(c)],'FontSize',20)
        n_cells(t) = sum(goodAll==1 &treatment==t& hasWn==1)
        text(2, 30, ['n = ' num2str(n_cells(t))],'FontSize',20)
    end
end


clear data
titles={'Saline', 'DOI', '5HT'}
data=(meanRdark)
figure
set(gcf,'Name','mean dark FR')
for t=1:3
    subplot(1,3,t)
    plot(meanRdark(~inhAll & treatment==t,1), meanRdark(~inhAll &treatment==t,2),'.','Markersize',12);hold on; axis square;
    set(gca,'FontSize',18)
    plot(meanRdark(inhAll  & treatment==t,1), meanRdark(inhAll &treatment==t,2),'r.','Markersize',12);hold on; axis square;
%     mdl_dark= fitlm(meanRdark(treatment==t,1), meanRdark(treatment==t,2),2)
%     rsquared_dark(t) = mdl_dark.Rsquared.Ordinary
    %xlim([0 20]);ylim([0 20]); %xlabel('Pre spikes/sec');ylabel('Post spikes/sec');
    mpre=nanmean(mean(meanRdark(goodAll==1 & treatment==t,1),2))
    mpost=nanmean(mean(meanRdark(goodAll==1 & treatment==t,2),2))
    plot(mpre,mpost,'+k','Markersize',12,'Linewidth',2)
    plot([0 30],[0 30],'Linewidth',2);
    pre= meanRdark(goodAll==1&treatment==t& hasWn==1,1)
    post=meanRdark(goodAll==1&treatment==t& hasWn==1,2)
%     c = corrcoef(pre,post)
%     c = c(2,1);
%     text(2, 28, ['c.c = ' num2str(c)],'FontSize',20)
    n_cells(t) = sum(goodAll==1 &treatment==t)
%     text(2, 28, ['n = ' num2str(n_cells(t))],'FontSize',18)
%     text(2, 25, ['r^2 = ' num2str(rsquared_dark(t))],'FontSize',18)
    title(titles{t});
end


%spontaneous (dark) Proportion of cells suppressed or facilitated pre/post
for t=1:3
        preRdark= meanRdark(goodAll==1 &treatment==t,1)
        postRdark = meanRdark(goodAll==1&treatment==t,2)
        deltaRdark = postRdark-preRdark
        fracSuppDark(t) = sum(deltaRdark<-2)/length(deltaRdark)
        fracIncDark(t) = sum(deltaRdark>2)/length(deltaRdark)
    end
figure
    plot(fracSuppDark,'-o','LineWidth', 2); hold on; 
    plot(fracIncDark,'-^','MarkerEdgeColor','r','Color','r','LineWidth', 2); % plot(1-fracInc-fracSupp);
    legend('Suppressed','Facilitated'); title('Spontanous(dark)')
    axis square; ylabel('Proportion of cells','FontSize',18); ylim([0 0.3]);
    Labels = {'Saline','DOI', '5HT'};set(gca, 'XTick', 1:3, 'XTickLabel', Labels,'FontSize',18);

% titles ={'Saline', 'DOI','5HT'};
% figure
% set(gcf,'Name','cv2 wn');
% for t=1:3
%     subplot(1,3,t)
%     mdl_cv2wn= fitlm(cv2Wn(goodAll==1&treatment==t,1), cv2Wn(treatment==t & goodAll==1,2),2)
%     rs_cv2wn(t) = mdl_cv2wn.Rsquared.Ordinary
%     plot(cv2Wn(treatment==t& goodAll==1,1),cv2Wn(treatment==t& goodAll==1,2),'.');axis square;hold on
%     title(titles{t});
%     plot([0 2], [0 2]);ylim([.5 2]);xlim([.5 2]);
%     n_cells(t) = sum(treatment==t& goodAll==1)
%     text(.55, 1.75, ['r^2 = ' num2str(rs_cv2wn(t))],'FontSize',18)
%     text(.55, 1.5, ['n = ' num2str(n_cells(t))],'FontSize',18)
% end
% 
% figure
% set(gcf,'Name','cv2 drift');
% for t=1:3
%     subplot(1,3,t)
%     mdl_cv2drift= fitlm(cv2Drift(treatment==t&goodAll==1,1), cv2Drift(treatment==t&goodAll==1,2),2)
%     rs_cv2drift(t) = mdl_cv2drift.Rsquared.Ordinary
%     plot(cv2Drift(treatment==t&goodAll==1,1),cv2Drift(treatment==t&goodAll==1,2),'.');axis square;hold on
%     plot([0 2], [0 2]);ylim([.5 2]);xlim([.5 2]);
%     title(titles{t});
%     n_cells(t) = sum(treatment==t& goodAll==1)
%     text(.55, 1.85, ['r^2 = ' num2str(rs_cv2drift(t))],'FontSize',18)
%     text(.55, 1.65, ['n = ' num2str(n_cells(t))],'FontSize',18)  
% end
% 
% figure
% set(gcf,'Name','cv2 dark');
% for t=1:3
%     subplot(1,3,t)
%     mdl_cv2dark= fitlm(cv2Dark(treatment==t&goodAll==1,1), cv2Dark(treatment==t&goodAll==1,2),2)
%     rs_cv2dark(t) = mdl_cv2dark.Rsquared.Ordinary
%     plot(cv2Dark(treatment==t,1),cv2Dark(treatment==t,2),'.');axis square;hold on
%     plot([0 2], [0 2]);ylim([.5 2]);xlim([.5 2]);
%     title(titles{t});
%     n_cells(t) = sum(treatment==t& goodAll==1)
%     text(.55, 1.85, ['r^2 = ' num2str(rs_cv2dark(t))],'FontSize',18)
%     text(.55, 1.65, ['n = ' num2str(n_cells(t))],'FontSize',18)   
% end


% figure
% for t=1:3
%     subplot(1,3,t)
%     set(gcf,'Name', 'latent')
%     pre = [cycR(find(data & treatment==t),:,1,1) cycR(find(data & treatment==t),:,2,1)];
%     post = [cycR(find(data & treatment==t),:,1,2) cycR(find(data & treatment==t),:,2,2)];
%     [coeff score latent] = pca(pre','algorithm','als','variableweight','variance');
%     plot(latent(1:10)/sum(latent));title(titles{t});
%     
% end
% 
% % how much does each cell contribute to component
% figure
% for t=1:8
%     subplot(2,4,t)
%     set(gcf,'Name', 'coeff')
%     pre = [cycR(find(data_wn & ~inhAll & treatment==t),:,1,1) cycR(find(data_wn & ~inhAll& treatment==t),:,2,1)];
%     post = [cycR(find(data_wn & ~inhAll& treatment==t),:,1,2) cycR(find(data_wn & ~inhAll& treatment==t),:,2,2)];
%     [coeff score latent] = pca(pre','algorithm','als','variableweight','variance');
%     imagesc(coeff(:,1:10));title(titles{t});
% end
% 
% for i = 1:10; %pre stationary then moving
%     subplot(10,1,i);
%     plot(score(:,i))
% end
% 
% %============
% %first two components plotted...stationary then moving
% figure
% plot(score(:,1)); hold on; plot(score(:,2),'r')
% figure
% plot(score(:,1)); hold on; plot(score(:,3),'r')
% 
% preS = pre'*coeff;
% postS = post'*coeff;
% 
% %using 8 components...can change to however many you want to look at
% figure
% for i = 1:8
%     subplot(8,1,i)
%     plot(preS(:,i)); hold on; plot(postS(:,i),'r');
% end
% 
% 
% figure
% plot(preS(1:20,1),preS(1:20,2),'b'); hold on;plot(preS(21:40,1),preS(21:40,2),'b');
% plot(postS(1:20,1),postS(1:20,2),'r'); plot(postS(21:40,1),postS(21:40,2),'r');

%%correlation for wn %%
for t=1:3
    preWnCorr = wnCorr(corrTreatment==t,1)
    postWnCorr=wnCorr(corrTreatment==t,2)
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
    plot(h,'-b');hold on;
    plot(h2,'-r')
    %xlim([-1 1]);
    title(titles{i});
    set(gcf,'Name','Wn Corr Hist')
end

%%correlation for drift%%

clear rsquared
for t=3
    preDriftCorr = driftCorr(corrTreatment==t,1)
    postDriftCorr=driftCorr(corrTreatment==t,2)
    mdl = fitlm(preDriftCorr,postDriftCorr)
    rsquared(t) = mdl.Rsquared.Adjusted
end

titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI','MGluR2', 'MGluR2 + DOI','Lisuride'};
figure
for i = 3
    subplot(1,3,i);
    plot(driftCorr(corrTreatment==i,1),driftCorr(corrTreatment==i,2),'.'); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
    xlabel('correlation coefficient pre'); ylabel('correlation coefficient post')
    text(-.25, .85, ['r^2 = ' num2str(rsquared(i),'%.2f')],'FontSize',20)
    set(gca, 'FontSize',20)
    set(gcf,'Name','Drift Corr')
end

titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI','MGluR2', 'MGluR2 + DOI','Lisuride'};
figure
for i = 1:3
    h= histc(driftCorr(corrTreatment==i,1),-1:.1:1);
    h2=histc(driftCorr(corrTreatment==i,2),-1:.1:1);
    subplot(3,1,i);
    plot(h,'-b');hold on;
    plot(h2,'-r')
    %xlim([-1 1]);
    title(titles{i});
    set(gcf,'Name','Drift Corr Hist')
end

% correlation for darkness
clear rsquared
for t=1:8
    preDarkCorr = darkCorr(corrTreatment==t,1)
    postDarkCorr=darkCorr(corrTreatment==t,2)
    mdl = fitlm(preDarkCorr,postDarkCorr)
    rsquared(t) = mdl.Rsquared.Adjusted
end
titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI','MGlur2','MGlur2 + DOI', 'Lisuride'};
figure
for i = 1:3
    subplot(1,3,i);
    plot(darkCorr(corrTreatment==i,1),darkCorr(corrTreatment==i,2),'.'); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); %title(titles{i});
  %  xlabel('pre dark corr'); ylabel('post')
    set(gca,'FontSize',20)
     text(-.25, .85, ['r^2 = ' num2str(rsquared(i),'%.2f')],'FontSize',20)
    set(gcf,'Name','Dark Corr')
end

titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI','MGluR2', 'MGluR2 + DOI','Lisuride'};
figure
for i = 1:3
    h= histc(darkCorr(corrTreatment==i,1),-1:.1:1);
    h2=histc(darkCorr(corrTreatment==i,2),-1:.1:1);
    subplot(3,1,i);
    plot(h,'-b');hold on;
    plot(h2,'-r')
    %xlim([-1 1]);
    title(titles{i});
    set(gcf,'Name','Drift Corr Hist')
end

%% sta exp var
titles={'Saline','DOI','5HT'};
figure
for t=1:2
subplot(1,2,t)
mdl_expV = fitlm(sta_exp_var(data_wn &(treatment==t),1),sta_exp_var(data_wn&(treatment==t),2))
rsquared_expV(t) = mdl_expV.Rsquared.Ordinary
plot(sta_exp_var(data_wn & treatment==t,1),sta_exp_var(data_wn &treatment==t,2),'.')
%plot(sta_exp_var(useSTA & (treatment==t)',1),sta_exp_var(useSTA &(treatment==t)',2),'.')
hold on; plot([0 1],[0 1]); axis square
title(titles{t}); xlabel('Pre exp var');ylabel('Post exp var');
text(.02, .95, ['r^2 = ' num2str(rsquared_expV(t))],'FontSize',18)
n_cells(t) = sum(data_wn & treatment==t)
text(.02 ,.85, ['n = ' num2str(n_cells(t))],'FontSize',18)
end

%% nx & ny prepost
figure
titles={'Saline','DOI','5HT'};
useSTA = sta_exp_var(:,1)>.6 & sta_exp_var(:,2)>.6
useSTA = useSTA & data_wn'% &(layerAll==4)'
for t=1:2 %5ht doesnt show up with exp var smaller than .4
    figure
   mdl_nx = fitlm(sta_nx(useSTA & (treatment==t)',1),sta_nx(useSTA&(treatment==t)',2))
   mdl_ny= fitlm(sta_ny(useSTA & (treatment==t)',1),sta_ny(useSTA&(treatment==t)',2))
   rsquared_nx(t) = mdl_nx.Rsquared.Ordinary
   rsquared_ny(t) = mdl_ny.Rsquared.Ordinary
%    subplot(2,2,t)
   subplot(1,2,1)
    plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==t)',2),'.','Markersize',10);
    hold on;
    plot(sta_nx(useSTA & (inhAll==1)'&(treatment==t)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==t)',2),'.r','Markersize',10);
    axis square;xlim([0 .7]); ylim([0 .7]); set(gca,'FontSize',10);
    text(.02, .63, ['r^2 = ' num2str(rsquared_nx(t))],'FontSize',10)
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
    text(.02, .55, ['r^2 = ' num2str(rsquared_ny(t))],'FontSize',18)
    mpre=nanmean(mean(sta_ny(useSTA  & (treatment==t)',1)))
    mpost=nanmean(mean(sta_ny(useSTA & (treatment==t)',2)))
    plot(mpre,mpost,'pg','Markersize',10)
    title(titles{t},'FontSize',30); xlabel('Pre ny');ylabel('Post ny');
    n_cells(t) = sum(useSTA &(treatment==t)')
    text(.02 ,.65, ['n = ' num2str(n_cells(t))],'FontSize',18)
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
subplot(3,3,7)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==2|3)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==2|3)',2),'r.','Markersize',10);
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
subplot(3,3,7)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==2|3)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==2|3)',2),'r.','Markersize',10);
xlabel('pre ny');ylabel('post ny')
subplot(3,3,8)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==4)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==4)',2),'r.','Markersize',10);
xlabel('pre ny');ylabel('post ny')
subplot(3,3,9)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==5)',2),'.','Markersize',10);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==5)',2),'r.','Markersize',10);
xlabel('pre ny');ylabel('post ny')

useSTA = sta_exp_var(:,1)>.6 & sta_exp_var(:,2)>.6
useSTA = useSTA & data_wn'% &(layerAll==4)'
figure
titles={'Saline', 'DOI','5HT'}
for t=1:2
subplot(2,3,t)
set(gcf,'Name', 'layer 2/3 prepost nx and ny')
mdl_nxny_pre = fitlm(sta_nx(useSTA & (treatment==t)' &(layerAll==2|3)',1),sta_ny(useSTA&(treatment==t)'&(layerAll==2|3)',1))
rsquared_nxny_pre(t) = mdl_nxny_pre.Rsquared.Ordinary
mdl_nxny_post = fitlm(sta_nx(useSTA & (treatment==t)' &(layerAll==2|3)',2),sta_ny(useSTA&(treatment==t)'&(layerAll==2|3)',2))
rsquared_nxny_post(t) = mdl_nxny_post.Rsquared.Ordinary
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==t)' & (layerAll==2|3)',1),'.','Markersize',10); hold on;
xlabel('nx pre');ylabel('ny pre')
plot(sta_nx(useSTA & (inhAll==1)' & (treatment==t)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==t)' & (layerAll==2|3)',1),'r.','Markersize',10);
hold on;plot([0 1],[0 1]); axis square;xlim([0 1])
text(.02, .8, ['r^2 = ' num2str(rsquared_nxny_pre(t))],'FontSize',18)
subplot(2,3,t+3)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)' & (layerAll==2|3)',2),sta_ny(useSTA&(inhAll==0)'&(treatment==t)' & (layerAll==2|3)',2),'.','Markersize',10); hold on;
xlabel('nx post');ylabel('ny post')
plot(sta_nx(useSTA & (inhAll==1)' & (treatment==t)' & (layerAll==2|3)',2),sta_ny(useSTA&(inhAll==1)'&(treatment==t)' & (layerAll==2|3)',2),'r.','Markersize',10);
plot([0 1],[0 1]);axis square;xlim([0 1])
title(titles{t});
text(.02, .8, ['r^2 = ' num2str(rsquared_nxny_post(t))],'FontSize',18)
end

figure
titles={'Saline', 'DOI','5HT'}
for t=1:2
subplot(2,3,t)
set(gcf,'Name', 'all layers prepost nx and ny')
mdl_nxny_pre = fitlm(sta_nx(useSTA & (treatment==t)',1),sta_ny(useSTA&(treatment==t)',1))
rsquared_nxny_pre(t) = mdl_nxny_pre.Rsquared.Ordinary
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==t)' ,1),'.','Markersize',10);
hold on;plot([0 1],[0 1]); axis square; xlabel('pre nx');ylabel('pre ny');
text(.02, .8, ['r^2 = ' num2str(rsquared_nxny_pre(t))],'FontSize',18)
plot(sta_nx(useSTA & (inhAll==1)' & (treatment==t)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==t)' ,1),'r.','Markersize',10);
subplot(2,3,t+3)
mdl_nxny_post = fitlm(sta_nx(useSTA & (treatment==t)',1),sta_ny(useSTA&(treatment==t)',2))
rsquared_nxny_post(t) = mdl_nxny_post.Rsquared.Ordinary
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)' ,2),sta_ny(useSTA&(inhAll==0)'&(treatment==t)' ,2),'.','Markersize',10);
hold on;plot([0 1],[0 1]);xlim([0 1]);axis square; xlabel('post nx');ylabel('post ny');
plot(sta_nx(useSTA & (inhAll==1)' & (treatment==t)' ,2),sta_ny(useSTA&(inhAll==1)'&(treatment==t)' ,2),'r.','Markersize',10);
text(.02, .8, ['r^2 = ' num2str(rsquared_nxny_post(t))],'FontSize',18)
title(titles{t});
end


useSTA = find((sta_exp_var(:,1)>.6 & sta_exp_var(:,2)>.6) &(treatment==doi)');
for prepost=1:2
    figure
      if prepost==1
        set(gcf,'Name','pre')
    else
         set(gcf,'Name','post');end
    for c=1:length(useSTA)
        subplot(10,7,c)
        imagesc(sta_all_img{useSTA(c),prepost});axis square;
        %imagesc(sta_all_fit{useSTA(c),prepost});axis square
        set(gca,'xticklabel',{[]},'yticklabel',{[]}) 
    end
end

clear useSTA
useSTA = find((sta_exp_var(:,1)>.6 & sta_exp_var(:,2)>.6) &(treatment==saline)' )%&(layerAll==2|3)');
for prepost=1:2
    figure
    if prepost==1
        set(gcf,'Name','pre')
    else
        set(gcf,'Name','post');end
    for c=1:length(useSTA)
        subplot(7,7,c)
        %imagesc(sta_all_img{useSTA(c),prepost});axis square
        imagesc(sta_all_fit{useSTA(c),prepost},[-40 40]);axis square; colormap jet
        set(gca,'xticklabel',{[]},'yticklabel',{[]}) 
    end
end




%%% plot white noise response functions for all units

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
         %   title(titles{l});

        end
    legend ('stationary pre','mv pre','stationary post', 'mv post')
end
end

% titles={'Saline','DOI','5HT'};
% figure
% for t=1:3
%     useN= data_wn & treatment==t
%     MI_crf = (wn_crf(useN,:,:,2)-wn_crf(useN,:,:,1))./(wn_crf(useN,:,:,2)+wn_crf(useN,:,:,1));
%     subplot(1,3,t)
%     h= hist(MI_crf,-1:.1:1);
%     bar(Mbins,h/sum(useN));axis square; title(titles{t}); ylim([0 .5]);
% end

% titles={'Saline','DOI','5HT'};
for t=1:3
    clear useN MI_crf h
    figure
    if t==1, set(gcf, 'Name', 'Saline');
    elseif t==2, set(gcf,'Name','DOI');
    else set(gcf,'Name', '5HT'); end
    
    for l=1:6
        useN= data_wn & treatment==t & layerAll==l
        MI_crf = (wn_crf(useN,:,:,2)-wn_crf(useN,:,:,1))./(wn_crf(useN,:,:,2)+wn_crf(useN,:,:,1));
        subplot(2,3,l)
        Mbins = -1:.2:1
        h= hist(MI_crf,Mbins);
        bar(Mbins,h/sum(useN),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .5]);xlim([-1.5 1.5]); axis square
        title(layerz{l});
    end
end
 
 
for t = 1:2
    useN = find(goodAll & treatment==t& hasWn==1 & (layerAll==2 | layerAll==3) )
    for i = 1:ceil(length((useN))/16)
        figure
        if t==1, set(gcf,'Name','saline wn CRF'),
        elseif t==2, set(gcf,'Name','doi wn CRF'),
        elseif set(gcf, 'Name','ht wn CRF'), %end
        elseif t==3, set(gcf,'Name','ht wn CRF'),
        elseif t==4, set(gcf,'Name','ketanserin wn CRF'),
        elseif t==5, set(gcf,'Name','ketanserin + doi wn CRF'),
        elseif t==6, set(gcf,'Name','MGluR2 wn CRF'),
        elseif t==7, set(gcf, 'Name','MGluR2+ DOI wn CRF'),
        else set(gcf,'Name','Lisuride wn CRF'),
        end
        for j= 1:min(16,length(useN)-(i-1)*16)
            subplot(2,8,j); hold on;
            plot(wn_crf(useN((i-1)*16+j),:,1,1),'Color',[0.5 0 0]);  plot(wn_crf(useN((i-1)*16+j),:,2,1),'Color',[0 0.5 0]); axis square
            plot(wn_crf(useN((i-1)*16+j),:,1,2),'Color',[1 0 0]);  plot(wn_crf(useN((i-1)*16+j),:,2,2),'Color',[0 1 0]);
            yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)]);
            if inhAll(useN((i-1)*16+j)) ==1 , xlabel('inh'); else  xlabel('exc');
            end
        end
    end
end


% % titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
% % figure
% % for t=1:3
% %     set(gcf,'Name', 'spont pre vs post stop')
% %     subplot(1,3,t)
% %     plot(spont(find(data & treatment==t),1,1),spont(find(data& treatment==t),1,2),'.'); title(titles{t}); hold on; axis square; plot([0 50],[0 50]);
% %     xlabel('pre'); ylabel('post');
% %     plot(spont(find(dataInh&inhAll& treatment==t),1,1),spont(find(dataInh & inhAll& treatment==t),1,2),'.r'); hold on; plot([0 50],[0 50]); axis([0 10 0 10])
% % end
% % 
% % figure
% % for t=1:3
% %     set(gcf,'Name', 'spont pre vs post mv')
% %     subplot(1,3,t)
% %     plot(spont(find(data& treatment==t),2,1),spont(find(data& treatment==t),2,2),'.'); title(titles{t}); hold on; plot([0 50],[0 50]); axis square
% %     xlabel('pre'); ylabel('post');
% %     plot(spont(find(dataInh&inhAll& treatment==t),2,1),spont(find(dataInh & inhAll& treatment==t),2,2),'.r'); hold on; plot([0 50],[0 50]); axis([0 15 0 15])
% % end
% % 
% % figure
% % for t=1:3
% %     set(gcf,'Name', 'evoked pre vs post stop')
% %     subplot(1,3,t)
% %     plot(evoked(find(data& treatment==t),1,1),evoked(find(data& treatment==t),1,2),'.'); title(titles{t}); hold on; plot([0 50],[0 50]); axis square
% %     xlabel('pre'); ylabel('post');
% %     plot(evoked(find(dataInh&inhAll& treatment==t),1,1),evoked(find(dataInh & inhAll& treatment==t),1,2),'.r'); hold on; plot([0 50],[0 50]); axis([-5 10 -5 10])
% % end

%%%%%%%%%%%%%%%%%%%%% drift %%%%%%%%%%%%%%%%%%%%%%%%%%   

data = goodAll & useResp' & ~inhAll & hasDrift==1 ;%which cells to include
dataInh =goodAll & useResp' & inhAll & hasDrift==1 ;
dataAll =goodAll & useResp' & hasDrift==1 ;

useN =dataAll;
for c= 1:length(useN)
   % for prepost =1:2
        [respmax oind] = max(mean(drift_orient(c ,:,1,:),4)); %stationary
       % [g h]=(max(respmax(:)));
        prefOri(c)=oind;
   % end
end

useN =dataAll;
for c= 1:length(useN)
   % for prepost =1:2
        [respmax oind] = max(mean(drift_sf(c ,:,1,:),4)); %stationary
       % [g h]=(max(respmax(:)));
        prefSF(c)=oind;
    %end
end

%%% drift_cond_tcourse(cells,move,prepost,orient,sf,time);

for c=1:length(useN)
    tcourse_pref(c,:,:,:) = squeeze(nanmean(drift_cond_tcourse(c,:,:,prefOri(c),:,:),5));
    tcourse_orth(c,:,:,:) = squeeze(nanmean(drift_cond_tcourse(c,:,:,mod(prefOri(c)+3-1,12)+1,:,:),5));    
    tcourse_all(c,:,:,:) = squeeze(mean(nanmean(drift_cond_tcourse(c,:,:,:,:,:),5),4));
end

% figure
% plot(squeeze(tcourse_pref(30,1,1,:)));
% figure
% plot(squeeze(nanmean(drift_cond_tcourse(30,1,1,11,:,:),5))')


% for c=1:length(useN)
%     sf_pref_tcourse(c,:,:,:) = squeeze(nanmean(drift_cond_tcourse(c,:,:,:,prefSF(c),:),4));
% end

clear max amp low
peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont
useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2;



%PSTH for each cell's preferred orientation - stationary

for i= 1:4
    figure
    for t=1:3
        if i==1
            used =goodAll & useResp' & ~inhAll &(hasDrift==1) & (layerAll==5) & treatment==t;  
            set(gcf,'Name','grating lyr5 stationary');
            mn = squeeze(mean(tcourse_pref(used,1,:,:),1))';
            sem = std(tcourse_pref(data & (hasDrift==1) & (layerAll==5) & treatment==t,1,:,:),1)...
                /sqrt(sum(data & (hasDrift==1) & (layerAll==5) & treatment==t))
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==5) & treatment==t);xlim([0 2.5]);
              ylim([0 8]); xlim([0 2.5]);
        elseif i==2
            set(gcf,'Name','grating lyr 2/3 stationary');
            mn = squeeze(mean(tcourse_pref(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,1,:,:),1))';
            sem = std(tcourse_pref(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,1,:,:),1)...
                /sqrt(sum(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t))
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==2| layerAll==3) & treatment==t);xlim([0 2.5]);
             ylim([0 8]); xlim([0 2.5]);
        elseif i ==3
            set(gcf,'Name','grating lyr 4 stationary');
            mn = squeeze(mean(tcourse_pref(data &(hasDrift==1) & layerAll==4 & treatment==t,1,:,:),1))';
            sem = std(tcourse_pref(data &(hasDrift==1) & layerAll==4 & treatment==t,1,:,:),1)...
                /sqrt(sum(data &(hasDrift==1) & layerAll==4 & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==4) & treatment==t);xlim([0 2.5]);
            %ylim([min(min(mn))-1 max(max(mn))+.5]); xlim([0 2.5]);
        elseif i==4
            set(gcf,'Name','grating inh stationary');
            mn = squeeze(mean(tcourse_pref(dataInh & (hasDrift==1) & treatment==t,1,:,:),1))';
            sem = std(tcourse_pref(dataInh &(hasDrift==1) & treatment==t,1,:,:),1)...
                /sqrt(sum(dataInh &(hasDrift==1) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(dataInh & (hasDrift==1) & treatment==t);xlim([0 2.5]);
              ylim([0 16]); xlim([0 2.5]);
        elseif i==5
            set(gcf,'Name','grating lyr6 stationary');
            mn = squeeze(mean(tcourse_pref(data & layerAll==6 & (hasDrift==1) & treatment==t,1,:,:),1))';
            sem = std(tcourse_pref(data & layerAll==6 & (hasDrift==1) & treatment==t,1,:,:),1)...
                /sqrt(sum(data & layerAll==6 & (hasDrift==1) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==6) & treatment==t);xlim([0 2.5]);
            % ylim([min(min(mn))-1 max(max(mn))+.5]); xlim([0 2.5]);
        end
        %         mn = mn - repmat(mn(1,:),[50 1]); %%% subtracts prestim spont from
        %        all time points
        mn = circshift(mn,10);
        sem = circshift(sem,10);
        subplot(1,3,t)
        x = ((1:length(mn)-5)*dt -dt/2);
        plot(x,mn(1:45,1),'LineWidth',2); hold on;set(gca,'fontsize', 18); hold on;
        plot(x,mn(1:45,2),'LineWidth',2);axis square; set(gca,'fontsize', 18); 
        preerr=shadedErrorBar(x,mn(1:45,1),sem(1:45,1)','b',1); axis square; set(gca,'fontsize', 18); hold on;
        posterr=shadedErrorBar(x,mn(1:45,2),sem(1:45,2)','r',1);%axis square; set(gca,'fontsize', 18)
        text(1.75,max(max((mn))+.75), ['n = ' num2str(ncells)],'FontSize',20)
        title([titles{t} ' n = ' num2str(ncells) ]);
        xlabel('time (s)');
        %   ylabel('spikes/sec');
        %ylim([0 max(mn(:))+1])
    end
end

%PSTH for each cell's preferred orientation - moving
for i= 1:5
    figure
    for t=1:3
        if i==1
            set(gcf,'Name','grating lyr5 mv');
            mn = squeeze(nanmean(tcourse_pref(data & (hasDrift==1) & (layerAll==5) & treatment==t,2,:,:),1))';
            sem = nanstd(tcourse_pref(data & (hasDrift==1) & (layerAll==5) & treatment==t,2,:,:),1)...
                /sqrt(sum(data & (hasDrift==1) & (layerAll==5) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==5) & treatment==t); xlim([0 2.5]);
             ylim([-2 7]);
        elseif i==2
            set(gcf,'Name','grating lyr 2/3 mv');
            mn = squeeze(nanmean(tcourse_pref(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,2,:,:),1))';
            sem = nanstd(tcourse_pref(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,2,:,:),1)...
                /sqrt(sum(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==2| layerAll==3) & treatment==t);xlim([0 2.5]);
              ylim([-2 10]);xlim([0 2.5]);
        elseif i ==3
            set(gcf,'Name','grating lyr 4 mv');
            mn = squeeze(nanmean(tcourse_pref(data &(hasDrift==1) & layerAll==4 & treatment==t,2,:,:),1))';
            sem = nanstd(tcourse_pref(data &(hasDrift==1) & layerAll==4 & treatment==t,2,:,:),1)...
                /sqrt(sum(data &(hasDrift==1) & layerAll==4 & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==4) & treatment==t);xlim([0 2.5]);
            %ylim([-.5 4]);xlim([0 2.5]);
        elseif i==4
            set(gcf,'Name','grating inh mv');
            mn = squeeze(nanmean(tcourse_pref(dataInh & (hasDrift==1) & treatment==t,2,:,:),1))';
             sem = nanstd(tcourse_pref(dataInh &(hasDrift==1) & treatment==t,2,:,:),1)...
                /sqrt(sum(dataInh &(hasDrift==1) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(dataInh & (hasDrift==1) & treatment==t);xlim([0 2.5]);
              ylim([-.5 20]);xlim([0 2.5]);
        else i==5
            set(gcf,'Name','grating lyr6 mv');
            mn = squeeze(nanmean(tcourse_pref(data & layerAll==6 & (hasDrift==1) & treatment==t,2,:,:),1))';
            sem = nanstd(tcourse_pref(data & layerAll==6 & (hasDrift==1) & treatment==t,2,:,:),1)...
                /sqrt(sum(data & layerAll==6 & (hasDrift==1) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==6) & treatment==t);xlim([0 2.5]);
            %  ylim([-.5 12]);xlim([0 2.5]);
        end
        mn = mn - repmat(mn(1,:),[50 1]); %%% subtracts prestim spont from
        %        all time points
        mn = circshift(mn,10);
        sem = circshift(sem,10);
        subplot(1,3,t)
        x = ((1:length(mn)-5)*dt -dt/2);
        plot(x,mn(1:45,1),'LineWidth',2); hold on;set(gca,'fontsize', 18); hold on;
        plot(x,mn(1:45,2),'LineWidth',2);axis square; set(gca,'fontsize', 18);
        preerr=shadedErrorBar(x,mn(1:45,1),(sem(1:45,1))','b',1); axis square; set(gca,'fontsize', 18); hold on;
        posterr=shadedErrorBar(x,mn(1:45,2),(sem(1:45,2))','r',1);%axis square; set(gca,'fontsize', 18)
        text(1.75,max(max((mn))+.75), ['n = ' num2str(ncells)],'FontSize',20)
        title(titles{t});
        xlabel('time (s)');
        %   ylabel('spikes/sec');
        %ylim([0 max(mn(:))+1])
    end
end



% prefOri_pre = prefOri(1,:)'
% prefOri_post = prefOri(2,:)'
% 
% mn_pref = squeeze(nanmean(drift_cond_tcourse,5)); %average over all 6 SFs
% 
% %ave PSTH for each orientation
% for t=1:3
%     figure
%     clear mn_pref_psth
%     for ori=1:12
%         pref = prefOri_pre==ori;
%         mn_pref_psth = squeeze(nanmean(mn_pref(data & treatment==t & layerAll==5 &  (pref==1)',:,:,ori,:),1));
% %         mn_pref_psth = mn_pref_psth - repmat(mn_pref_psth(1,1,:),[50 1]);%subtract prestim activity
%         mn_pref_psth = circshift(mn_pref_psth,10);
%         subplot(3,4,ori)
%         plot((1:length(mn_pref_psth)-5)*dt-dt/2,squeeze(mn_pref_psth(1,:,1:45)),'LineWidth',2);axis square; set(gca,'fontsize', 18);
%         axis square; set(gca,'fontsize',18);
%         title(titles{t}); xlabel('time(s)')
%         ylabel('spikes/sec');
%         mn_pref_all(:,:,:,ori,t)= mn_pref_psth
%     end
% end


% test2 =repmat(mn_pref_psth(1,1,:),[1 50]);

% figure
% for t=1:3
%     subplot(1,3,t)
%     R = normrnd(0,.2,1,sum(data&treatment==t))
%     plot(prefOri(1,data & treatment==t)+R,prefOri(2,data & treatment==t)+R,'.'); axis square; xlim([0 13]);ylim([0 13]);
%     hold on;  plot([0 13],[0 13])
% end
% %+ normrand same size as number of points w/mean 0 std .2
% %or heatmap
% 
% meanOri = round(mean(prefOri(:,:),1))
% titles ={'Saline','DOI','5HT'};
% figure
% for t=1:3
%     useOri_pre= prefOri_pre(treatment==t)
%     useOri_post= prefOri_post(treatment==t)
%     %     useOri_post= find(prefOri_post & (treatment==t)')
%     subplot(2,3,t)
%     hist(useOri_pre);xlabel('pre orientation pref'); ylabel('proportion of cells')
%     xlim([0 12]); ylim([0 .03]);
%     Mbins= 0:1:12
%     hpre= hist(useOri_pre,0:1:12);
%     bar(Mbins,hpre/sum(useOri_pre)); xlim([0 13]);ylim([0 .035]);axis square;
%     title(titles{t});
%     subplot(2,3,t+3)
%     hpost=hist(useOri_post,0:1:12);
%     bar(Mbins,hpost/sum(useOri_post))
%     xlabel('post orientation preference');ylabel('proportion of cells')
%     axis square
%     xlim([0 13]);ylim([0 .035]);
% %     [p h]=ranksum(prefOri_pre(treatment==t),prefOri_post(treatment==t))
% %     text(6, .0275, ['p = ' num2str(p)],'FontSize',18)
% end

useOsi_pre = drift_osi(:,1,1)>0
useOsi_post = drift_osi(:,1,2)>0
useOsi=(useOsi_pre & useOsi_post)'


useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2;
data = goodAll & useResp' & ~inhAll & hasDrift==1 & useOsi==1  ;%which cells to include
dataInh =goodAll & useResp' & inhAll & hasDrift==1 & useOsi==1;
dataAll =goodAll & useResp' & hasDrift==1 &useOsi==1;
for mv=1:2
clear R rs_pf
figure
if mv==1, set(gcf,'Name', 'stationary');
else
set(gcf,'Name', 'mv');
end
for t=1:3
%  R = normrnd(0,.2,sum(data&treatment==t),1)
%  R2 = normrnd(0,.2,sum(data&treatment==t),1)
%  R3 = normrnd(0,.2,sum(dataInh&treatment==t),1)
%  R4 = normrnd(0,.2,sum(dataInh&treatment==t),1)
 R = 0
 R2 = 0
 R3 = 0
 R4 = 0
 pre=drift_pref_theta(dataAll&treatment==t,mv,1); post= drift_pref_theta(dataAll&treatment==t,mv,2)
post(pre-post > pi/2) = post(pre-post>pi/2)+pi; pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;

 %mdl_pf= fitlm(pre,post)
 c=corrcoef(pre,post)
 rs_pf(t) = c(2,1);
 
 [r m b] =regression(pre',post')
 rs_pf(t) = r;
 %rs_pf(t) = mdl_pf.Rsquared.Ordinary
 subplot(1,3,t)
  pre=drift_pref_theta(data&treatment==t,mv,1); post= drift_pref_theta(data&treatment==t,mv,2)
post(pre-post > pi/2) = post(pre-post>pi/2)+pi; pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;
 plot(pre+R,post+R2,'.','Markersize',18)
 hold on; axis square;plot([0 3.25], [0 3.25]);xlim([0 3.25]);ylim([0 3.25]);
 
  pre=drift_pref_theta(dataInh&treatment==t,mv,1); post= drift_pref_theta(dataInh&treatment==t,mv,2)
post(pre-post > pi/2) = post(pre-post>pi/2)+pi; pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi; 
 plot(pre+R3,post+R4,'r.','Markersize',18)
 text(.2, 3, ['r^2 = ' num2str(rs_pf(t),'%.2f')],'FontSize',18)
  n_cells(t) = sum(dataAll&treatment==t) %% not calculating
%  correctly...showing fewer units on plot
  text(.2, 2.7, ['n = ' num2str(n_cells(t))],'FontSize',20)
 %set(gca,'xtick',0:pi/2:pi,'xticklabel',1:180,'Fontsize',14)
 set(gca,'xtick',([pi/4 pi/2 .75*pi pi]))
 % set(gca,'xticklabels',({'\pi/4','\pi/2','3/4\pi','\pi'}),'FontSize',20)
 set(gca,'xticklabels',({'45','90','135','180'}),'FontSize',20)
 set(gca,'ytick',([0 pi/4 pi/2 .75*pi pi]))
 %set(gca,'yticklabels',({'0','\pi/4','\pi/2','3/4\pi','\pi'}),'FontSize',20)
 set(gca,'yticklabels',({'0','45','90','135','180'}),'FontSize',20)
end
end

% figure
% for t=1:3
%     useN=dataAll & treatment==t & hasDrift==1
%     miDrift= (drift_pref_theta(useN &(layerAll==2|layerAll==3),1,2))-(drift_pref_theta(useN&(layerAll==2|layerAll==3),1,1))...
%         ./drift_pref_theta(useN&(layerAll==2|layerAll==3),1,2)+ (drift_pref_theta(useN&(layerAll==2|layerAll==3),1,1));
%     h= hist(miDrift,-1:.1:1);
%     Mbins=-1:.1:1
%     % hist(drift_pref_theta(dataAll & treatment==2,1,1))
%     subplot(1,3,t)
%     bar(Mbins,h/sum(useN &(layerAll==2|layerAll==3)),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .4]); xlim([-1.5 1.5]);axis square; title('layer 2/3')
%     
% end
% 
% figure
% for t=1:3
%     clear useN
%     useN = data & treatment==t
%     for c = 1:sum(useN)
%         mn = squeeze(mean(drift_cond_tcourse(useN(c),1,:,meanOri(useN(c),:,:))))
%         % mn = squeeze(drift_cond_tcourse(useN(c),1,:,meanOri(useN(c)),:,:))
%         mn = mn - repmat(mn(1,:),[50 1]);
%         mn = circshift(mn,10);
%         subplot(1,3,t)
%         %    if prepost==1
%         plot((1:length(mn)-5)*dt -dt/2,mn(1:45,:),'b');axis square; set(gca,'fontsize', 18);hold on
%         % else
%         %   plot((1:length(mn)-5)*dt -dt/2,mn(1:45,:),'r')
%         %end
%     end
% end
useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2;
data = goodAll & useResp' & ~inhAll & hasDrift==1 ;%which cells to include
dataInh =goodAll & useResp' & inhAll & hasDrift==1 ;
dataAll =goodAll & useResp' & hasDrift==1 ;

titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
dt = 0.05;
for i= 1:4
    figure
    for t=1:3
        if i==1
            set(gcf,'Name','grating lyr5');
            mn = squeeze(mean(mnPsth(data & (hasDrift==1) & (layerAll==5) & treatment==t,:,:),1))';
            sem = std(mnPsth(data & (hasDrift==1) & (layerAll==5) & treatment==t,:,:),1)...
                /sqrt(sum(data & (hasDrift==1) & (layerAll==5) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==5) & treatment==t);xlim([0 2.5]);
            %    ylim([-1 1]); xlim([0 2.5]);
        elseif i==2
            set(gcf,'Name','grating lyr 2/3');
            mn = squeeze(mean(mnPsth(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,:,:),1))';
            sem = std(mnPsth(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,:,:),1)...
                /sqrt(sum(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==2| layerAll==3) & treatment==t);xlim([0 2.5]);
            %  ylim([-.5 5]);xlim([0 2.5]);
        elseif i ==3
            set(gcf,'Name','grating lyr 4');
            mn = squeeze(mean(mnPsth(data &(hasDrift==1) & layerAll==4 & treatment==t,:,:),1))';
            sem = std(mnPsth(data &(hasDrift==1) & layerAll==4 & treatment==t,:,:),1)...
                /sqrt(sum(data &(hasDrift==1) & layerAll==4 & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==4) & treatment==t);xlim([0 2.5]);
            %  ylim([-.5 4]);xlim([0 2.5]);
        elseif i==4
            set(gcf,'Name','grating inh');
            mn = squeeze(mean(mnPsth(dataInh & (hasDrift==1) & treatment==t,:,:),1))';
            sem = std(mnPsth(dataInh &(hasDrift==1) & treatment==t,:,:),1)...
                /sqrt(sum(dataInh &(hasDrift==1) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(dataInh & (hasDrift==1) & treatment==t);xlim([0 2.5]);
            %  ylim([-.5 12]);xlim([0 2.5]);
        elseif i==5
            set(gcf,'Name','grating lyr6');
            mn = squeeze(mean(mnPsth(data & layerAll==6 & (hasDrift==1) & treatment==t,:,:),1))';
            sem = std(mnPsth(data & layerAll==6 & (hasDrift==1) & treatment==t,:,:),1)...
                /sqrt(sum(data & layerAll==6 & (hasDrift==1) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==6) & treatment==t);xlim([0 2.5]);
          %  ylim([-.5 12]);xlim([0 2.5]);
        end
    %   mn = mn - repmat(mn(1,:),[50 1]); %%% subtracts prestim spont
        mn = circshift(mn,10);
        sem = circshift(sem,10);
        subplot(1,3,t)
        
        x = ((1:length(mn)-5)*dt -dt/2);
        plot(x,mn(1:45,1),'LineWidth',2); hold on;set(gca,'fontsize', 18); hold on;
        plot(x,mn(1:45,2),'LineWidth',2);axis square; set(gca,'fontsize', 18); 
        preerr=shadedErrorBar(x,mn(1:45,1),sem(1:45,1)','b',1); axis square; set(gca,'fontsize', 18); hold on;
        posterr=shadedErrorBar(x,mn(1:45,2),sem(1:45,2)','r',1);%axis square; set(gca,'fontsize', 18)
        text(1.75,max(max((mn))+.75), ['n = ' num2str(ncells)],'FontSize',20);
        
        %plot((1:length(mn)-5)*dt -dt/2,mn(1:45,:),'LineWidth',2);axis square; set(gca,'fontsize', 18);
        title(titles{t}); xlabel('time (s)');
        ylabel('spikes/sec');
        %ylim([0 max(mn(:))+1])
    end
end

mnPsth_mv = squeeze(nanmean(nanmean(drift_cond_tcourse,5),4))
titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
dt = 0.05;
for i= 1:4
    figure
    for t=1:3
        if i==1
            set(gcf,'Name','grating lyr5 stationary');
            mn = squeeze(mean(mnPsth_mv(data & (hasDrift==1) & (layerAll==5) & treatment==t,1,:,:),1))';
            sem = nanstd(mnPsth_mv(data & (hasDrift==1) & (layerAll==5) & treatment==t,1,:,:),1)...
                /sqrt(sum(data & (hasDrift==1) & (layerAll==5) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==5) & treatment==t);xlim([0 2.5]);
            % ylim([-1 1]); xlim([0 2.5]);
        elseif i==2
            set(gcf,'Name','grating lyr 2/3 stationary');
            mn = squeeze(mean(mnPsth_mv(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,1,:,:),1))';
            sem = nanstd(mnPsth_mv(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,1,:,:),1)...
                /sqrt(sum(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==2| layerAll==3) & treatment==t);xlim([0 2.5]);
            %  ylim([-.5 5]);xlim([0 2.5]);
        elseif i ==3
            set(gcf,'Name','grating lyr 4 stationary');
            mn = squeeze(mean(mnPsth_mv(data &(hasDrift==1) & layerAll==4 & treatment==t,1,:,:),1))';
            sem = nanstd(mnPsth_mv(data &(hasDrift==1) & layerAll==4 & treatment==t,1,:,:),1)...
                /sqrt(sum(data &(hasDrift==1) & layerAll==4 & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==4) & treatment==t);xlim([0 2.5]);
            %  ylim([-.5 4]);xlim([0 2.5]);
        elseif i==4
            set(gcf,'Name','grating inh stationary');
            mn = squeeze(mean(mnPsth_mv(dataInh & (hasDrift==1) & treatment==t,1,:,:),1))';
            sem = nanstd(mnPsth_mv(dataInh &(hasDrift==1) & treatment==t,1,:,:),1)...
                /sqrt(sum(dataInh &(hasDrift==1) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(dataInh & (hasDrift==1) & treatment==t);xlim([0 2.5]);
            %  ylim([-.5 12]);xlim([0 2.5]);
        elseif i==5
            set(gcf,'Name','grating lyr6 stationary');
            mn = squeeze(mean(mnPsth_mv(data & layerAll==6 & (hasDrift==1) & treatment==t,1,:,:),1))';
            sem = nanstd(mnPsth_mv(data & layerAll==6 & (hasDrift==1) & treatment==t,1,:,:),1)...
                /sqrt(sum(data & layerAll==6 & (hasDrift==1) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==6) & treatment==t);xlim([0 2.5]);
            %  ylim([-.5 12]);xlim([0 2.5]);
        end
        %        mn = mn - repmat(mn(1,:),[50 1]); %%% subtracts prestim spont from
        %        all time points
        mn = circshift(mn,10);
        sem = circshift(sem,10);
        subplot(1,3,t)
        %  plot((1:length(mn)-5)*dt -dt/2,mn(1:45,:),'LineWidth',2);axis square; set(gca,'fontsize', 18);
        x = ((1:length(mn)-5)*dt -dt/2);
        plot(x,mn(1:45,1),'LineWidth',2); hold on;set(gca,'fontsize', 18); hold on;
        plot(x,mn(1:45,2),'LineWidth',2);axis square; set(gca,'fontsize', 18);
        preerr=shadedErrorBar(x,mn(1:45,1),sem(1:45,1)','b',1); axis square; set(gca,'fontsize', 18); hold on;
        posterr=shadedErrorBar(x,mn(1:45,2),sem(1:45,2)','r',1);%axis square; set(gca,'fontsize', 18)
        text(1.75,max(max((mn))+.75), ['n = ' num2str(ncells)],'FontSize',20)
        title(titles{t});
        xlabel('time (s)');
        %   ylabel('spikes/sec');
        %ylim([0 max(mn(:))+1])
    end
end


mnPsth_mv = squeeze(nanmean(nanmean(drift_cond_tcourse,5),4))
titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
dt = 0.05;
for i= 1:4
    figure
    for t=1:3
        if i==1
            set(gcf,'Name','grating lyr5 moving');
            mn = squeeze(mean(mnPsth_mv(data & (hasDrift==1) & (layerAll==5) & treatment==t,2,:,:),1))';
            sem = nanstd(mnPsth_mv(data & (hasDrift==1) & (layerAll==5) & treatment==t,2,:,:),1)...
                /sqrt(sum(data & (hasDrift==1) & (layerAll==5) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==5) & treatment==t);xlim([0 2.5]);
            % ylim([-1 1]); xlim([0 2.5]);
        elseif i==2
            set(gcf,'Name','grating lyr 2/3 moving');
            mn = squeeze(mean(mnPsth_mv(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,2,:,:),1))';
            sem = nanstd(mnPsth_mv(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,2,:,:),1)...
                /sqrt(sum(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==2| layerAll==3) & treatment==t);xlim([0 2.5]);
            %  ylim([-.5 5]);xlim([0 2.5]);
        elseif i ==3
            set(gcf,'Name','grating lyr 4 moving');
            mn = squeeze(mean(mnPsth_mv(data &(hasDrift==1) & layerAll==4 & treatment==t,2,:,:),1))';
            sem = nanstd(mnPsth_mv(data &(hasDrift==1) & layerAll==4 & treatment==t,2,:,:),1)...
                /sqrt(sum(data &(hasDrift==1) & layerAll==4 & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==4) & treatment==t);xlim([0 2.5]);
            %  ylim([-.5 4]);xlim([0 2.5]);
        elseif i==4
            set(gcf,'Name','grating inh moving');
            mn = squeeze(mean(mnPsth_mv(dataInh & (hasDrift==1) & treatment==t,2,:,:),1))';
            sem = nanstd(mnPsth_mv(dataInh &(hasDrift==1) & treatment==t,2,:,:),1)...
                /sqrt(sum(dataInh &(hasDrift==1) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(dataInh & (hasDrift==1) & treatment==t);xlim([0 2.5]);
            %  ylim([-.5 12]);xlim([0 2.5]);
        elseif i==5
            set(gcf,'Name','grating lyr6 moving');
            mn = squeeze(mean(mnPsth_mv(data & layerAll==6 & (hasDrift==1) & treatment==t,2,:,:),1))';
            sem = nanstd(mnPsth_mv(data & layerAll==6 & (hasDrift==1) & treatment==t,2,:,:),1)...
                /sqrt(sum(data & layerAll==6 & (hasDrift==1) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==6) & treatment==t);xlim([0 2.5]);
            %  ylim([-.5 12]);xlim([0 2.5]);
        end
        %        mn = mn - repmat(mn(1,:),[50 1]); %%% subtracts prestim spont from
        %        all time points
        mn = circshift(mn,10);
        sem = circshift(sem,10);
        subplot(1,3,t)
        %plot((1:length(mn)-5)*dt -dt/2,mn(1:45,:),'LineWidth',2);axis square; set(gca,'fontsize', 18);
        x = ((1:length(mn)-5)*dt -dt/2);
        plot(x,mn(1:45,1),'LineWidth',2); hold on;set(gca,'fontsize', 18); hold on;
        plot(x,mn(1:45,2),'LineWidth',2);axis square; set(gca,'fontsize', 18);
        preerr=shadedErrorBar(x,mn(1:45,1),sem(1:45,1)','b',1); axis square; set(gca,'fontsize', 18); hold on;
        posterr=shadedErrorBar(x,mn(1:45,2),sem(1:45,2)','r',1);%axis square; set(gca,'fontsize', 18)
        text(1.75,max(max((mn))+.75), ['n = ' num2str(ncells)],'FontSize',20);
        title(titles{t});
        xlabel('time (s)');
        %   ylabel('spikes/sec');
        %ylim([0 max(mn(:))+1])
    end
end



% titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
% dt = 0.05;
% for i= 1:5
%     figure
%     for t=1:3
%         if i==1
%             normA = mnPsth(data & layerAll==5 & treatment==t,:,:) - min(mnPsth(data & layerAll==5 & treatment==t,:,:),1)
%             normA = normA ./ max(normA,1)  
%             mn=squeeze(mean(normA))'  
%             set(gcf,'Name','grating lyr5');
%             %mn = squeeze(mean(mnPsth(data & (hasDrift==1) & (layerAll==5) & treatment==t,:,:),1))';
%             ylim([0 1]); xlim([0 2.5]);
%         elseif i==2
%             set(gcf,'Name','grating lyr 2/3');
%             normA = mnPsth(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,:,:);
%             normA = normA ./ max(normA,1)
%             mn=squeeze(mean(normA))'
%             %mn = squeeze(mean(mnPsth(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,:,:),1))';
%             ylim([0 1]);xlim([0 2.5]);
%         elseif i ==3
%             normA = mnPsth(data &(hasDrift==1) & layerAll==4 & treatment==t,:,:);
%             normA = normA ./ max(normA,1)
%             mn=squeeze(mean(normA))'
%             set(gcf,'Name','grating lyr 4');
% %             mn = squeeze(mean(mnPsth(data &(hasDrift==1) & layerAll==4 & treatment==t,:,:),1));
%             ylim([0 1]);xlim([0 2.5]);
%         elseif i==4
%             normA = mnPsth(dataInh & (hasDrift==1) & treatment==t,:,:);
%             normA = normA ./ max(normA,1)
%             mn=squeeze(mean(normA))'
%             set(gcf,'Name','grating inh');
%            % mn = squeeze(mean(mnPsth(dataInh & (hasDrift==1) & treatment==t,:,:),1))';
%             ylim([0 1]);xlim([0 2.5]);
%         elseif i==5
%             normA = mnPsth(data & layerAll==6 & (hasDrift==1) & treatment==t,:,:);
%             set(gcf,'Name','grating lyr6');
%             %mn = squeeze(mean(mnPsth(data & layerAll==6 & (hasDrift==1) & treatment==t,:,:),1))';
%             ylim([0 1]);xlim([0 2.5]);
%         end
%        mn = mn - repmat(mn(1,:),[50 1]); %%% subtracts prestim spont
%         mn = circshift(mn,10);
%         subplot(1,3,t)
%         plot((1:length(mn)-5)*dt -dt/2,mn(1:45,:),'LineWidth',2);axis square; set(gca,'fontsize', 18);
%      
%         title(titles{t}); xlabel('time (s)');
%         ylabel('spikes/sec');
%         %ylim([0 max(mn(:))+1])
%     end
% end



for l=1:6
figure
useN= find(data & hasDrift & treatment==ht &layerAll==l)
length(useN)
for c=1:length(useN)/2
   % set (gcf,'Name','Lyr 2/3 psth/unit') % set up in loop
    subplot(11,11,c)
    plot(squeeze(mnPsth(useN(c),1,:)));hold on;
    plot(squeeze(mnPsth(useN(c),2,:)),'r')
end
end


titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
figure
for t=1:3
    mnpre = squeeze(mean(mnPsth(data & layerAll==5 & treatment==t,1,:),1))';
    mnpost = squeeze(mean(mnPsth(data & layerAll==5 &treatment==t,2,:),1))';
    mnpre = mnpre - mnpre(1);
    mnpre = circshift(mnpre,10);
    mnpost = mnpost - mnpost(1);
    scalefact(t)=mean(mnpre)./mean(mnpost)
    mnpost_scale = scalefact(t)*(circshift(mnpost,10));
    mnpost=circshift(mnpost,10);
    subplot(1,3,t)
    plot((1:length(mnpre)-5)*dt -dt/2,mnpre(:,1:45));axis square; title(titles{t}); xlabel('secs'); ylabel('sp/sec'); hold on;
    plot((1:length(mnpost)-5)*dt -dt/2,mnpost(:,1:45),'r')
    plot((1:length(mnpost_scale)-5)*dt -dt/2,mnpost_scale(:,1:45),'g'); ylim([0 4.5]);
    text(.02, .63, ['scale factor = ' num2str(scalefact(t))],'FontSize',12)
end

useOsi_pre = drift_osi(:,1,1)>0
useOsi_post = drift_osi(:,1,2)>0
useOsi=(useOsi_pre & useOsi_post)'
useResp = peakresp(:,1,1)>2 & peakresp(:,1,2)>2;
% data = goodAll & useResp' &  useOsi==1 ;%which cells to include
dataAll = goodAll & useResp' & hasDrift==1 & useOsi ==1 ;

titles = {'Saline','DOI','5HT'};
%osi all layers
for mv=1:2
    if mv ==1
        useResp = peakresp(:,1,1)>2 & peakresp(:,1,2)>2;
    else
        useResp = peakresp(:,2,1)>2 & peakresp(:,2,2)>2;
    end
    dataAll = goodAll & useResp' & hasDrift==1 & useOsi ==1 ;
    figure
if mv==1 set(gcf,'Name','stationary mean OSI'), else set(gcf,'Name','mv mean OSI'), end;
for t=1:3
    osi(1,:) = nanmean(abs(drift_osi(dataAll&treatment==t,mv,1)))
    osierr(1,:) = squeeze(nanstd(drift_osi(dataAll & treatment==t,mv,1)/sqrt(sum(drift_osi(dataAll & treatment==t,mv,1),1))));
    osi(2,:) = nanmean(abs(drift_osi(dataAll&treatment==t,mv,2)))
    osierr(2,:) = squeeze(nanstd(drift_osi(dataAll & treatment==t,mv,2)/sqrt(sum(drift_osi(dataAll & treatment==t,mv,2),1))));
    subplot(1,3,t)
    barweb(osi,osierr); ylim([0 1]); axis square
%     [p h]=ranksum(drift_osi(dataAll&treatment==t,mv,1),drift_osi(dataAll&treatment==t,mv,2))
%     text(.55, .9, ['p = ' num2str(p,'%.2f')],'FontSize',18)
    n_cells(t) = sum(treatment==t& dataAll)
    text(.55, .8, ['n = ' num2str(n_cells(t))],'FontSize',18)
    title(titles{t})
    set(gca,'xtick',[1 2])
    set(gca,'xticklabel',['Pre      Post'])
end
end

for mv=1:2
figure
  if mv ==1
        useResp = peakresp(:,1,1)>2 & peakresp(:,1,2)>2; set(gcf,'Name','stationary OSI');
    else
        useResp = peakresp(:,2,1)>2 & peakresp(:,2,2)>2;
       set(gcf,'Name','mv OSI'), end
       dataAll = goodAll & useResp' & hasDrift==1 & useOsi ==1 ;

for t=1:3
    subplot(1,3,t)
    mdl_osi= fitlm(drift_osi(dataAll&treatment==t & hasDrift==1,mv,1),drift_osi(dataAll&treatment==t &hasDrift==1,mv,2),2)
    rs_osi(t) = mdl_osi.Rsquared.Ordinary
    plot(drift_osi(dataAll&treatment==t & ~inhAll&hasDrift==1,mv,1),drift_osi(dataAll&treatment==t & ~inhAll&hasDrift==1,mv,2),'.','Markersize',10);hold on;
    plot(drift_osi(dataAll&treatment==t & inhAll==1&hasDrift==1,mv,1),drift_osi(dataAll&treatment==t & inhAll==1&hasDrift==1,mv,2),'r.','Markersize',10);
    axis square;hold on; plot([0 1],[0 1]);
   % n_cells(t) = sum(treatment==t& dataAll&hasDrift==1)
    text(.1, .9, ['r^2 = ' num2str(rs_osi(t))],'FontSize',18)
   % text(.1, .8, ['n = ' num2str(n_cells(t))],'FontSize',18) 
    title(titles{t})
end
end


%osi specific layers
figure
for t=1:3
%     if mv ==1
%         useResp = peakresp(:,1,1)>2 & peakresp(:,1,2)>2;
%     else
%         useResp = peakresp(:,2,1)>2 & peakresp(:,2,2)>2;
%     end
    dataAll = goodAll & useResp' & hasDrift==1 & useOsi ==1 ;
    figure
if t == 1, set(gcf,'Name','Saline')
elseif t==2, set(gcf,'Name','DOI')
else set(gcf,'Name','5HT'), end
    for l=1:6
    osi(1,:) = nanmean(abs(drift_osi(data&treatment==t&layerAll==l,1,1)))
    osierr(1,:) = squeeze(nanstd(drift_osi(data & treatment==t&layerAll==l,1,1)/sqrt(sum(drift_osi(data & treatment==t&layerAll==l,1,1),1))));
    osi(2,:) = nanmean(abs(drift_osi(data&treatment==t&layerAll==l,1,2)))
    osierr(2,:) = squeeze(nanstd(drift_osi(data & treatment==t&layerAll==l,1,2)/sqrt(sum(drift_osi(data & treatment==t&layerAll==l,1,2),1))));
    subplot(2,3,l)
    barweb(osi,osierr); ylim([0 1]); axis square
%     [p h]=ranksum(drift_osi(data&treatment==t&layerAll==l,1,1),drift_osi(data&treatment==t&layerAll==l,1,2))
%     text(.9, 1, ['p = ' num2str(p)],'FontSize',18)
    end
end

%%% plot orientation tuning curves for all units
for t = 2
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
    for i = 1:length(useN)/4
        np = ceil(sqrt(length(useN)/4));
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

titlesInh = {'Saline Inh','DOI Inh','5HT Inh'};
for mv=1:2
figure
 if mv==1, set(gcf,'Name', 'mean SF stationary')
    else set(gcf,'Name', 'mean SF mv'), end
for t=1:3
SFpre= nanmean(drift_sf(data & treatment==t & hasDrift==1,2:7,mv,1),1)
stderrorPre = nanstd(drift_sf(data & treatment==t & hasDrift==1,2:7,mv,1)) / sqrt(sum(data & treatment==t & hasDrift==1));
SFpost = nanmean(drift_sf(data & treatment==t & hasDrift==1,2:7,mv,2),1)
stderrorPost = nanstd(drift_sf(data & treatment==t & hasDrift==1,2:7,mv,2)) / sqrt(sum(data & treatment==t & hasDrift==1));
subplot(2,3,t)
errorbar(SFpre,stderrorPre) %'-^','MarkerEdgeColor','r','Color','r','LineWidth', 2);
hold on; 
errorbar(SFpost,stderrorPost)
%plot(SFpost,'r','LineWidth',2);axis square
title(titles{t}); xlabel('SF'); ylabel('sp/sec');xlim([1 6]); set(gca, 'XTickLabels', [2 3 4 5 6 7]);
SFpreInh= nanmean(drift_sf(dataInh & treatment==t & hasDrift==1,2:7,mv,1),1)
stderrorPreInh = nanstd(drift_sf(dataInh & treatment==t & hasDrift==1,2:7,mv,1)) / sqrt(sum(dataInh & treatment==t & hasDrift==1));

SFpostInh = nanmean(drift_sf(dataInh & treatment==t & hasDrift==1,2:7,mv,2),1)
stderrorPostInh = nanstd(drift_sf(dataInh & treatment==t & hasDrift==1,2:7,mv,2)) / sqrt(sum(dataInh & treatment==t & hasDrift==1));

subplot(2,3,t+3)
errorbar(SFpreInh,stderrorPreInh)%'LineWidth',2);
hold on; errorbar(SFpostInh,stderrorPostInh)%'r','LineWidth',2);
axis square; xlim([1 6]); set(gca, 'XTickLabels', [2 3 4 5 6 7]);
title(titlesInh{t}); xlabel('SF'); ylabel('sp/sec');
end
end
%%%%

% titlesInh = {'Saline Inh','DOI Inh','5HT Inh'};
% for mv=1:2
% figure
%  if mv==1, set(gcf,'Name', 'mean LOW SF stationary')
%     else set(gcf,'Name', 'mean LOW SF mv'), end
% for t=1:3
% lowSFpre= nanmean(drift_sf(data & treatment==t & hasDrift==1,2:4,mv,1),1)
% lowSFpost = nanmean(drift_sf(data & treatment==t & hasDrift==1,2:4,mv,2),1)
% subplot(2,3,t)
% plot(lowSFpre,'LineWidth',2);hold on; plot(lowSFpost,'r','LineWidth',2);axis square
% title(titles{t}); xlabel('SF'); ylabel('sp/sec');set(gca, 'XTickLabels', [2 3 4]);
% lowSFpreInh= nanmean(drift_sf(dataInh & treatment==t & hasDrift==1,2:4,mv,1),1)
% lowSFpostInh = nanmean(drift_sf(dataInh & treatment==t & hasDrift==1,2:4,mv,2),1)
% subplot(2,3,t+3)
% plot(lowSFpreInh,'LineWidth',2);hold on; plot(lowSFpostInh,'r','LineWidth',2);axis square
% title(titlesInh{t}); xlabel('SF'); ylabel('sp/sec');set(gca, 'XTickLabels', [2 3 4]);
% end
% end
% 
% titlesInh = {'Saline Inh','DOI Inh','5HT Inh'};
% for mv=1:2
% figure
%  if mv==1, set(gcf,'Name', 'mean High SF stationary')
%     else set(gcf,'Name', 'mean High SF mv'), end
% for t=1:3
% highSFpre= nanmean(drift_sf(data & treatment==t & hasDrift==1,5:7,mv,1),1)
% highSFpost = nanmean(drift_sf(data & treatment==t & hasDrift==1,5:7,mv,2),1)
% subplot(2,3,t)
% plot(highSFpre,'LineWidth',2);hold on; plot(highSFpost,'r','LineWidth',2);axis square
% title(titles{t}); xlabel('SF'); ylabel('sp/sec');set(gca, 'XTickLabels', [5 6 7]);
% highSFpreInh= nanmean(drift_sf(dataInh & treatment==t & hasDrift==1,4:6,mv,1),1)
% highSFpostInh = nanmean(drift_sf(dataInh & treatment==t & hasDrift==1,4:6,mv,2),1)
% subplot(2,3,t+3)
% plot(highSFpreInh,'LineWidth',2);hold on; plot(highSFpostInh,'r','LineWidth',2);axis square
% title(titlesInh{t}); xlabel('SF'); ylabel('sp/sec');set(gca, 'XTickLabels', [5 6 7]);
% 
% end
end

%plot spatial frequency tuning curves for all units
for t = 2
    figure
    if t==1, set(gcf,'Name','saline SF'),
    elseif t==2, set(gcf,'Name','doi SF'),
    elseif t==3, set(gcf,'Name','ht SF'),
    elseif t==4, set(gcf,'Name','ketanserin SF')
    elseif t==5, set(gcf,'Name','ketanserin + doi SF')
    elseif t==6, set(gcf,'Name','mglur2 SF')
    elseif t==7, set(gcf,'Name','mglur2 + DOI SF')
    else set(gcf,'Name','lisuride SF'),end
    useN = find(data & treatment==t)
    for i = 1:length(useN)
        np = ceil(sqrt(length(useN)));
        subplot(np,np,i);
        hold on
        plot(drift_sf(useN(i),:,1,1),'Color',[0.5 0 0]);  plot(drift_sf(useN(i),:,2,1),'Color',[0 0.5 0]); %pre, mv=2
        plot(drift_sf(useN(i),:,1,2),'Color',[1 0 0]);plot(drift_sf(useN(i),:,2,2),'Color',[0 1 0]); %post
        plot([1 7], [1 1]*drift_spont(useN(i),1,1),':','Color',[0.5 0 0]);  plot([1 7], [1 1]*drift_spont(useN(i),2,1),':','Color',[0 0.5 0]);
        plot([1 7], [1 1]*drift_spont(useN(i),1,2),':','Color',[1 0 0]);  plot([1 7], [1 1]*drift_spont(useN(i),2,2),':','Color',[0 1 0]);
        yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)]); xlim([0.5 8.5])
%         if inhAll(useN(i)) ==1 , xlabel('inh'); else  xlabel('exc');
%         end
    end
end

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
        subplot(5,4,s)
        imagesc(squeeze(LFPallCh(use(s),:,:,1,prepost)))
        xlabel('frequency (Hz)');ylabel('tetrode #')
    end
end

use =find(salineSessions==1)
for prepost=1:2
    figure
    if prepost==1, set(gcf,'Name','Saline LFP PRE');
    else set(gcf,'Name','Saline LFP POST'); end
    for s =1:sum(salineSessions)
        subplot(5,4,s)
        imagesc(squeeze(LFPallCh(use(s),:,:,1,prepost)))
        xlabel('frequency (Hz)');ylabel('tetrode #')
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
        xlabel('frequency (Hz)');ylabel('tetrode #')
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
% figure
% for t=1:4
% subplot(2,2,t)
% set(gcf, 'Name', 'evoked LFP')
% plot(squeeze(LFPall(salineSessions==1,:,1,1)),'b');hold on;
% plot(squeeze(LFPall(salineSessions==1,:,1,2)),'r'); xlabel 'Frequency (Hz)'; ylabel 'normalized power';
% end
% %
layerTet = layerSites(:,2:4:64)
use = layerTet(:,:)==4


use =find(salineSessions==1)
for prepost=1:2
    figure
    if prepost==1, set(gcf,'Name','DOI LFP PRE');
    else set(gcf,'Name','DOI LFP POST'); end
    for s =1:sum(salineSessions)
        %subplot(4,4,s)
        imagesc(squeeze(LFPallChDrift(use(s),:,:,1,prepost)))
    end
end

clear s use
use = find(doiSessions==1)
for mv = 1:2
    figure
    if mv==1, set(gcf,'Name', 'moving');
    else set(gcf,'Name', 'stationary')
    end
    for prepost =1:2
        for s =  1:length(use)
            subplot(4,5,s)
         %   plot(squeeze(LFPfreqDrift(:,prepost)),
            plot(squeeze(nanmean(LFPallChDrift(use(s), layerTet(s,:)==4,:,mv,prepost))));
            hold on
        end
    end
end


figure
set(gcf, 'Name', 'all ch darkness LFP')
for t=1:3
subplot(1,3,t)
plot(squeeze(mean(LFPallDark(sessionTreatment==t,:,1,1))),'b');hold on;
plot(squeeze(mean(LFPallDark(sessionTreatment==t,:,1,2))),'r');% xlabel 'Frequency (Hz)'; ylabel 'normalized power';
title(titles{t}); axis square
end

C = {[1 0 0],[.5 0 0]}; %red = mv=1, light =pre, dark =post
D = {[0 1 0],[0 .5 0]}; %green = mv=2
figure
for prepost =1:2
for tet=1:16
%subplot(4,4,tet)
plot(squeeze(LFPallChDark(9,tet,:,1,prepost)),'Color',C{prepost});hold on; xlim([0 100]);
plot(squeeze(LFPallChDark(9,tet,:,2,prepost)),'Color',D{prepost}); xlim([0 100]);
xlabel 'Frequency (Hz)'; ylabel 'normalized power'; %mv
% title(sprintf('layer %d',layerTet(tet)));
end
end

figure
subplot(1,2,1)
imagesc(squeeze(LFPallChDark(1,layerTet==3,:,1,1)));
subplot(1,2,2)
imagesc(squeeze(LFPallChDark(1,layerTet==3,:,1,2))); xlabel 'Frequency (Hz)'; ylabel 'normalized power';




