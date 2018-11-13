clear all
close all
dbstop if error

batchDOIephys_filtered; %%% load batch file
set(groot,'defaultFigureVisible','off') %disable figure plotting
% set(groot,'defaultFigureVisible','on')

%%% select the sessions you want based on filters
%use =  find(strcmp({files.notes},'good data'))%useSess = use;
use =  find( strcmp({files.treatment},'Saline') & strcmp({files.notes},'good data'))
%use =  find( strcmp({files.treatment},'Saline')|strcmp({files.treatment},'DOI') & strcmp({files.notes},'good data'))


%for specific experiment:
%use =  find(strcmp({files.notes},'good data') & strcmp({files.expt},'080715'))
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
    
    %     if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2})
    %         % get pre/post running speed
    %         for prepost = 1:2
    %             spd = getSpeed(clustfile,afile,files(use(i)).blockWn{prepost},0);
    %             speedHistWn(i,:,prepost) = hist(spd.v,0.5:1:100)/length(spd.v);
    %             speedTraceWn{i,prepost}=spd.v;
    %         end
    %     else
    %         speedHistWn(i,:,1:2) = nan;
    %     %    speedTraceWn{i, 1:2} = nan;
    %     end
    
    if ~isempty(files(use(i)).blockDrift{1}) & ~isempty(files(use(i)).blockDrift{2})
        for prepost = 1:2
            spd = getSpeed(clustfile,afile,files(use(i)).blockDrift{prepost},0);
            speedHistDrift(i,:,prepost) = hist(spd.v,0.5:1:100)/length(spd.v);
            speedTraceDrift{i,prepost}=spd.v;
        end
    else
        speedHistDrift(i,:,1:2) = nan;
        %  speedTraceDrift{i, 1:2} = nan;
        
        
    end
    if ~isempty(files(use(i)).blockDark{1}) & ~isempty(files(use(i)).blockDark{2})
        %  hasDarkSess(sessionTreatment) = 1;
        for prepost = 1:2
            spd = getSpeed(clustfile,afile,files(use(i)).blockDark{prepost},0);
            speedHistDark(i,:,prepost) = hist(spd.v,0.5:1:100)/length(spd.v);
            speedTraceDark{i,prepost}=spd.v;
        end
    else
        % hasDarkSess(i) = treatment(cellrange(1))==1;
        speedHistDark(i,:,1:2) = nan;
        % speedTraceDark{i, 1:2} = nan;
    end
    
    figure
    %     for prepost = 1:2
    %         subplot(2,1,prepost);
    %         plot(speedTraceDark{i,prepost}); hold on
    %         plot([1 length(speedTraceDark{i,prepost})] , [ 0.5 0.5], ':'); ylabel('speed'); ylim([0 10])
    %         if prepost==2,  title(sprintf('%s post %s',files(use(i)).expt, files(use(i)).treatment)); end
    %     end
    
    if ~isempty(files(use(i)).blockDrift{1}) & ~isempty(files(use(i)).blockDrift{2})
        %%% get grating responses
        % try
        hasDrift(cellrange)=1;
        for prepost = 1:2
            drift = getDrift_mv(clustfile,afile,files(use(i)).blockDrift{prepost},1);
            drift_spikeT{i,:,prepost} = drift.drift_spikeT;
            drift_orient(cellrange,:,:,prepost)=drift.orient_tune;
            drift_sf(cellrange,:,:,prepost) = drift.sf_tune;
            drift_spont(cellrange,:,prepost) = drift.interSpont;
            drift_fl_spont(cellrange,:,prepost) = drift.spont;
            drift_osi(cellrange,:,prepost) = drift.osi;
            drift_osiCV(cellrange,:,prepost) = drift.cv_osi;
            drift_osiCV2(cellrange,:,prepost) = drift.cv_osi2; %weird values
            %drift_F1F0(:,cellrange,prepost)= drift.F1F0;
            drift_ot_tune(cellrange,:,:,prepost)=drift.orient_tune;
            drift_sf_trial{i,prepost} = drift.trialSF;
            drift_orient_trial{i,prepost} = drift.trialOrient;
            drift_pref_theta(cellrange,:,prepost) = drift.pref_theta;
            drift_pref_thetaCV(cellrange,:,prepost) = drift.pref_thetaCV;
            drift_dsi_theta(cellrange,:,prepost) = drift.dsi_theta;
            drift_dsi(cellrange,:,prepost) = drift.dsi;
            for c = 1:length(cellrange)
                drift_trial_psth{cellrange(c),:,prepost} = squeeze(drift.trialPsth(c,:,:));
                mnPsth(cellrange(c),prepost,:) = squeeze(mean(drift.trialPsth(c,:,:),2));
                drift_spikeTrial{cellrange(c),:,prepost} = drift.drift_spikeT;               
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
            
            
            for c = 1:length(cellrange)
                for cond = 1:72
                    for mv=1:2
                        if mv ==1
                            tr = find(drift.trialOrient(1:end-1) == ceil(cond/6) & drift.trialSF(1:end-1)==mod(cond-1,6)+1 & drift.frameSpd<1);
                        else
                            tr = find(drift.trialOrient(1:end-1) == ceil(cond/6) & drift.trialSF(1:end-1)==mod(cond-1,6)+1 & drift.frameSpd>1);
                        end
                        drift_spikeCounts{cellrange(c),cond,mv,prepost,:} = drift.trialPsth(c,tr,:);
                    end
                end
            end
        end
        
        %             for prepost=1:2
        %                 lfpMoveDrift = getLfpMovement(clustfile,afile,files(use(i)).blockDrift{prepost},0);
        %                 LFPallDrift(i,:,:,prepost) =squeeze(nanmedian(lfpMoveDrift.meanSpect,1));
        %                 if size(lfpMoveDrift.meanSpect,1)==16
        %                     LFPallChDrift(i,:,:,:,prepost) = lfpMoveDrift.meanSpect;
        %                     LFPfreqDrift(:,:) = lfpMoveDrift.freq;
        %                      %LFPtDrift(i,:,prepost) = lfpMoveDrift.t;
        %                     display('good lfp');
        %
        %                 else
        %                     display('lfp wrong size')
        %                 end
        %             end
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
        drift_sf_trial{i,prepost} = NaN;
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
        %
        
        lfpMoveDrift = NaN;
        LFPallDrift(i,300,1:2,prepost) =NaN;
        LFPallChDrift(i,16,300,1:2,prepost) = NaN;
        LFPfreqDrift(1,300) = NaN;
        
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
            %            wn_frameR(cellrange,:,:,prepost) = downsamplebin(wn.frameR,2,15)/15;
            wn_reliability(1,cellrange,prepost) = wn.reliability;
            
        end
    else
        hasWn(cellrange)=0;
        %    wn = NaN;
        wn_crf(cellrange,1:20,1:2,prepost)= NaN;
        wn_spont(cellrange,1:2,prepost)= NaN;
        wn_evoked(cellrange,1:2,prepost)= NaN;
        wn_gain(1,cellrange,prepost)= NaN;
        %   wn_frameR(cellrange,:,1:2,prepost) = NaN;
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
    
    %     if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2})
    %         %% lfp power
    %         %     %%%(right now averages over all sites, should use layer info)
    %         %  try
    %         for prepost=1:2
    %             lfpMove = getLfpMovement(clustfile,afile,files(use(i)).blockWn{prepost},0);
    %             LFPall(i,:,:,prepost) =squeeze(nanmedian(lfpMove.meanSpect,1));
    %             if size(lfpMove.meanSpect,1)==16
    %                 LFPallCh(i,:,:,:,prepost) = lfpMove.meanSpect;
    %                 LFPfreq(:,:) = lfpMove.freq;
    %                 %LFPt(i,:,prepost) = lfpMove.t
    %                 display('good lfp');
    %                 lfpMove
    %             else
    %                 display('lfp wrong size')
    %             end
    %         end
    %         %   catch end
    %     else
    %         lfpMove = NaN;
    %         LFPall(i,300,1:2,prepost) =NaN;
    %         LFPallCh(i,16,300,1:2,prepost) = NaN;
    %         LFPfreq(1,300) = NaN;
    %     end
    
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
                LFPallDark(i,300,1:2,prepost) = NaN;
                
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
        for prepost=1:2
            % cv2Wn(cellrange,:) = cv2;
            for c = 1:length(cellrange)
                spWn{cellrange(c),:,prepost} = R;
            end
        end
    else wnCorr(corrRange,1) = NaN;
        wnCorr(corrRange,2)=NaN;
        % cv2Wn(cellrange,:) = NaN;
        for prepost = 1:2
            for c = 1:length(cellrange)
                spWn{cellrange(c),:,prepost} = NaN;
            end
        end
    end
    
    if ~isempty(files(use(i)).blockDrift{1}) & ~isempty(files(use(i)).blockDrift{2})
        for prepost=1:2
            dt = 1;
            [preCorr postCorr cv2 R eigs] =  prepostDOIdarkness(clustfile,afile,files(use(i)).blockDrift,dt,0);
            driftCorr(corrRange,1) = preCorr(:); driftCorr(corrRange,2)=postCorr(:);
            %cv2Drift(cellrange,:) = cv2;
            for c = 1:length(cellrange)
                spDrift{cellrange(c),:,prepost} = R;
            end
        end
    else
        for prepost = 1:2
            %   driftCorr(corrRange,1) = NaN;
            %  driftCorr(corrRange,2)=NaN;
            % cv2Drift(cellrange,:) = NaN;
            for c = 1:length(cellrange)
                spDrift{cellrange(c),:,prepost}=NaN;
            end
            
        end
        
    end
    
    
    %%%% prepost correlation in darkness
    if ~isempty(files(use(i)).blockDark{1}) & ~isempty(files(use(i)).blockDark{2})
        
        dt = 1;
        [preCorr postCorr cv2 R eigs] = prepostDOIdarkness(clustfile,afile,files(use(i)).blockDark,dt,0);
        darkCorr(corrRange,1) = preCorr(:); darkCorr(corrRange,2)=postCorr(:);
        cv2Dark(cellrange,:) = cv2;
        meanRdark(cellrange,:) = mean(R,2);
        for prepost = 1:2
            for c = 1:length(cellrange)
                spDark{cellrange(c),:,prepost} = R;
            end
        end
    else
        for prepost=1:2
            cv2Dark(cellrange,:) = NaN;
            meanRdark(cellrange,:) = NaN;
            darkCorr(corrRange,1) = NaN;
            darkCorr(corrRange,2)= NaN;
            for c = 1:length(cellrange)
                spDark{cellrange(c),:,prepost} = NaN;
            end
        end
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
%%

%for theoretical people:
 save('compileSalineEphys_101118','layerAll','drift_spikeT','drift_sf_trial','drift_orient_trial', 'drift_trial_psth', 'drift_cond_tcourse', 'speedTraceDrift', 'treatment', 'drift', 'drift_spikeCounts','inhAll')

%%% plot speed histogram
titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
figure
for t =1:2
    subplot(1,2,t)
    spdMn = squeeze(nanmean(speedHistWn(sessionTreatment==t,1:10,:),1));
    spdErr = squeeze(nanstd(speedHistWn(sessionTreatment==t,1:10,:) / sqrt(length(speedHistWn(sessionTreatment==t)))));
    prop = shadedErrorBar(.5:1:10,spdMn(:,1)',spdErr(:,1)','k',1); hold on
    proppost = shadedErrorBar(.5:1:10,spdMn(:,2)',spdErr(:,2)','r',1)
    
    title(titles{t}); xlabel('speed (cm/sec)');axis square; ylim ([0 1]);
end

[h p]=ttest(nanmean(speedHistWn(sessionTreatment==t,:,1),1),nanmean(speedHistWn(sessionTreatment==t,:,2),1))

test1= kstest(speedHistWn(sessionTreatment==t,:,1))

legend('pre','post')
%%
titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
figure
for t =1:2
    subplot(1,2,t)
    spdMn = squeeze(nanmean(speedHistDrift(sessionTreatment==t,1:10,:),1));
    spdErr= squeeze(nanstd(speedHistDrift(sessionTreatment==t,1:10,:)/sqrt(length(speedHistDrift(sessionTreatment==t)))));
    prop = shadedErrorBar(.5:1:10,spdMn(:,1)',spdErr(:,1)','k',1); hold on
    proppost = shadedErrorBar(.5:1:10,spdMn(:,2)',spdErr(:,2)','r',1)
    
    title(titles{t}); xlabel('speed (cm/sec)');axis square; ylim ([0 1]); xlim([0 6])
end
legend('pre','post')



%%
titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
figure
for t =1:2
    subplot(1,2,t)
    spdMn = squeeze(nanmean(speedHistDark(sessionTreatment==t,1:10,:),1));
    spdErr = squeeze(nanstd(speedHistDark(sessionTreatment==t,1:10,:)/sqrt(length(speedHistDark(sessionTreatment==t)))));
    prop = shadedErrorBar(.5:1:10,spdMn(:,1)',spdErr(:,1)','k',1); hold on
    proppost = shadedErrorBar(.5:1:10,spdMn(:,2)',spdErr(:,2)','r',1)
    
    title(titles{t}); xlabel('speed (cm/sec)');axis square; ylim ([0 1]);
end
legend('pre','post')
%% filter data before plotting%%%
clear max_wn
% low_wn = squeeze(min(wn_crf,[],2));
mn_wn_ev = squeeze(nanmean(wn_evoked,2));
mn_wn_spont = squeeze(nanmean(wn_spont,2));
max_wn = (mn_wn_ev-mn_wn_spont);%
% max_wn = (wn_evoked-wn_spont);%
amp_wn = max_wn;  %max_wn-low_wn
useResp_wn = amp_wn(:,1)>2 | amp_wn(:,2)>2; %| amp_wn(:,2,1)>3 | amp_wn(:,2,2)>3;

% useResp_wn = amp_wn(:,1,1)>2 | amp_wn(:,1,2)>2; %| amp_wn(:,2,1)>3 | amp_wn(:,2,2)>3;
data_wn = goodAll==1 & useResp_wn';

clear max amp low
peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont
%peakresp = squeeze(median(drift_orient,2))-squeeze(drift_spont);%subtract driftspont

peakresp = squeeze(mean(drift_orient,2))-squeeze(drift_spont);%subtract driftspont
%
% max_wn = wn_evoked-wn_spont% squeeze(max(wn_evoked,[],2))-(squeeze(wn_spont))
% amp_wn = max_wn  %max_wn-low_wn
% useResp = amp_wn(:,1,1)>3 | amp_wn(:,1,2)>3 %| amp_wn(:,2,1)>3 | amp_wn(:,2,2)>3;
% data_wn = goodAll==1 & useResp';

%%

%evoked FR from drift
titles ={'Saline', 'DOI'};
%find proportions between -1:9
for mv =1
    figure
    if mv==1
        set(gcf,'Name','mean stat drift FR')
    else
        set(gcf,'Name','mean mv drift FR')
    end
    for t=1:2
        subplot(1,2,t)
        used = goodAll==1 & treatment==t & hasDrift==1;
        plot(peakresp(inhAll==0 & used,mv,1),...
            peakresp(inhAll==0 & used,mv,2),'bo','Markersize',10,'linewidth',1)%,'MarkerFaceColor',[.2 0.5 .8]);
        hold on;
        set(gca,'FontSize',25,'linewidth',2);
        xlabel('Pre (spikes/sec)');ylabel('Post (spikes/sec)');
        title(titles{t})
        axis square; xlim([-1 10]);ylim([-1 10]);
        %title(titles{t}); xlabel('Pre spikes/sec');ylabel('Post spikes/sec');
        % ylim ([-5 20]);xlim([-5 20]);
        % ylim ([-8.5 24]);xlim([-8.5 24]);
        % ylim ([-5 5]);xlim([-5 5]);
        mpre=nanmean(peakresp(used,mv,1))
        mpost=nanmean(peakresp(used,mv,2))
       % p(t) = ranksum(peakresp(used,mv,1),(peakresp(used,mv,2)))
       % plot(mpre,mpost,'+k','Markersize',15,'Linewidth',4)
        plot([-10 37], [-10 37])
        plot(peakresp(inhAll==1  &used,mv,1),...
            peakresp(inhAll==1 & used,mv,2),'ro','Markersize',10,'linewidth',1)%,'MarkerFaceColor',[.9 .2 .5]);
        plot(mpre,mpost,'k+','Markersize',20,'Linewidth',4)
        preP = peakresp(used,mv,1); postP = peakresp(used,mv,2);
        propCropped(t) =(sum(preP>10 | preP<-1 |postP>10 | postP<-1)/length(preP)) 
        n_cells(t) = sum(used)
       % text(-5, 32, ['n = ' num2str(n_cells(t))],'FontSize',20)
        for i = 0:1
            clear pre post
            pre = peakresp(used & inhAll==i,mv,1); post = peakresp(used & inhAll==i,mv,2);
            c = corrcoef(pre,post,'rows','pairwise'); c = c(2,1);
            mdl= fitlm(pre,post);

           % pinh(t)= ranksum(pre,post)
           % rsquared(t) = mdl.Rsquared.Ordinary;

             if i ==0, % p(t) = coefTest(mdl);
                [h p(t)] = ttest(pre,post);
                rsquared(t) = mdl.Rsquared.Ordinary;

                 text(-5, 30, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
                text(-5, 28, ['c.c = ' num2str(c,'%.2f')],'FontSize',20)
             dim = [.15 .6 .2 .2];

                 str = {['r^2 = ' num2str(rsquared(t),'%.2f')],...
                    ['c.c = ' num2str(c,'%.2f')],['% not shown = ' num2str(propCropped(t),'%.2f')]...
                    ,['n = ' num2str(n_cells(t))],['p = ' num2str(p(t))]};
                annotation('textbox',dim,'String',str,'FitBoxToText','on')
            else    %pinh(t) = coefTest(mdl)
               [h pinh(t)] = ttest(pre,post)
                   rsquared(t) = mdl.Rsquared.Ordinary;
                text(2, 8, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
                text(2, 6, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);
             text(2,4, ['p Inh = ' num2str(pinh(t))],'FontSize',20);
             end

        end
        axis square;
    end
end

%%
for mv =1
    figure
    if mv==1
        set(gcf,'Name','mean stat drift FR')
    else
        set(gcf,'Name','mean mv drift FR')
    end
    for t=1:2
        usedAll = goodAll==1 & treatment==t & hasDrift==1;
        subplot(1,2,t)
         sn = (unique(sessionNum(treatment==t)))
       for sess= 1:length(unique(sessionNum(treatment==t)))      
           used = goodAll==1 & treatment==t & hasDrift==1 &sessionNum==sn(sess);
           plot(peakresp(inhAll==0 & used,mv,1),...
               peakresp(inhAll==0 & used,mv,2),'bo','Markersize',6,'linewidth',1)%,'MarkerFaceColor',[.2 0.5 .8]);
           hold on;
           premn(sess) = nanmean(peakresp(used,mv,1));  postmn(sess) = nanmean(peakresp(used,mv,2));
           plot(premn(sess),postmn(sess),'r+','Markersize',12,'linewidth',2)
           set(gca,'FontSize',25,'linewidth',2);
           xlabel('Pre (spikes/sec)');ylabel('Post (spikes/sec)');
           title(titles{t})
           axis square; xlim([-1 10]);ylim([-1 10]);
          % animean(sess) = [premn postmn]
       % p(t) = ranksum(peakresp(used,mv,1),(peakresp(used,mv,2)))    
        plot([-10 37], [-10 37])
        plot(peakresp(inhAll==1  &used,mv,1),...
            peakresp(inhAll==1 & used,mv,2),'ro','Markersize',6,'linewidth',1)%,'MarkerFaceColor',[.9 .2 .5]);
       % plot(mpre,mpost,'k+','Markersize',20,'Linewidth',4)
        preP = peakresp(used,mv,1); postP = peakresp(used,mv,2);
        propCropped(t) =(sum(preP>10 | preP<-1 |postP>10 | postP<-1)/length(preP)) 
        n_cells(t) = sum(used)
       end
        mnpre = nanmean(premn);mnpost=nanmean(postmn)
%         mpre=nanmean(peakresp(usedAll,mv,1))
%         mpost=nanmean(peakresp(usedAll,mv,2))  
        plot(mnpre,mnpost,'+k','Markersize',15,'Linewidth',4)
        [h p(t)] = ttest(premn,postmn)
    end
end

%%
% normalized FR
% for mv =1:2
%     figure
%     if mv==1
%         set(gcf,'Name','mean stat drift FR')
%     else
%         set(gcf,'Name','mean mv drift FR')
%     end
%     for t=1:2
%         clear allR minpre maxpre pre post
%         subplot(1,2,t)
%         used = goodAll==1 & treatment==t & hasDrift==1;
%         allR = peakresp(used&inhAll==0,mv,:);
%      if t==1
%         maxpre = max(allR(:,1,1)); minpre = min(allR(:,1,1));
%         allR = allR(:);
%      else t==2
%        % maxpre = max(allR(:,1,1)); minpre = min(allR(:,1,1));
%         [m indPre] = max(allR(:,1,1)); [m2 indPost] = max(allR(:,1,2)); %to remove outlier
%         allRpre=allR([1:indPre-1 indPre+1:end],1,1)
%         allRpost=allR([1:indPost-1 indPost+1:end],1,2)
%         allR= [allRpre;allRpost];
%         maxpre = max(allR(:,1,1)); minpre = min(allR(:,1,1));
% %         allR = allR(:);
%        maxpre = max(allR(1:length(allR)/2)); maxpost = max(allR(length(allR)/2+1:length(allR))) %find new max after removing outlier
%       % allR=allR([1:ind-1 (ind+1:(length(allR)/2)+ind-1) (length(allR)/2)+ind+1:allR(end)])
%      end
%        % normR=(allR-min(allR))/((max(allR)-min(allR)))
%         normR=(allR-min(allR))/(maxpre-minpre)
%         pre = normR(1:length(normR)/2); post = normR(length(pre)+1:length(normR))
%         plot(pre,post, 'ko','Markersize',15,'linewidth',2,'MarkerFaceColor',[.2 0.5 .8]);
%         hold on;
%         clear pre post allR minpre maxpre ind m
%         allR = peakresp(used&inhAll==1,mv,:);
%         maxpre = max(allR(:,1,1)); minpre = min(allR(:,1,1));
%         allR = allR(:);
%         normR=(allR-min(allR))/(max(allR)-min(allR));
%         normR=(allR-min(allR))/(maxpre-minpre);
%         pre = normR(1:length(normR)/2); post = normR(length(pre)+1:length(normR))
%         plot(pre,post, 'ko','Markersize',15,'linewidth',2,'MarkerFaceColor',[.9 .2 .5]);
%
%         set(gca,'FontSize',25,'linewidth',2);
%        % xlabel('Pre (spikes/sec)');ylabel('Post (spikes/sec)');
%         title(titles{t})
%         axis square; xlim([0 1]);ylim([0 1]);
%         title(titles{t}); xlabel('pre normalized firing rate');ylabel('post normalized firing rate');
%
% %         mpre=nanmean(peakresp(used,mv,1))
% %         mpost=nanmean(peakresp(used,mv,2))
% %
% %         plot(mpre,mpost,'+k','Markersize',15,'Linewidth',4)
%          plot([0 1], [0 1])
%       %  plot(mpre,mpost,'+k','Markersize',20,'Linewidth',4)
%
%         n_cells(t) = sum(used)
%         text(-5, 32, ['n = ' num2str(n_cells(t))],'FontSize',20)
%         for i = 0:1
%             pre = peakresp(used & inhAll==i,mv,1); post = peakresp(used & inhAll==i,mv,2);
%             c = corrcoef(pre,post,'rows','pairwise'); c = c(2,1);
%             mdl= fitlm(pre,post);
%             rsquared(t) = mdl.Rsquared.Ordinary;
%             if i ==0, text(-5, 30, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
%                 text(-5, 28, ['c.c = ' num2str(c,'%.2f')],'FontSize',20)
%             else text(-5, 26, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
%                 text(-5, 24, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
%         end
%         axis square;
%     end
% end
%%
peakresp_mn = squeeze(mean(peakresp,2))
titles ={'Saline', 'DOI','5HT'};
for mv =2
    figure();
    if mv==1
        set(gcf,'Name','mean stat drift FR')
    else
        set(gcf,'Name','mean mv drift FR')
    end
    for t=1:2
        subplot(1,2,t)
        used = goodAll==1 & treatment==t & hasDrift==1;
        plot(peakresp_mn(inhAll==0 & used,1), peakresp_mn(inhAll==0 & used,2),'.','Markersize',20);
        hold on;
        set(gca,'FontSize',18);
        xlabel('Pre (spikes/sec)');ylabel('Post (spikes/sec)');
        title(titles{t})
        axis square; %xlim([0 20]);ylim([0 20]);
        %title(titles{t}); xlabel('Pre spikes/sec');ylabel('Post spikes/sec');
        ylim ([-2.5 33]);xlim([-2.5 33]);
        mpre=nanmean(peakresp_mn(inhAll==1 & used,1)); mpost=nanmean(peakresp_mn(inhAll==1 & used,2))
        plot(mpre,mpost,'+k','Markersize',24,'Linewidth',4)
        plot([-10 35], [-10 35])
        plot(peakresp_mn(used & inhAll==1,1), peakresp_mn(used & inhAll==1,2),'r.','Markersize',20);
        n_cells(t) = sum(goodAll==1 &treatment==t& hasDrift==1)
        text(-5, 19, ['n = ' num2str(n_cells(t))],'FontSize',20)
        for i = 0:1
            pre = peakresp_mn(used & inhAll==i,1); post = peakresp_mn(used & inhAll==i,2);
            c = corrcoef(pre,post,'rows','pairwise'); c = c(2,1);
            mdl= fitlm(pre,post)
            rsquared(t) = mdl.Rsquared.Ordinary
            if i ==0, text(-5, 17, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
                text(-5, 15, ['c.c = ' num2str(c,'%.2f')],'FontSize',20)
            else text(-5, 13, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
                text(-5, 11, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
        end
    end
end

%%
%Drift-evoked Proportion of cells suppressed or facilitated pre/post
clear fracSupp stderrorSupp fracInc stderrorInc deltaR preR postR
clear evoked Proportion of cells suppressed or facilitated pre/post
figure
for mv = 1:2
    for t=1:2
        preR = mean(drift_orient(goodAll==1&treatment==t& hasDrift==1,:,mv,1),2);
        postR = mean(drift_orient(goodAll==1  &treatment==t& hasDrift==1,:,mv,2),2);
        deltaR = postR-preR;
        stderror = nanstd(deltaR) / sqrt(length(deltaR));
        fracSupp(t) = sum(deltaR<-1)/length(deltaR);
        stderrorSupp(t) = nanstd(fracSupp) / sqrt(length(fracSupp));
        fracInc(t) = sum(deltaR>1)/length(deltaR);
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
    axis square; ylabel('Proportion of cells','FontSize',18); ylim([0 0.25]); xlim([.5 2.5]);
    Labels = {'Saline','DOI'};set(gca, 'XTick', 1:2, 'XTickLabel', Labels,'FontSize',18);
    
end

%% Modulation indices for drift, stationary
useResp = peakresp(:,1,1)>2| peakresp(:,1,2)>2
layerz={'l1','l2','l3','layer 4','layer 5','L6'};
clear p23 pInh p4 p5
for t=1:2
    figure
    if t == 1
        set(gcf,'Name', 'hist saline');
    else
        set(gcf,'Name', 'mi hist doi');
    end
    useN = goodAll==1 & treatment==t & hasDrift==1& inhAll==0 &useResp';
    useInh = goodAll==1 & treatment==t & hasDrift==1& inhAll==1 &useResp';
    pre = peakresp(useN & (layerAll==2|layerAll==3),1,1); post = peakresp(useN & (layerAll==2|layerAll==3),1,2);
    miDrift23 = (post-pre)./(post+pre);
  %  miDrift23(find(miDrift23<-1))=-1; miDrift23(find(miDrift23>1))=1
  %  m23 = mean(abs(miDrift23));
    m23(t) = mean(abs(miDrift23));
    s23(t) = std(miDrift23);
    h23 = hist(miDrift23,-1:.2:1);
    t23(:,t) = hist(miDrift23,-1:.2:1);
    Mbins = -1:.2:1;
    subplot(1,4,1)
    bar(Mbins,h23/length(pre),'FaceColor',[0 .5 .5],'Linewidth',.5);
    hold on; plot(m23,.25,'p'); plot([0 0],[0 1],'--');
    ylim([0 .5]);xlim([-1.5 1.5]); axis square
    title('layer 2/3')
    clear miDrift pre post
    for i=4:5
        pre = peakresp(useN & layerAll==i, 1,1); post = peakresp(useN & layerAll==i,1,2);
        miDrift = (post-pre)./(post+pre)
       % miDrift(find(miDrift<-1))=-1;miDrift(find(miDrift>1))=1
      %  m45(i) = mean(abs(miDrift));
        m45(i,t) = mean(abs(miDrift));
        s45(i,t)=std(miDrift);
        h = hist(miDrift,-1:.2:1);
        t45(t,:,i) = hist(miDrift,-1:.2:1); 
        subplot(1,4,i-2)
        Mbins = -1:.2:1;
        bar(Mbins,h/sum(useN&layerAll==i),'FaceColor',[0 .5 .5],'Linewidth',1);
        ylim([0 .5]);xlim([-1.5 1.5]); axis square
        title(layerz{i});
        hold on;plot(m45(i),.25,'p');plot([0 0],[0 1],'--');
    end
    clear miDrift 
    subplot(1,4,4)
    preInh= peakresp(useInh,1,1);postInh = peakresp(useInh,1,2);
    miDriftInh = (postInh-preInh)./(postInh+preInh);
   % miDriftInh(find(miDriftInh<-1))=-1; miDriftInh(find(miDriftInh>1))=1
   % mnInh = mean(abs(miDriftInh));
    mnInh(t) = mean(abs(miDriftInh));
    sInh(t) = std(miDriftInh);
    hInh = hist(miDriftInh,-1:.2:1);
    tInh(:,t) = hist(miDriftInh,-1:.2:1);
    bar(Mbins,hInh/sum(useInh),'FaceColor',[0 .5 .5],'Linewidth',1); 
    hold on; plot(mnInh,.25,'p', 'linewidth',1);plot([0 0],[0 1],'--');   
    ylim([0 .5]); xlim([-1.5 1.5]); axis square
    title('inhibitory');
   end
    
% [k pInh] = kstest2(tInh(:,1),tInh(:,2));
% [k2 p23] = kstest2(t23(:,1),t23(:,2));
% [k4 p4] = kstest2(t45(1,:,4),t45(2,:,4));
% [k5 p5] = kstest2(t45(1,:,5),t45(2,:,5));

[k p]=ttest(miDrift23)

disp(m23);disp(m45);disp(mnInh);
table (m23',s23', m45(1,4:5)', s45(4:5,:),mnInh',sInh')

%% Modulation indices for drift, stationary
layerz={'l1','l2','l3','layer 4','layer 5','L6'};
for t=1:2
    figure
    if t == 1
        set(gcf,'Name', 'hist saline');
    elseif t==2
        set(gcf,'Name', 'mi hist doi');
    else
        set(gcf,'Name', 'mi hist 5ht');
    end
    useN = goodAll==1 & treatment==t & hasDrift==1& inhAll==0;
    useInh = goodAll==1 & treatment==t & hasDrift==1& inhAll==1;
    
    miDrift = (mean(peakresp(useN &(layerAll==2|layerAll==3),1,2),3)-mean(peakresp(useN&(layerAll==2|layerAll==3),1,1),3))./...
        (mean(peakresp(useN&(layerAll==2|layerAll==3),1,2),3)+mean(peakresp(useN&(layerAll==2|layerAll==3),1,1),3));
    h = hist(miDrift,-1:.2:1);
    Mbins = -1:.2:1;
    subplot(1,4,1)
    bar(Mbins,h/sum(useN &(layerAll==2|layerAll==3)),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .5]);xlim([-1.5 1.5]); axis square
    title('layer 2/3')
    clear miDrift
    for i=4:5
        miDrift= (mean(peakresp(useN &(layerAll==i),1,2),3)-mean(peakresp(useN&(layerAll==i),1,1),3))./...
            (mean(peakresp(useN&(layerAll==i),1,2),3)+mean(peakresp(useN&(layerAll==i),1,1),3));
        h = hist(miDrift,-1:.2:1);
        subplot(1,4,i-2)
        Mbins = -1:.2:1;
        bar(Mbins,h/sum(useN&layerAll==i),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .5]);xlim([-1.5 1.5]); axis square
        title(layerz{i});
    end
    subplot(1,4,4)
    useN = goodAll==1 & treatment==t & hasDrift==1& inhAll==1;
    miDrift= (mean(peakresp(useInh,1,2),3)-mean(peakresp(useInh,1,1),3))./...
        (mean(peakresp(useInh,1,2),3)+mean(peakresp(useInh,1,1),3));
    h = hist(miDrift,-1:.2:1);
    bar(Mbins,h/sum(useInh),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .5]);xlim([-1.5 1.5]); axis square
    title('inhibitory');
    
end

%%

%%MI for evoked drift, move%%
layerz={'l1','l2','l3','layer 4','layer 5'};
for t=1:2
    figure
    if t == 1
        set(gcf,'Name', 'hist saline');
    elseif t==2
        set(gcf,'Name', 'mi hist doi');
    else
        set(gcf,'Name', 'mi hist 5ht');
    end
    useN=goodAll & treatment==t & ~inhAll & hasDrift==1;
    miDrift= (mean(peakresp(useN &(layerAll==2|layerAll==3),2,2),3)-mean(peakresp(useN&(layerAll==2|layerAll==3),2,1),3))./...
        (mean(peakresp(useN&(layerAll==2|layerAll==3),2,2),3)+mean(peakresp(useN&(layerAll==2|layerAll==3),2,1),3));
    h= hist(miDrift,-1:.2:1);
    Mbins=-1:.2:1;
    subplot(1,4,1)
    bar(Mbins,h/sum(useN &(layerAll==2|layerAll==3)),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .5]);xlim([-1.5 1.5]); axis square
    title('layer 2/3')
    for i=4:5
        miDrift= (mean(peakresp(useN &(layerAll==i),2,2),3)-mean(peakresp(useN&(layerAll==i),2,1),3))./...
            (mean(peakresp(useN&(layerAll==i),2,2),3)+mean(peakresp(useN&(layerAll==i),2,1),3));
        h= hist(miDrift,-1:.2:1);
        subplot(1,4,i-2)
        Mbins=-1:.2:1;
        bar(Mbins,h/sum(useN&layerAll==i),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .5]);xlim([-1.5 1.5]); axis square
        title(layerz{i});
        
    end
    subplot(1,4,4)
    useN = goodAll==1 & treatment==t & hasDrift==1& inhAll==1;
    miDrift= (mean(peakresp(useInh,2,2),3)-mean(peakresp(useInh,2,1),3))./...
        (mean(peakresp(useInh,2,2),3)+mean(peakresp(useInh,2,1),3));
    h = hist(miDrift,-1:.2:1);
    bar(Mbins,h/sum(useInh),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .5]);xlim([-1.5 1.5]); axis square
    title('inhibitory');
end

%%
%MI for evoked drift, avg mv and stat%%
mnpeakresp = squeeze(mean(peakresp,2));

layerz={'l1','l2','l3','layer 4','layer 5'};
for t=1:2
    figure
    if t == 1
        set(gcf,'Name', 'hist saline');
    elseif t==2
        set(gcf,'Name', 'mi hist doi');
    else
        set(gcf,'Name', 'mi hist 5ht');
    end
    useN=goodAll & treatment==t & ~inhAll & hasDrift==1;
    miDrift= (mean(mnpeakresp(useN &(layerAll==2|layerAll==3),2),2)-mean(peakresp(useN&(layerAll==2|layerAll==3),1),2))./...
        (mean(peakresp(useN&(layerAll==2|layerAll==3),2),2)+mean(peakresp(useN&(layerAll==2|layerAll==3),1),2));
    h= hist(miDrift,-1:.2:1);
    Mbins=-1:.2:1;
    subplot(1,4,1)
    bar(Mbins,h/sum(useN &(layerAll==2|layerAll==3)),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .5]);xlim([-1.5 1.5]); axis square
    title('layer 2/3')
    for i=4:5
        miDrift= (mean(mnpeakresp(useN &(layerAll==i),2),2)-mean(mnpeakresp(useN&(layerAll==i),1),2))./...
            (mean(mnpeakresp(useN&(layerAll==i),2),2)+mean(mnpeakresp(useN&(layerAll==i),1),2));
        h= hist(miDrift,-1:.2:1);
        subplot(1,4,i-2)
        Mbins=-1:.2:1;
        bar(Mbins,h/sum(useN&layerAll==i),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .4]);xlim([-1.5 1.5]); axis square
        title(layerz{i});
    end
    subplot(1,4,4)
    useInh=goodAll & treatment==t & inhAll==1 & hasDrift==1;
    miDrift= (mean(mnpeakresp(useInh,2),2)-mean(mnpeakresp(useInh,1),2))./...
        (mean(mnpeakresp(useInh,2),2)+mean(mnpeakresp(useInh,1),2));
    h= hist(miDrift,-1:.2:1);
    bar(Mbins,h/sum(useN&layerAll==i),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .4]);xlim([-1.5 1.5]); axis square
    title ('inhibitory');
    
end

%%
%spont FR from drift
clear mdl rsquared c

titles ={'Saline','DOI','5HT'};
for mv =1:2
    figure
    if mv==1
        set(gcf,'Name','mean stationary')
    else
        set(gcf,'Name','mv')
    end
    for t=1:2
        subplot(1,2,t)
        used = goodAll==1 &treatment==t & hasDrift==1;
        plot(drift_spont(used & inhAll==0,mv,1),drift_spont(used &inhAll==0,mv,2),'ko','Markersize',15,'linewidth',2,'MarkerFaceColor',[.2 0.5 .8]);  hold on;
        set(gca,'FontSize',18);
        xlabel('Pre (spikes/sec)'); ylabel('Post (spikes/sec)');title(titles{t})
        ylim([0 16]); xlim([0 16]);
        set(gca,'FontSize',25,'linewidth',2);
        axis square
        plot([0 35], [0 35])
        plot(drift_spont(used &inhAll==1,mv,1),drift_spont(used &inhAll==1,mv,2),'ko','Markersize',15,'linewidth',2,'MarkerFaceColor',[.9 .2 .5]);
        n_cells(t) = sum(used)
        plot(nanmean(drift_spont(used&inhAll==1, mv,1)),nanmean(drift_spont(used&inhAll==1, mv, 2)),'+y','Markersize',20,'Linewidth',6)
        text(1, 15, ['n = ' num2str(n_cells(t))],'FontSize',18)
        for i = 0:1
            pre = drift_spont(used & inhAll==i,mv,1); post = drift_spont(used & inhAll==i,mv,2);
            c = corrcoef(pre,post,'rows','pairwise'); c = c(2,1);
            mdl= fitlm(pre,post)
            rsquared(t) = mdl.Rsquared.Ordinary
            if i ==0, text(1, 14, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
                text(1, 13, ['c.c = ' num2str(c,'%.2f')],'FontSize',20)
            else text(1, 12, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
                text(1, 11, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
        end
        axis square
    end
end

%%

%MI for drift spont mv
layerz={'l1','l2','l3','layer 4','layer 5'};
for t=1:2
    figure
    if t == 1
        set(gcf,'Name', 'hist saline');
    elseif t==2
        set(gcf,'Name', 'mi hist doi');
    else
        set(gcf,'Name', 'mi hist 5ht');
    end
    useN=goodAll & treatment==t
    miDrift= (mean(drift_spont(useN &(layerAll==2|layerAll==3),2,2),2)-...
        mean(drift_spont(useN&(layerAll==2|layerAll==3),2,1),2))./...
        (mean(drift_spont(useN&(layerAll==2|layerAll==3),2,2),2)+mean(drift_spont(useN&(layerAll==2|layerAll==3),2,1),2));
    h = hist(miDrift,-1:.2:1);
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


%% drift spont mean of mv and stat

clear mdl_dSpont
%spont FR from drift
drift_spont_all = squeeze(mean(drift_spont,2));
titles ={'Saline', 'DOI','5HT'};
figure
for t=1:2
    subplot(1,2,t)
    used = goodAll==1 &treatment==t & hasDrift==1;
    plot(drift_spont_all(used & inhAll==0,1),drift_spont_all(used& inhAll==0,2),'.','Markersize',20);  hold on;
    set(gca,'FontSize',18);
    xlabel('Pre (spikes/sec)'); ylabel('Post (spikes/sec)');title(titles{t})
    % ylim([0 35]); xlim([0 35]);
    axis square
    plot([0 33], [0 33])
    plot(drift_spont_all(used & inhAll==1,1),drift_spont_all(used & inhAll==1,2),'r.','Markersize',20);
    n_cells(t) = sum(used)
    text(2, 32, ['n = ' num2str(n_cells(t))],'FontSize',18)
    for i = 0:1
        pre = drift_spont_all(used & inhAll==i,1); post = drift_spont_all(used & inhAll==i,2);
        c = corrcoef(pre,post,'rows','pairwise');c = c(2,1);
        mdl= fitlm(pre,post)
        rsquared(t) = mdl.Rsquared.Ordinary
        if i ==0, text(2, 30, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
            text(2, 28, ['c.c = ' num2str(c,'%.2f')],'FontSize',20)
        else text(2, 26, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
            text(2, 24, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
    end
end


%%
clear cycR
for f = 1:20;
    cycR(:,f,:,:) = nanmean(wn_frameR(:,f:20:end,:,:),2);
end
spont = squeeze(mean(cycR(:,[1 2 19 20],:,:),2));
evoked = squeeze(mean(cycR(:,9:12,:,:),2)) - spont;

clear max_wn
% max_wn = wn_evoked-wn_spont; % squeeze(max(wn_evoked,[],2))-(squeeze(wn_spont))
% amp_wn = max_wn;
amp_wn = wn_evoked;
useResp = amp_wn(:,1,1)>2 | amp_wn(:,1,2)>2 %| amp_wn(:,2,1)>3 | amp_wn(:,2,2)>3;
data_wn = goodAll==1 & hasWn==1 %& useResp';

clear max_wn
max_wn = wn_evoked-wn_spont% squeeze(max(wn_evoked,[],2))-(squeeze(wn_spont))
amp_wn = max_wn  %max_wn-low_wn
useResp = amp_wn(:,1,1)>3 | amp_wn(:,1,2)>3 %| amp_wn(:,2,1)>3 | amp_wn(:,2,2)>3;
data_wn = goodAll==1 %& useResp';

%%
%wn evoked mv and stationary
titles = {'Saline','DOI'}
for mv=1:2
    figure
    if mv==1, set(gcf,'Name', 'wn evoked stationary'); else set(gcf,'Name','wn evoked mv');end
    for t=1:2
        used = goodAll==1 & treatment==t & inhAll==0 & hasWn==1;
        usedInh = goodAll==1 & treatment==t & inhAll==1 & hasWn==1;
        usedAll = goodAll==1&treatment==t& hasWn==1;
        subplot(1,2,t)
        plot(amp_wn(used,mv,1),amp_wn(used,mv,2),'.','Markersize',20);
        hold on;axis square;
        %xlim([-1 25]); ylim([-1 25]);
        plot([-30 30],[-30 30]);title(titles{t})
        plot(amp_wn(usedInh,mv,1),amp_wn(usedInh,mv,2),'.r','Markersize',20);
        n_cells(t) = sum(usedAll)
        text(-20, 33, ['n = ' num2str(n_cells(t))],'FontSize',20)
        
        for i = 0:1
            pre = amp_wn(usedAll& inhAll==i,mv,1); post = amp_wn(usedAll & inhAll==i,mv,2);
            c = corrcoef(pre,post,'rows','pairwise');c = c(2,1);
            mdl= fitlm(pre,post);
            rsquared(t) = mdl.Rsquared.Ordinary;
            if i ==0, text(-25, 29, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20);
                text(-25, 25, ['c.c = ' num2str(c,'%.2f')],'FontSize',20);
            else text(-25, 21, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20);
                text(-25, 17, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
        end
    end
end

%%
%%

%spont FR from wn
clear mdl rsquared c
titles ={'Saline', 'DOI','5HT'};
for mv=1:2
    figure
    if mv==1, set(gcf,'Name', 'wn spont stationary');
    else set(gcf,'Name', 'wn spont mv');end
    for t=1:2
        subplot(1,2,t)
        used = goodAll==1 & treatment==t & hasWn==1;
        plot(wn_spont(used & inhAll==0,mv,1), wn_spont(used & inhAll==0,mv,2),'.', 'Markersize', 20);
        hold on;axis square;plot([-30 37],[-30 37]);
        plot(wn_spont(used&inhAll==1,mv,1),wn_spont(used&inhAll==1,mv,2),'.r', 'Markersize', 20);
        xlim([0 37]); ylim([0 37]);
        title(titles{t});
        n_cells(t) = sum(used);
        text(2, 33, ['n = ' num2str(n_cells(t))],'FontSize',20);
        for i = 0:1
            pre = wn_spont(used & inhAll==i,mv,1); post = wn_spont(used & inhAll==i,mv,2);
            c = corrcoef(pre,post,'rows','pairwise');c = c(2,1);
            mdl= fitlm(pre,post)
            rsquared(t) = mdl.Rsquared.Ordinary
            if i ==0, text(2, 29, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
                text(2, 25, ['c.c = ' num2str(c,'%.2f')],'FontSize',20)
            else text(2, 21, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
                text(2, 17, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
        end
    end
end

%% wn_spont, mean of stationary & running
mn_wn_spont = squeeze(mean(wn_spont,2));
titles = {'Saline','DOI','5HT'};
figure
for t=1:2
    subplot(1,2,t)
    used = goodAll==1 & treatment==t & hasWn==1;
    plot(mn_wn_spont(used & inhAll==0,1), mn_wn_spont(used & inhAll==0,2),'.', 'Markersize',20);
    hold on;axis square;plot([-30 30],[-30 30]);
    %     xlim([min(wn_spont(data_wn,mv,1))-.5 max(wn_spont(data_wn,mv,1))+.5]); ylim([min(wn_spont(data_wn,mv,1))-.5 max(wn_spont(data_wn,mv,1))+.5]);
    plot(mn_wn_spont(used&inhAll==1,1),mn_wn_spont(used&inhAll==1,2),'.r', 'Markersize',20);
    xlim([0 42]); ylim([0 42]);
    title(titles{t});
    n_cells(t) = sum(used);
    text(2, 40, ['n = ' num2str(n_cells(t))],'FontSize',20);
    for i = 0:1
        pre = mn_wn_spont(used & inhAll==i,1); post = mn_wn_spont(used & inhAll==i,2);
        c = corrcoef(pre,post,'rows','pairwise');c = c(2,1);
        mdl= fitlm(pre,post)
        rsquared(t) = mdl.Rsquared.Ordinary
        if i ==0, text(2, 38, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
            text(2, 36, ['c.c = ' num2str(c,'%.2f')],'FontSize',20)
        else text(2, 34, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
            text(2, 32, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
    end
end


%%
%evoked calculcated from cyc avg
clear n_cells
for mv = 1:2
    figure
    if mv==1, set(gcf,'Name', 'evoked pre vs post stationary')
    else set(gcf,'Name', 'evoked pre vs post mv'), end
    for t=1:2
        subplot(1,2,t)
        used = goodAll==1 & treatment==t & hasWn==1;
        plot(evoked(used & inhAll==0,mv,1),evoked(used & inhAll==0,mv,2),'.','MarkerSize',20); hold on
        plot(evoked(used & inhAll==1,mv,1),evoked(used & inhAll==1,mv,2),'.r','MarkerSize',20);
        axis square; plot([-30 30],[-30 30]);
        xlim([-24 24]);ylim([-24 24])
        n_cells(t) = sum(used);
        text(-25, 26, ['n = ' num2str(n_cells(t))],'FontSize',20); set(gca,'FontSize',20);
        for i = 0:1
            pre = evoked(used & inhAll==i,mv,1); post = evoked(used & inhAll==i,mv,2);
            c = corrcoef(pre,post,'rows','pairwise'); c = c(2,1);
            mdl= fitlm(pre,post)
            rsquared(t) = mdl.Rsquared.Ordinary
            if i ==0, text(-25, 22, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
                text(-25, 18, ['c.c = ' num2str(c,'%.2f')],'FontSize',20)
            else text(-25, 14, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
                text(-25, 10, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
        end
    end
end

%%
%wn evoked from cyc avg, avg of stat and running
evoked_mn = squeeze(mean(evoked,2));
clear n_cells
figure
set(gcf,'Name', 'evoked mean of stationary and running')
for t=1:2
    subplot(1,2,t)
    used = goodAll==1 & treatment==t & hasWn==1;
    plot(evoked_mn(used & inhAll==0,1),evoked_mn(used & inhAll==0,2),'.','MarkerSize',20);
    title(titles{t});
    hold on; plot([-50 50],[-50 50]); axis xy
    %xlabel('pre'); ylabel('post');% xlim([min(evoked(data_wn,mv,1))-.5 max(evoked(data_wn,mv,1))+.5]); ylim([min(evoked(data_wn,mv,1))-.5 max(evoked(data_wn,mv,1))+.5]);
    plot(evoked_mn(used & inhAll==1,1),evoked_mn(used & inhAll==1,2),'.r','MarkerSize',20); axis square;
    xlim([-24 24]);ylim([-24 24])
    n_cells(t) = sum(used);
    text(-25, 26, ['n = ' num2str(n_cells(t))],'FontSize',20); set(gca,'FontSize',20);
    for i = 0:1
        pre = evoked_mn(used & inhAll==i,1); post = evoked_mn(used & inhAll==i,2);
        c = corrcoef(pre,post,'rows','pairwise'); c = c(2,1);
        mdl= fitlm(pre,post)
        rsquared(t) = mdl.Rsquared.Ordinary
        if i ==0, text(-25, 22, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
            text(-25, 18, ['c.c = ' num2str(c,'%.2f')],'FontSize',20)
        else text(-25, 14, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
            text(-25, 10, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
    end
end

%%
%spont calculcated from cyc avg
clear used
titles = {'Saline','DOI','5HT'};
for mv = 1:2
    figure
    if mv==1, set(gcf,'Name', 'wn spont stationary')
    else set(gcf,'Name', 'wn spont mv'), end
    for t=1:2
        subplot(1,2,t)
        used = goodAll==1 & treatment==t & hasWn==1;
        plot(spont(used & inhAll==0,mv,1),spont(used & inhAll==0,mv,2),'.','Markersize',20);
        title(titles{t}); hold on; axis square; plot([0 50],[0 50]);
        xlabel('pre'); ylabel('post');
        plot(spont(used & inhAll==1,mv,1),spont(used & inhAll==1,mv,2),'.r','Markersize',20);
        hold on; xlim([0 16]); ylim([0 16]);
        %     axis([0 36 0 36]);
        n_cells(t) = sum(used)
        text(1, 19, ['n = ' num2str(n_cells(t))],'FontSize',20)
        for i = 0:1
            pre = spont(used & inhAll==i,mv,1); post = spont(used & inhAll==i,mv,2);
            c = corrcoef(pre,post,'rows','pairwise');c = c(2,1);
            mdl= fitlm(pre,post)
            rsquared(t) = mdl.Rsquared.Ordinary
            if i ==0, %text(2, 33, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
                text(1, 17, ['c.c = ' num2str(c,'%.2f')],'FontSize',20)
            else %text(2, 16, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
                text(1, 15, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
        end
    end
end

%%

data_wn = goodAll==1 & hasWn==1 & inhAll==0;
dataInh = goodAll==1& hasWn==1 & inhAll==1;

titles = {'saline','doi','ht','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
figure
for t=1:2
    set(gcf,'Name', 'spont pre vs post stop')
    subplot(1,2,t)
    plot(spont(find(data_wn & treatment==t),1,1),spont(find(data_wn& treatment==t),1,2),'.','Markersize',20); title(titles{t}); hold on; axis square; plot([-50 50],[-50 50]);
    xlabel('pre'); ylabel('post');
    plot(spont(find(dataInh&inhAll& treatment==t),1,1),spont(find(dataInh & inhAll& treatment==t),1,2),'.r','Markersize',20); hold on; plot([-50 50],[-50 50]);
    axis([0 10 0 10])
end

figure
for t=1:2
    set(gcf,'Name', 'spont pre vs post mv')
    subplot(1,2,t)
    plot(spont(find(data_wn& treatment==t),2,1),spont(find(data_wn& treatment==t),2,2),'.','Markersize',20); title(titles{t}); hold on; plot([-50 50],[-50 50]); axis square
    xlabel('pre'); ylabel('post');
    plot(spont(find(dataInh&inhAll& treatment==t),2,1),spont(find(dataInh & inhAll& treatment==t),2,2),'.r','Markersize',20); hold on; plot([-50 50],[-50 50]); axis([0 15 0 15])
end

%%
figure
for t=1:2
    set(gcf,'Name', 'evoked pre vs post stop')
    subplot(1,2,t)
    plot(evoked(find(data_wn& treatment==t),1,1),evoked(find(data_wn& treatment==t),1,2),'.','Markersize',20); title(titles{t}); hold on; plot([-50 50],[-50 50]); axis square
    xlabel('pre'); ylabel('post');
    plot(evoked(find(dataInh&inhAll& treatment==t),1,1),evoked(find(dataInh & inhAll& treatment==t),1,2),'.r','Markersize',20); hold on; plot([-50 50],[-50 50]); axis([-10 10 -10 10])
end
figure
for t=1:2
    set(gcf,'Name', 'evoked pre vs post mv')
    subplot(1,2,t)
    plot(evoked(find(data_wn& treatment==t),2,1),evoked(find(data_wn& treatment==t),2,2),'.','Markersize',20); title(titles{t}); hold on; plot([-50 50],[-50 50]); axis square
    xlabel('pre'); ylabel('post');
    plot(evoked(find(dataInh&inhAll& treatment==t),2,1),evoked(find(dataInh & inhAll& treatment==t),2,2),'.r','Markersize',20); hold on; plot([-50 50],[-50 50]); axis([-10 10 -10 10])
end

%%
mnspont=squeeze(mean(spont,2))
titles = {'Saline','DOI','5HT'};
figure
for t=1:2
    clear used
    subplot(1,2,t)
    used = goodAll==1 & treatment==t & hasWn==1;
    plot(mnspont(used & inhAll==0,1),mnspont(used & inhAll==0,2),'.','Markersize',20);
    title(titles{t}); hold on; axis square; plot([0 50],[0 50]);
    xlabel('pre'); ylabel('post');
    plot(mnspont(used & inhAll==1,1),mnspont(used & inhAll==1,2),'.r','Markersize',20);
    hold on; plot([0 50],[0 50]);
    axis([0 40.5 0 40.5]);
    n_cells(t) = sum(used)
    text(2, 39, ['n = ' num2str(n_cells(t))],'FontSize',20)
    for i = 0:1
        pre = spont(used & inhAll==i,1); post = spont(used & inhAll==i,2);
        c = corrcoef(pre,post,'rows','pairwise');c = c(2,1);
        mdl= fitlm(pre,post)
        rsquared(t) = mdl.Rsquared.Ordinary
        if i ==0, text(2, 37, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
            text(2, 35, ['c.c = ' num2str(c,'%.2f')],'FontSize',20)
        else text(2, 33, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
            text(2, 31, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
    end
end

%%
layerz={'l1','l2','l3','layer 4','layer 5','L6'};
for t=1:2
    figure
    if t == 1
        set(gcf,'Name', 'hist saline');
    elseif t==2
        set(gcf,'Name', 'mi hist doi');
    else
        set(gcf,'Name', 'mi hist 5ht');
    end
    useN = goodAll==1 & treatment==t & hasWn==1& inhAll==0;
    miDrift = (mean(spont(useN &(layerAll==2|layerAll==3),1,2),3)-mean(spont(useN&(layerAll==2|layerAll==3),1,1),3))./...
        (mean(spont(useN&(layerAll==2|layerAll==3),1,2),3)+mean(spont(useN&(layerAll==2|layerAll==3),1,1),3));
    h = hist(miDrift,-1:.2:1);
    Mbins = -1:.2:1;
    subplot(1,4,1)
    bar(Mbins,h/sum(useN &(layerAll==2|layerAll==3)),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .3]);xlim([-1.5 1.5]); axis square
    title('layer 2/3')
    clear miDrift
    for i=4:5
        miDrift= (mean(spont(useN &(layerAll==i),1,2),3)-mean(spont(useN&(layerAll==i),1,1),3))./...
            (mean(spont(useN&(layerAll==i),1,2),3)+mean(spont(useN&(layerAll==i),1,1),3));
        h = hist(miDrift,-1:.2:1);
        subplot(1,4,i-2)
        Mbins = -1:.2:1;
        bar(Mbins,h/sum(useN&layerAll==i),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .3]);xlim([-1.5 1.5]); axis square
        title(layerz{i});
    end
    subplot(1,4,4)
    useInh = goodAll==1 & treatment==t & hasWn==1& inhAll==1;
    miDrift= (mean(spont(useInh,1,2),3)-mean(spont(useInh,1,1),3))./...
        (mean(spont(useInh,1,2),3)+mean(spont(useInh,1,1),3));
    h = hist(miDrift,-1:.2:1);
    bar(Mbins,h/sum(useInh),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .3]);xlim([-1.5 1.5]); axis square
    title('inhibitory');
    
end


%% mean of wn and drift
peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont
%peakresp = squeeze(mean(drift_orient,2))-squeeze(drift_spont);
amp_wn = wn_evoked-wn_spont;
allEv1 = [peakresp(:,1,:)  amp_wn(:,1,:)]; %2nd dimension is 1 = peakresp or 2 = amp_wn
allEv_stat = squeeze(nanmean(allEv1,2));
allEv2 = [peakresp(:,2,:)  amp_wn(:,2,:)];allEv_mv = squeeze(nanmean(allEv2,2));

peakresp_mn = squeeze(mean(peakresp,2)); amp_wn_mn = squeeze(mean(amp_wn,2));
allEv_pre = [peakresp_mn(:,1) amp_wn_mn(:,1)]; allEv_pre = mean(allEv_pre,2);
allEv_post = [peakresp_mn(:,2) amp_wn_mn(:,2)]; allEv_post = mean(allEv_pre,2);
allEv = [allEv_pre allEv_post];

titles = {'Saline','DOI'}
figure
for t=1:2
    set(gcf,'Name', 'drift_wn_mean_mv_mean');
    used = goodAll==1 & treatment==t & (hasWn==1 |hasDrift==1);
    subplot(1,2,t)
    plot(allEv_mv(used &inhAll==0,1),allEv_mv(used&inhAll==0,2),'.','Markersize',20);
    hold on;axis square;
    %xlim([-1 25]); ylim([-1 25]);
    plot([-30 30],[-30 30]);title(titles{t})
    plot(allEv_mv(used & inhAll==1,1),allEv_mv(used &inhAll==1,2),'.r','Markersize',20);
    n_cells(t) = sum(used);
    text(-20, 33, ['n = ' num2str(n_cells(t))],'FontSize',20)
    for i = 0:1
        pre = allEv_mv(used & inhAll==i,1); post = allEv_mv(used & inhAll==i,2);
        c = corrcoef(pre,post,'rows','pairwise');c = c(2,1);
        mdl= fitlm(pre,post);
        rsquared(t) = mdl.Rsquared.Ordinary;
        if i ==0, text(-25, 29, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20);
            text(-25, 25, ['c.c = ' num2str(c,'%.2f')],'FontSize',20);
        else text(-25, 21, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20);
            text(-25, 17, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
    end
end

%%
%peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont
peakresp = squeeze(mean(drift_orient,2))-squeeze(drift_spont);
amp_wn = wn_evoked-wn_spont;

peakresp_mn = squeeze(mean(peakresp,2)); amp_wn_mn = squeeze(mean(amp_wn,2));
allEv_pre = [peakresp_mn(:,1) amp_wn_mn(:,1)]; allEv_pre = mean(allEv_pre,2);
allEv_post = [peakresp_mn(:,2) amp_wn_mn(:,2)]; allEv_post = mean(allEv_post,2);
allEv = [allEv_pre allEv_post];

titles = {'Saline','DOI','5HT'}
figure
for t=1:2
    set(gcf,'Name', 'drift_wn_mean of stat and running');
    used = goodAll==1 & treatment==t & (hasWn==1 |hasDrift==1);
    subplot(1,2,t)
    plot(allEv(used &inhAll==0,1),allEv(used&inhAll==0,2),'.','Markersize',20);
    hold on;axis square;
    %xlim([-1 25]); ylim([-1 25]);
    plot([-30 30],[-30 30]);title(titles{t})
    plot(allEv(used & inhAll==1,1),allEv(used &inhAll==1,2),'.r','Markersize',20);
    n_cells(t) = sum(used);
    text(-20, 33, ['n = ' num2str(n_cells(t))],'FontSize',20)
    for i = 0:1
        pre = allEv(used & inhAll==i,1); post = allEv(used & inhAll==i,2);
        c = corrcoef(pre,post,'rows','pairwise');c = c(2,1);
        mdl= fitlm(pre,post);
        rsquared(t) = mdl.Rsquared.Ordinary;
        if i ==0, text(-25, 29, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20);
            text(-25, 25, ['c.c = ' num2str(c,'%.2f')],'FontSize',20);
        else text(-25, 21, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20);
            text(-25, 17, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
    end
end

%%
clear data
titles={'Saline', 'DOI', '5HT'}
figure
set(gcf,'Name','mean dark FR')

for t=1:2
    used =  goodAll==1 & treatment==t;
    subplot(1,2,t)
    plot(meanRdark(inhAll==0 & used,1), meanRdark(inhAll==0 & used,2),'.','Markersize',20);hold on; axis square;
    set(gca,'FontSize',18)
    plot(meanRdark(inhAll==1  & used,1), meanRdark(inhAll==1 &used,2),'r.','Markersize',20);hold on; axis square;
    xlim([0 25]);ylim([0 25]); %xlabel('Pre spikes/sec');ylabel('Post spikes/sec');
    mpre=nanmean(mean(meanRdark(used,1),2))
    mpost=nanmean(mean(meanRdark(used,2),2))
    plot(mpre,mpost,'+k','Markersize',12,'Linewidth',2)
    plot([0 30],[0 30],'Linewidth',2);
    n_cells(t) = sum(used); title(titles{t});
    text(2, 24, ['n = ' num2str(n_cells(t))],'FontSize',18)
    for i = 0:1
        pre = meanRdark(used & inhAll==i,1); post = meanRdark(used & inhAll==i,2);
        c = corrcoef(pre,post,'rows','pairwise'); c = c(2,1);
        mdl= fitlm(pre,post)
        rsquared(t) = mdl.Rsquared.Ordinary
        if i ==0, text(2, 22, ['r^2 = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
            text(2, 20, ['c.c = ' num2str(c,'%.2f')],'FontSize',20)
        else text(2, 18, ['r^2 inh = ' num2str(rsquared(t),'%.2f')],'FontSize',20)
            text(2, 16, ['c.c inh = ' num2str(c,'%.2f')],'FontSize',20);end
    end
end

%%

use = ~isnan(meanRdark(:,1))
for t=1:3
    subplot(1,3,t)
    plot(meanRdark(~inhAll & treatment==t,1), meanRdark(~inhAll &treatment==t,2),'.','Markersize',12);hold on; axis square;
    set(gca,'FontSize',18)
    plot(meanRdark(inhAll  & treatment==t,1), meanRdark(inhAll &treatment==t,2),'r.','Markersize',12);hold on; axis square;
    
    
    allPre = meanRdark(goodAll==1 & (use==1)' & treatment==t,1)
    allPost = meanRdark(goodAll==1 & (use==1)' & treatment==t,2)
    c = corrcoef(allPre,allPost)
    c=c(2)
    %     mdl_dark= fitlm(meanRdark(treatment==t,1), meanRdark(treatment==t,2),2)
    %     rsquared_dark(t) = mdl_dark.Rsquared.Ordinary
    %xlim([0 20]);ylim([0 20]); %xlabel('Pre spikes/sec');ylabel('Post spikes/sec');
    mpre=nanmean(mean(meanRdark(goodAll==1 & treatment==t,1),2))
    mpost=nanmean(mean(meanRdark(goodAll==1 & treatment==t,2),2))
    plot(mpre,mpost,'+k','Markersize',12,'Linewidth',2)
    plot([0 30],[0 30],'Linewidth',2);
    n_cells(t) = sum(goodAll==1 &treatment==t)
    text(2, 28, ['n = ' num2str(n_cells(t))],'FontSize',18)
    %     text(2, 25, ['r^2 = ' num2str(rsquared_dark(t))],'FontSize',18)
    text(2, 25, ['cc = ' num2str(c)],'FontSize',18)
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


%%
titles ={'Saline', 'DOI','5HT'};
figure
set(gcf,'Name','cv2 wn');
for t=1:2
    subplot(1,2,t)
    mdl_cv2wn= fitlm(cv2Wn(goodAll==1&treatment==t,1), cv2Wn(treatment==t & goodAll==1,2),2)
    rs_cv2wn(t) = mdl_cv2wn.Rsquared.Ordinary
    plot(cv2Wn(treatment==t& goodAll==1,1),cv2Wn(treatment==t& goodAll==1,2),'.','Markersize',18);axis square;hold on
    title(titles{t});
    plot([0 2], [0 2]);ylim([.5 2]);xlim([.5 2]);
    n_cells(t) = sum(treatment==t& goodAll==1)
    text(.55, 1.75, ['r^2 = ' num2str(rs_cv2wn(t))],'FontSize',18)
    text(.55, 1.5, ['n = ' num2str(n_cells(t))],'FontSize',18)
end

%%
figure
set(gcf,'Name','cv2 drift');
for t=1:2
    subplot(1,2,t)
    mdl_cv2drift= fitlm(cv2Drift(treatment==t&goodAll==1,1), cv2Drift(treatment==t&goodAll==1,2),2)
    rs_cv2drift(t) = mdl_cv2drift.Rsquared.Ordinary
    plot(cv2Drift(treatment==t&goodAll==1,1),cv2Drift(treatment==t&goodAll==1,2),'.','Markersize',18);axis square;hold on
    plot([0 2], [0 2]);ylim([.5 2]);xlim([.5 2]);
    title(titles{t});
    n_cells(t) = sum(treatment==t& goodAll==1)
    text(.55, 1.85, ['r^2 = ' num2str(rs_cv2drift(t))],'FontSize',18)
    text(.55, 1.65, ['n = ' num2str(n_cells(t))],'FontSize',18)
end
%%
figure
set(gcf,'Name','cv2 dark');
for t=1:2
    subplot(2,2,t)
    mdl_cv2dark= fitlm(cv2Dark(treatment==t&goodAll==1,1), cv2Dark(treatment==t&goodAll==1,2),2)
    rs_cv2dark(t) = mdl_cv2dark.Rsquared.Ordinary
  %  plot(cv2Dark(treatment==t,1),cv2Dark(treatment==t,2),'.','Markersize',18);axis square;hold on

  hist(cv2Dark(treatment==t,1));xlim([0 2]);
 subplot(2,2,t+2)
  hist(cv2Dark(treatment==t,2))

   % plot([0 2], [0 2]);
 %ylim([.5 2]);
xlim([0 2]);
    title(titles{t});
    n_cells(t) = sum(treatment==t& goodAll==1)
   % text(.55, 1.85, ['r^2 = ' num2str(rs_cv2dark(t))],'FontSize',18)
   % text(.55, 1.65, ['n = ' num2str(n_cells(t))],'FontSize',18)
end

%%
figure
set(gcf,'Name','cv2 drift');
for t=1:3
    subplot(1,3,t)
    mdl_cv2drift= fitlm(cv2Drift(treatment==t&goodAll==1,1), cv2Drift(treatment==t&goodAll==1,2),2)
    rs_cv2drift(t) = mdl_cv2drift.Rsquared.Ordinary
    plot(cv2Drift(treatment==t&goodAll==1,1),cv2Drift(treatment==t&goodAll==1,2),'.','Markersize',18);axis square;hold on
    plot([0 2], [0 2]);ylim([.5 2]);xlim([.5 2]);
    title(titles{t});
    n_cells(t) = sum(treatment==t& goodAll==1)
    text(.55, 1.85, ['r^2 = ' num2str(rs_cv2drift(t))],'FontSize',18)
    text(.55, 1.65, ['n = ' num2str(n_cells(t))],'FontSize',18)
end

figure
set(gcf,'Name','cv2 dark');
for t=1:3
    subplot(1,3,t)
    mdl_cv2dark= fitlm(cv2Dark(treatment==t&goodAll==1,1), cv2Dark(treatment==t&goodAll==1,2),2)
    rs_cv2dark(t) = mdl_cv2dark.Rsquared.Ordinary
    plot(cv2Dark(treatment==t,1),cv2Dark(treatment==t,2),'.','markersize',18);axis square;hold on
    plot([0 2], [0 2]);ylim([.5 2]);xlim([.5 2]);
    title(titles{t});
    n_cells(t) = sum(treatment==t& goodAll==1)
    text(.55, 1.85, ['r^2 = ' num2str(rs_cv2dark(t))],'FontSize',18)
    text(.55, 1.65, ['n = ' num2str(n_cells(t))],'FontSize',18)
end

%%

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
%
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

%% correlation for wn %%
figure
for t=1:2
    preWnCorr = wnCorr(corrTreatment==t,1)
    postWnCorr=wnCorr(corrTreatment==t,2)
    mdl = fitlm(preWnCorr,postWnCorr)
    rsquared(t) = mdl.Rsquared.Adjusted
    cc = corrcoef(preWnCorr,postWnCorr,'rows','pairwise'); cc=cc(2,1);
    cc(t) = cc;
    
    subplot(2,2,t);
    plot(wnCorr(corrTreatment==t,1),wnCorr(corrTreatment==t,2),'.','Markersize',12); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{t});
    xlabel('pre wn corr'); ylabel('post')
    set(gcf,'Name','Wn Corr')
    text(-.25, .85, ['r^2 = ' num2str(rsquared(t),'%.2f')])
    text(-.25, .95, ['cc = ' num2str(cc(t),'%.2f')])
    h= hist(wnCorr(corrTreatment==t,1),-1:.1:1);
    h2=hist(wnCorr(corrTreatment==t,2),-1:.1:1);
    subplot(2,2,t+2);
    plot(h,'-b');hold on; axis square;
    plot(h2,'-r')
    %xlim([-1 1]);
    %title(titles{i});
end
%% corr for drift
figure
for t=1:2
    predriftCorr = driftCorr(corrTreatment==t,1)
    postdriftCorr = driftCorr(corrTreatment==t,2)
    mdl = fitlm(predriftCorr,postdriftCorr)
    rsquared(t) = mdl.Rsquared.Adjusted
    cc = corrcoef(predriftCorr,postdriftCorr,'rows','pairwise'); cc=cc(2,1);
    cc(t) = cc;
    
    subplot(2,2,t);
    plot(driftCorr(corrTreatment==t,1),driftCorr(corrTreatment==t,2),'.','Markersize',12); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{t});
    xlabel('pre drift corr'); ylabel('post')
    set(gcf,'Name','Drift Corr')
    text(-.25, .85, ['r^2 = ' num2str(rsquared(t),'%.2f')])
    text(-.25, .95, ['cc = ' num2str(cc(t),'%.2f')])
    h= hist(driftCorr(corrTreatment==t,1),-1:.1:1);
    h2=hist(driftCorr(corrTreatment==t,2),-1:.1:1);
    subplot(2,2,t+2);
    plot(h,'-b');hold on; axis square;
    plot(h2,'-r')
    %xlim([-1 1]);
    %title(titles{t});
end

%% corr for dark
figure
for t=1:2
    predarkCorr = darkCorr(corrTreatment==t,1)
    postdarkCorr=darkCorr(corrTreatment==t,2)
    mdl = fitlm(predarkCorr,postdarkCorr)
    rsquared(t) = mdl.Rsquared.Adjusted
    cc = corrcoef(predarkCorr,postdarkCorr,'rows','pairwise'); cc=cc(2,1);
    cc(t) = cc;
    subplot(2,2,t);
    plot(darkCorr(corrTreatment==t,1),darkCorr(corrTreatment==t,2),'.'); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{t});
    xlabel('pre wn corr'); ylabel('post')
    set(gcf,'Name','Wn Corr')
    text(-.25, .85, ['r^2 = ' num2str(rsquared(t),'%.2f')])
    text(-.25, .95, ['cc = ' num2str(cc(t),'%.2f')])
    h= hist(darkCorr(corrTreatment==t,1),-1:.1:1);
    h2=hist(darkCorr(corrTreatment==t,2),-1:.1:1);
    subplot(2,2,t+2);
    plot(h,'-b');hold on; axis square;
    plot(h2,'-r')
    %xlim([-1 1]);
    title(titles{t});
end

%% sta exp var
titles={'Saline','DOI','5HT'};
figure
for t=1:2
    subplot(1,2,t)
    mdl_expV = fitlm(sta_exp_var(data_wn &(treatment==t),1),sta_exp_var(data_wn&(treatment==t),2))
    rsquared_expV(t) = mdl_expV.Rsquared.Ordinary
    plot(sta_exp_var(data_wn & treatment==t&inhAll==0,1),sta_exp_var(data_wn &treatment==t&inhAll==0,2),'.','MarkerSize',20)
    %plot(sta_exp_var(useSTA & (treatment==t)',1),sta_exp_var(useSTA &(treatment==t)',2),'.')
    hold on; plot([0 1],[0 1]);
    plot(sta_exp_var(data_wn & treatment==t & inhAll==1,1),sta_exp_var(data_wn &treatment==t&inhAll==1,2),'.r','MarkerSize',20)
    
    axis square
    title(titles{t}); xlabel('Pre exp var');ylabel('Post exp var');
    text(.02, .85, ['r^2 = ' num2str(rsquared_expV(t),'%.2f')],'FontSize',18)
    c = corrcoef(sta_exp_var(data_wn & treatment==t,1),sta_exp_var(data_wn &treatment==t,2));
    c= c(2,1);
    n_cells(t) = sum(data_wn & treatment==t)
    text(.02 ,.95, ['n = ' num2str(n_cells(t))],'FontSize',18)
    text(.02 ,.75, ['cc = ' num2str(c,'%.2f')],'FontSize',18)
    
end

%% nx & ny prepost
titles={'Saline','DOI','5HT'};


clear max_wn
% low_wn = squeeze(min(wn_crf,[],2));
amp_wn = wn_evoked-wn_spont;%
useResp_wn = amp_wn(:,1,1)>0 | amp_wn(:,1,2)>0; %| amp_wn(:,2,1)>3 | amp_wn(:,2,2)>3;
data_wn = goodAll==1 &hasWn==1& useResp_wn';

useSTA = sta_exp_var(:,1)>.6 & sta_exp_var(:,2)>.6
useSTA = useSTA & data_wn'% &(layerAll==4)'
for t=1:2 %5ht doesnt show up with exp var smaller than .4
    figure
    mdl_nx = fitlm(sta_nx(useSTA & (treatment==t)',1),sta_nx(useSTA&(treatment==t)',2))
    mdl_ny= fitlm(sta_ny(useSTA & (treatment==t)',1),sta_ny(useSTA&(treatment==t)',2))
    rsquared_nx(t) = mdl_nx.Rsquared.Ordinary
    rsquared_ny(t) = mdl_ny.Rsquared.Ordinary
    c_nx = corrcoef(sta_nx(useSTA & (treatment==t)',1),sta_nx(useSTA&(treatment==t)',2));c_nx=c_nx(2,1);
    c_ny = corrcoef(sta_ny(useSTA & (treatment==t)',1),sta_ny(useSTA&(treatment==t)',2));c_ny=c_ny(2,1);
    subplot(1,2,1)
    plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==t)',2),'.','Markersize',20);
    hold on;
    plot(sta_nx(useSTA & (inhAll==1)'&(treatment==t)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==t)',2),'.r','Markersize',20);
    axis square;xlim([0 .7]); ylim([0 .7]); set(gca,'FontSize',18);
    text(.02, .63, ['r^2 = ' num2str(rsquared_nx(t),'%.2f')],'FontSize',18)
    text(.02, .53, ['cc = ' num2str(c_nx,'%.2f')],'FontSize',18)
    
    mpre=nanmean(mean(sta_nx(useSTA & (treatment==t)',1)))
    mpost=nanmean(mean(sta_nx(useSTA & (treatment==t)',2)))
    plot(mpre,mpost,'pg','Markersize',10);hold on
    title(titles{t},'FontSize',25); xlabel('Pre nx');ylabel('Post nx');
    plot([0 1],[0 1]);
    % subplot(2,2,t+2)
    subplot(1,2,2)
    plot(sta_ny(useSTA & (inhAll==0)'& (treatment==t)',1),sta_ny(useSTA &(inhAll==0)'&(treatment==t)',2),'.','Markersize',20);
    hold on;
    plot(sta_ny(useSTA & (inhAll==1)'& (treatment==t)',1),sta_ny(useSTA &(inhAll==1)'&(treatment==t)',2),'.r','Markersize',20);
    axis square;xlim([0 .7]);ylim([0 .7]);set(gca,'FontSize',18);
    hold on; plot([0 1],[0 1]);
    text(.02, .55, ['r^2 = ' num2str(rsquared_ny(t),'%.2f')],'FontSize',18)
    text(.02, .45, ['cc = ' num2str(c_ny,'%.2f')],'FontSize',18)
    
    mpre=nanmean(mean(sta_ny(useSTA  & (treatment==t)',1)))
    mpost=nanmean(mean(sta_ny(useSTA & (treatment==t)',2)))
    plot(mpre,mpost,'pg','Markersize',10)
    title(titles{t},'FontSize',25); xlabel('Pre ny');ylabel('Post ny');
    n_cells(t) = sum(useSTA &(treatment==t)')
    text(.02 ,.65, ['n = ' num2str(n_cells(t))],'FontSize',18)
end

%%
figure
set(gcf,'Name', 'nx prepost all layers')
subplot(2,3,1)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==1)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==1)' & (layerAll==2|3)',2),'.','Markersize',20);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==1)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==1)' & (layerAll==2|3)',2),'r.','Markersize',20);
title('Layer 2/3');xlabel('pre nx');ylabel('post nx');
subplot(2,3,2)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==1)' & (layerAll==4)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==1)' & (layerAll==4)',2),'.','Markersize',20);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==1)' & (layerAll==4)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==1)' & (layerAll==4)',2),'r.','Markersize',20);
title('Layer 4');xlabel('pre nx');ylabel('post nx');
subplot(2,3,3)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==1)' & (layerAll==5)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==1)' & (layerAll==5)',2),'.','Markersize',20);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==1)' & (layerAll==5)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==1)' & (layerAll==5)',2),'r.','Markersize',20);
title('Layer 5');xlabel('pre nx');ylabel('post nx');
subplot(2,3,4)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==2)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==2)' & (layerAll==2|3)',2),'.','Markersize',20);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==2)' & (layerAll==2|3)',2),'r.','Markersize',20);
xlabel('pre nx');ylabel('post nx');
subplot(2,3,5)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==2)' & (layerAll==4)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==2)' & (layerAll==4)',2),'.','Markersize',20);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==4)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==2)' & (layerAll==4)',2),'r.','Markersize',20);
xlabel('pre nx');ylabel('post nx');
subplot(2,3,6)
plot(sta_nx(useSTA & (inhAll==0)' & (treatment==2)' & (layerAll==5)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==2)' & (layerAll==5)',2),'.','Markersize',20);
hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==5)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==2)' & (layerAll==5)',2),'r.','Markersize',20);
xlabel('pre nx');ylabel('post nx');
% subplot(2,3,7)
% plot(sta_nx(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==2|3)',2),'.','Markersize',20);
% hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==2|3)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==2|3)',2),'r.','Markersize',20);
% xlabel('pre nx');ylabel('post nx')
% subplot(2,3,8)
% plot(sta_nx(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==4)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==4)',2),'.','Markersize',20);
% hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==4)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==4)',2),'r.','Markersize',20);
% xlabel('pre nx');ylabel('post nx');
% subplot(2,3,9)
% plot(sta_nx(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==5)',1),sta_nx(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==5)',2),'.','Markersize',20);
% hold on;plot([0 1],[0 1]); plot(sta_nx(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==5)',1),sta_nx(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==5)',2),'r.','Markersize',20);
% xlabel('pre nx');ylabel('post nx');

%%

figure
set(gcf,'Name', 'ny prepost all layers')
subplot(2,3,1)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==1)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==1)' & (layerAll==2|3)',2),'.','Markersize',20);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==1)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==1)' & (layerAll==2|3)',2),'r.','Markersize',20);
title('Layer 2/3');xlabel('pre ny');ylabel('post ny')
subplot(2,3,2)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==1)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==1)' & (layerAll==4)',2),'.','Markersize',20);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==1)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==1)' & (layerAll==4)',2),'r.','Markersize',20);
title('Layer 4');xlabel('pre ny');ylabel('post ny')
subplot(2,3,3)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==1)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==1)' & (layerAll==5)',2),'.','Markersize',20);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==1)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==1)' & (layerAll==5)',2),'r.','Markersize',20);
title('Layer 5');xlabel('pre ny');ylabel('post ny')
subplot(2,3,4)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==2)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==2)' & (layerAll==2|3)',2),'.','Markersize',20);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==2)' & (layerAll==2|3)',2),'r.','Markersize',20);
xlabel('pre ny');ylabel('post ny')
subplot(2,3,5)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==2)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==2)' & (layerAll==4)',2),'.','Markersize',20);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==2)' & (layerAll==4)',2),'r.','Markersize',20);
xlabel('pre ny');ylabel('post ny')
subplot(2,3,6)
plot(sta_ny(useSTA & (inhAll==0)' & (treatment==2)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==2)' & (layerAll==5)',2),'.','Markersize',20);
hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==2)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==2)' & (layerAll==5)',2),'r.','Markersize',20);
xlabel('pre ny');ylabel('post ny')
%subplot(2,3,7)
% plot(sta_ny(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==2|3)',2),'.','Markersize',20);
% hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==2|3)',2),'r.','Markersize',20);
% xlabel('pre ny');ylabel('post ny')
% subplot(2,3,8)
% plot(sta_ny(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==4)',2),'.','Markersize',20);
% hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==4)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==4)',2),'r.','Markersize',20);
% xlabel('pre ny');ylabel('post ny')
% subplot(2,3,9)
% plot(sta_ny(useSTA & (inhAll==0)' & (treatment==3)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==3)' & (layerAll==5)',2),'.','Markersize',20);
% hold on;plot([0 1],[0 1]); plot(sta_ny(useSTA & (inhAll==1)' & (treatment==3)' & (layerAll==5)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==3)' & (layerAll==5)',2),'r.','Markersize',20);
% xlabel('pre ny');ylabel('post ny')

%%
useSTA = sta_exp_var(:,1)>.6 & sta_exp_var(:,2)>.6
useSTA = useSTA & data_wn'% &(layerAll==4)'
figure
titles={'Saline', 'DOI','5HT'}
for t=1:2
    subplot(2,2,t)
    set(gcf,'Name', 'layer 2/3 prepost nx and ny')
    mdl_nxny_pre = fitlm(sta_nx(useSTA & (treatment==t)' &(layerAll==2|3)',1),sta_ny(useSTA&(treatment==t)'&(layerAll==2|3)',1))
    rsquared_nxny_pre(t) = mdl_nxny_pre.Rsquared.Ordinary
    mdl_nxny_post = fitlm(sta_nx(useSTA & (treatment==t)' &(layerAll==2|3)',2),sta_ny(useSTA&(treatment==t)'&(layerAll==2|3)',2))
    rsquared_nxny_post(t) = mdl_nxny_post.Rsquared.Ordinary
    plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==t)' & (layerAll==2|3)',1),'.','Markersize',20); hold on;
    xlabel('nx pre');ylabel('ny pre')
    plot(sta_nx(useSTA & (inhAll==1)' & (treatment==t)' & (layerAll==2|3)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==t)' & (layerAll==2|3)',1),'r.','Markersize',20);
    hold on;plot([0 1],[0 1]); axis square;xlim([0 1])
    text(.02, .8, ['r^2 = ' num2str(rsquared_nxny_pre(t))],'FontSize',18)
    subplot(2,2,t+2)
    plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)' & (layerAll==2|3)',2),sta_ny(useSTA&(inhAll==0)'&(treatment==t)' & (layerAll==2|3)',2),'.','Markersize',20); hold on;
    xlabel('nx post');ylabel('ny post')
    plot(sta_nx(useSTA & (inhAll==1)' & (treatment==t)' & (layerAll==2|3)',2),sta_ny(useSTA&(inhAll==1)'&(treatment==t)' & (layerAll==2|3)',2),'r.','Markersize',20);
    plot([0 1],[0 1]);axis square;xlim([0 1])
    title(titles{t});
    text(.02, .8, ['r^2 = ' num2str(rsquared_nxny_post(t))],'FontSize',18)
end

%%
figure
titles={'Saline', 'DOI','5HT'}
for t=1:2
    subplot(2,2,t)
    set(gcf,'Name', 'all layers prepost nx and ny')
    mdl_nxny_pre = fitlm(sta_nx(useSTA & (treatment==t)',1),sta_ny(useSTA&(treatment==t)',1))
    rsquared_nxny_pre(t) = mdl_nxny_pre.Rsquared.Ordinary
    plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)',1),sta_ny(useSTA&(inhAll==0)'&(treatment==t)' ,1),'.','Markersize',20);
    hold on;plot([0 1],[0 1]); axis square; xlabel('pre nx');ylabel('pre ny');
    text(.02, .8, ['r^2 = ' num2str(rsquared_nxny_pre(t))],'FontSize',18)
    plot(sta_nx(useSTA & (inhAll==1)' & (treatment==t)',1),sta_ny(useSTA&(inhAll==1)'&(treatment==t)' ,1),'r.','Markersize',20);
    subplot(2,2,t+2)
    mdl_nxny_post = fitlm(sta_nx(useSTA & (treatment==t)',1),sta_ny(useSTA&(treatment==t)',2))
    rsquared_nxny_post(t) = mdl_nxny_post.Rsquared.Ordinary
    plot(sta_nx(useSTA & (inhAll==0)' & (treatment==t)' ,2),sta_ny(useSTA&(inhAll==0)'&(treatment==t)' ,2),'.','Markersize',20);
    hold on;plot([0 1],[0 1]);xlim([0 1]);axis square; xlabel('post nx');ylabel('post ny');
    plot(sta_nx(useSTA & (inhAll==1)' & (treatment==t)' ,2),sta_ny(useSTA&(inhAll==1)'&(treatment==t)' ,2),'r.','Markersize',20);
    text(.02, .8, ['r^2 = ' num2str(rsquared_nxny_post(t))],'FontSize',18)
    title(titles{t});
end
%%


useSTA = find((sta_exp_var(:,1)>.6 | sta_exp_var(:,2)>.6) &(treatment==doi)');
for prepost=1:2
    figure
    if prepost==1
        set(gcf,'Name','pre')
    else
        set(gcf,'Name','post');end
    for c=8%1:length(useSTA)
        %subplot(10,7,c)
        imagesc(sta_all_img{useSTA(c),prepost});axis square;
        %imagesc(sta_all_fit{useSTA(c),prepost});axis square
        set(gca,'xticklabel',{[]},'yticklabel',{[]})
    end
end
%%
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

%%

amp_wn = wn_evoked-wn_spont;%
useResp_wn = amp_wn(:,2,1)>1 | amp_wn(:,2,2)>1; %| amp_wn(:,2,1)>3 | amp_wn(:,2,2)>3;
data_wn = goodAll==1 &hasWn==1& useResp_wn';

clear useSTA
thing = {'saline','doi'};
mycol = {'k','r'};
grpRF = {}
figure;
for t=1
%useSTA = find((sta_exp_var(:,1)>.6 & sta_exp_var(:,2)>.6) &(treatment==doi)' )%&(layerAll==2|3)');
useSTA = find(treatment==t &data_wn)%&(layerAll==2|layerAll==3));
% useSTA = find(treatment==saline &data_wn& layerAll==1);
RFcorr = [];
for c=1:length(useSTA)
    RFcorr(c) = corr2(sta_all_img{useSTA(c),1},sta_all_img{useSTA(c),2});
end
grpRF{t} = RFcorr;
h(t,:)=hist(RFcorr, -1:.2:1);
hold on
plot(-1:0.2:1,(h(t,:)/length(useSTA)));%ylim([0 13])
axis square
end
legend(thing)

[h0,p,kstat] = kstest2(h(1,:),h(2,:)) %saline vs doi
[h0,p] = ttest2(grpRF{1},grpRF{2})% ttest
%% plot white noise response functions for all units

titles = {'layer 1','layer 2','layer 3','layer 4','layer 5','layer 6'};
C = {[1 0 0],[.5 0 0]} %bright = pre
D = {[0 1 0],[0 .5 0]}
for t=1:2
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

%% mean CRF all layers
for t=1:2
    figure
    if t == 1
        set(gcf,'Name', 'mean SALINE CRF');
    else t==2
        set(gcf,'Name', 'mean DOI CRF');
    end
    use = data_wn & treatment==t & hasWn==1
    for prepost = 1:2
        mn_wn = squeeze(mean(wn_crf(use,:,1,prepost),1)); %stationary prepost
        sem_stat = std(wn_crf(use,:,1,prepost),1)/sqrt(sum(use))
        mn_wn_mv = squeeze(mean(wn_crf(data_wn & treatment==t,:,2,prepost),1));
        sem_mv = std(wn_crf(use,:,2,prepost),1)/sqrt(sum(use))
        
        
        %ncells = sum(use);xlim([0 2.5]);
        if prepost ==1
            err_stat = shadedErrorBar(1:20,mn_wn,sem_stat,'b',1)
            hold on;
            err_mv = shadedErrorBar(1:20,mn_wn_mv,sem_mv,'g',1)
        else
            err_stat = shadedErrorBar(1:20,mn_wn,sem_stat,'r',1) %post stat
            hold on;
            err_mv = shadedErrorBar(1:20,mn_wn_mv,sem_mv,'k',1) % post mv
        end
    end
end

% titles={'Saline','DOI','5HT'};
for t=1:2
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
        bar(Mbins,h/sum(useN),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .6]);xlim([-1.5 1.5]); axis square
        title(layerz{l});
    end
end

%% to do: seperate to have multiple figures for each treatment
for t = 1:2
    figure
    if t==1, set(gcf,'Name','saline wn CRF'),
    else t==2, set(gcf,'Name','doi wn CRF'),end
    
    useN = find(goodAll & treatment==t & hasWn==1)%& layerAll==4)
    for i = 1:ceil(length(useN))/6
        np = ceil(sqrt(length(useN)/6));
        subplot(np,np,i);
        hold on
        plot(wn_crf(useN(i),:,1,1),'Color',[0.5 0 0]);  plot(wn_crf(useN(i),:,2,1),'Color',[0 0.5 0]);
        plot(wn_crf(useN(i),:,1,2),'Color',[1 0 0]);  plot(wn_crf(useN(i),:,2,2),'Color',[0 1 0]);
        yl = get(gca,'Ylim'); ylim([0 max(yl(2),10)]);
        if inhAll(useN(i)) ==1 , xlabel('inh'); else  xlabel('exc');
        end
    end
end
%%
%%%%%%%%%%%%%%%%%%% drift %%%%%%%%%%%%%%%%%%%%%%%%%%
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
    [respmax oind] = max(mean(drift_sf(c ,:,:,1),3)); %pre
    % [g h]=(max(respmax(:)));
    prefSF(c)=oind;
    %end
end

useN =dataAll;
for c= 1:length(useN)
    % for prepost =1:2
    [respmax oind] = max(mean(drift_sf(c ,:,:,2),3)); %post
    % [g h]=(max(respmax(:)));
    prefSFPost(c)=oind;
    %end
end



%%% drift_cond_tcourse(cells,move,prepost,orient,sf,time);
for c=1:length(useN)
    tcourse_pref(c,:,:,:) = squeeze(nanmean(drift_cond_tcourse(c,:,:,prefOri(c),:,:),5));
    tcourse_orth(c,:,:,:) = squeeze(nanmean(drift_cond_tcourse(c,:,:,mod(prefOri(c)+3-1,12)+1,:,:),5));
    tcourse_all(c,:,:,:) = squeeze(mean(nanmean(drift_cond_tcourse(c,:,:,:,:,:),5),4));
end

clear max amp low
peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont
useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2;
data = goodAll & useResp' & ~inhAll & hasDrift==1 ;%which cells to include
dataInh =goodAll & useResp' & inhAll & hasDrift==1 ;

%%

%PSTH for each cell's preferred orientation - stationary
titles={'Saline', 'DOI'};

for i= 1:3
    figure
    for t=1:2
        if i==1
            used = data & (layerAll==5) & treatment==t;
            set(gcf,'Name','grating lyr5 stationary');
            mn = squeeze(mean(tcourse_pref(used,1,:,:),1))';
            sem = std(tcourse_pref(used,1,:,:),1)...
                /sqrt(sum(used))
            sem=squeeze(sem);sem=sem';
            ncells = sum(used);xlim([0 2.5]);
            xlim([0 2.5]);%ylim([0 8]);
        elseif i==2
            set(gcf,'Name','grating lyr 2/3 stationary');
            used = data & (layerAll==2 |layerAll==3) & treatment==t;
            mn = squeeze(mean(tcourse_pref(used,1,:,:),1))';
            sem = std(tcourse_pref(used,1,:,:),1)...
                /sqrt(sum(used))
            sem=squeeze(sem);sem=sem';
            ncells = sum(used);xlim([0 2.5]);
            xlim([0 2.5]);%ylim([0 8]);
        elseif i ==4
            set(gcf,'Name','grating lyr 4 stationary');
            used = data & (layerAll==4) & treatment==t;
            
            mn = squeeze(mean(tcourse_pref(used,1,:,:),1))';
            sem = std(tcourse_pref(used,1,:,:),1)...
                /sqrt(sum(used));
            sem=squeeze(sem);sem=sem';
            ncells = sum(used);xlim([0 2.5]);
            %ylim([min(min(mn))-1 max(max(mn))+.5]);
            xlim([0 2.5]);
        elseif i==3
            set(gcf,'Name','grating inh stationary');
            used = dataInh  & treatment==t;
            
            mn = squeeze(mean(tcourse_pref(used,1,:,:),1))';
            sem = std(tcourse_pref(used,1,:,:),1)...
                /sqrt(sum(used));
            sem=squeeze(sem);sem=sem';
            ncells = sum(used);xlim([0 2.5]);
            xlim([0 2.5]);%ylim([0 16]);
        elseif i==5
            set(gcf,'Name','grating lyr6 stationary');
            used = data & layerAll==6  & treatment==t;
            
            mn = squeeze(mean(tcourse_pref(used,1,:,:),1))';
            sem = std(tcourse_pref(used,1,:,:),1)...
                /sqrt(sum(used));
            sem=squeeze(sem);sem=sem';
            ncells = sum(used);xlim([0 2.5]);
            % ylim([min(min(mn))-1 max(max(mn))+.5]);
            xlim([0 2.5]);
        end
        %         mn = mn - repmat(mn(1,:),[50 1]); %%% subtracts prestim spont from
        %        all time points
        mn = circshift(mn,10);
        sem = circshift(sem,10);
        subplot(1,2,t)
        x = ((1:length(mn)-5)*dt -dt/2);
        %         plot(x,mn(1:45,1),'LineWidth',2); hold on;set(gca,'fontsize', 18); hold on;
        %         plot(x,mn(1:45,2),'LineWidth',2);axis square; set(gca,'fontsize', 18);
        preerr=shadedErrorBar(x,mn(1:45,1),sem(1:45,1)','b',1); axis square;
        set(gca,'fontsize', 18); hold on;
        posterr=shadedErrorBar(x,mn(1:45,2),sem(1:45,2)','r',1);%axis square; set(gca,'fontsize', 18)
        text(1.75,max(max((mn))+.75), ['n = ' num2str(ncells)],'FontSize',20)
        title([titles{t} ' n = ' num2str(ncells) ]);
        xlabel('time (ms)');
        %   ylabel('spikes/sec');
        %ylim([0 max(mn(:))+1])
    end
end
%%
%PSTH for each cell's preferred orientation - moving
for i= 1:4
    figure
    for t=1:2
        if i==1
            set(gcf,'Name','grating lyr5 mv');
            mn = squeeze(nanmean(tcourse_pref(data & (hasDrift==1) & (layerAll==5) & treatment==t,2,:,:),1))';
            sem = nanstd(tcourse_pref(data & (hasDrift==1) & (layerAll==5) & treatment==t,2,:,:),1)...
                /sqrt(sum(data & (hasDrift==1) & (layerAll==5) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==5) & treatment==t);% xlim([0 2.5]);
            ylim([-2 7]);
        elseif i==2
            set(gcf,'Name','grating lyr 2/3 mv');
            mn = squeeze(nanmean(tcourse_pref(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,2,:,:),1))';
            sem = nanstd(tcourse_pref(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t,2,:,:),1)...
                /sqrt(sum(data & (hasDrift==1) & (layerAll==2 |layerAll==3) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==2| layerAll==3) & treatment==t);%xlim([0 2.5]);
            ylim([-2 10]);%xlim([0 2.5]);
        elseif i ==3
            set(gcf,'Name','grating lyr 4 mv');
            mn = squeeze(nanmean(tcourse_pref(data &(hasDrift==1) & layerAll==4 & treatment==t,2,:,:),1))';
            sem = nanstd(tcourse_pref(data &(hasDrift==1) & layerAll==4 & treatment==t,2,:,:),1)...
                /sqrt(sum(data &(hasDrift==1) & layerAll==4 & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==4) & treatment==t);%xlim([0 2.5]);
            %ylim([-.5 4]);xlim([0 2.5]);
        elseif i==4
            set(gcf,'Name','grating inh mv');
            mn = squeeze(nanmean(tcourse_pref(dataInh & (hasDrift==1) & treatment==t,2,:,:),1))';
            sem = nanstd(tcourse_pref(dataInh &(hasDrift==1) & treatment==t,2,:,:),1)...
                /sqrt(sum(dataInh &(hasDrift==1) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(dataInh & (hasDrift==1) & treatment==t);%xlim([0 2.5]);
            ylim([-.5 20]);%xlim([0 2.5]);
        else i==5
            set(gcf,'Name','grating lyr6 mv');
            mn = squeeze(nanmean(tcourse_pref(data & layerAll==6 & (hasDrift==1) & treatment==t,2,:,:),1))';
            sem = nanstd(tcourse_pref(data & layerAll==6 & (hasDrift==1) & treatment==t,2,:,:),1)...
                /sqrt(sum(data & layerAll==6 & (hasDrift==1) & treatment==t));
            sem=squeeze(sem);sem=sem';
            ncells = sum(data & (hasDrift==1) & (layerAll==6) & treatment==t);%xlim([0 2.5]);
            %  ylim([-.5 12]);xlim([0 2.5]);
        end
        mn = mn - repmat(mn(1,:),[50 1]); %%% subtracts prestim spont from
        %   all time points
        mn = circshift(mn,10);
        sem = circshift(sem,10);
        subplot(1,2,t)
        x = ((1:length(mn)-5)*dt -dt/2);
        %         plot(x,mn(1:45,1),'LineWidth',2); hold on;set(gca,'fontsize', 18); hold on;
        %         plot(x,mn(1:45,2),'LineWidth',2);axis square; set(gca,'fontsize', 18);
        preerr=shadedErrorBar(x,mn(1:45,1),(sem(1:45,1))','b',1); axis square; set(gca,'fontsize', 18); hold on;
        posterr=shadedErrorBar(x,mn(1:45,2),(sem(1:45,2))','r',1);%axis square; set(gca,'fontsize', 18)
        text(1.75,max(max((mn))+.75), ['n = ' num2str(ncells)],'FontSize',20)
        title(titles{t});
        xlabel('time (s)');
        %   ylabel('spikes/sec');
        %ylim([0 max(mn(:))+1])
    end
end




%%%%

% peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont
peakresp = squeeze(mean(drift_orient,2))-squeeze(drift_spont);%subtract driftspont

amp=peakresp;
useResp = amp(:,1,1)>1 | amp(:,1,2)>1;% |amp(:,2,1)>2 | amp(:,2,2)>2;
useResp = amp(:,1,1)>2 | amp(:,1,2)>2% |amp(:,2,1)>2 | amp(:,2,2)>2;
data = goodAll & useResp' & ~inhAll & hasDrift==1 ;%which cells to include
dataInh =goodAll & useResp' & inhAll & hasDrift==1 ;
dataAll =goodAll & useResp' & hasDrift==1 ;


titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
dt = 0.05;
for i= 1:4
    figure
    for t=1:2
        if i==1
            set(gcf,'Name','grating lyr5');
            used =data& layerAll==5 & treatment==t;
            mn = squeeze(mean(mnPsth(used,:,:),1))';
            sem = std(mnPsth(used,:,:),1)...
                /sqrt(sum(used));
            sem=squeeze(sem);sem=sem';
            ncells = sum(used);xlim([0 2.5]);
            %    ylim([-1 1]); xlim([0 2.5]);
        elseif i==2
            set(gcf,'Name','grating lyr 2/3');
            used =data& (layerAll==2 |layerAll==3) & treatment==t;
            mn = squeeze(mean(mnPsth(used,:,:),1))';
            sem = std(mnPsth(used,:,:),1)...
                /sqrt(sum(used));
            sem=squeeze(sem);sem=sem';
            ncells = sum(used);xlim([0 2.5]);
            %  ylim([-.5 5]);xlim([0 2.5]);
        elseif i ==4
            set(gcf,'Name','grating lyr 4');
            used =data& (layerAll==4) & treatment==t;
            
            mn = squeeze(mean(mnPsth(used,:,:),1))';
            
            sem = std(mnPsth(used,:,:),1)...
                /sqrt(sum(used));
            sem=squeeze(sem);sem=sem';
            ncells = sum(used);xlim([0 2.5]);
            %  ylim([-.5 4]);xlim([0 2.5]);
        elseif i==3
            set(gcf,'Name','grating inh');
            used =dataInh&  treatment==t;
            
            mn = squeeze(mean(mnPsth(used,:,:),1))';
            
            sem = std(mnPsth(used,:,:),1)...
                /sqrt(sum(used));
            sem=squeeze(sem);sem=sem';
            ncells = sum(used);xlim([0 2.5]);
            %  ylim([-.5 12]);xlim([0 2.5]);
        elseif i==5
            set(gcf,'Name','grating lyr6');
            used =data& layerAll==6 & treatment==t;
            
            mn = squeeze(mean(mnPsth(used,:,:),1))';
            sem = std(mnPsth(used,:,:),1)...
                /sqrt(sum(used));
            sem=squeeze(sem);sem=sem';
            ncells = sum(used);xlim([0 2.5]);
            %  ylim([-.5 12]);xlim([0 2.5]);
        end
        %   mn = mn - repmat(mn(1,:),[50 1]); %%% subtracts prestim spont
        mn = circshift(mn,10);
        sem = circshift(sem,10);
        subplot(1,2,t)
        x = ((1:length(mn)-5)*dt -dt/2);
        %         plot(x,mn(1:45,1),'LineWidth',2); hold on;set(gca,'fontsize', 18); hold on;
        %         plot(x,mn(1:45,2),'LineWidth',2);axis square; set(gca,'fontsize', 18);
        preerr=shadedErrorBar(x,mn(1:45,1),sem(1:45,1)','b',1); axis square; set(gca,'fontsize', 18); hold on;
        posterr=shadedErrorBar(x,mn(1:45,2),sem(1:45,2)','r',1);%axis square; set(gca,'fontsize', 18)
        text(1.75,max(max((mn))+.75), ['n = ' num2str(ncells)],'FontSize',20);
        xlim([0 2.5]);
        %plot((1:length(mn)-5)*dt -dt/2,mn(1:45,:),'LineWidth',2);axis square; set(gca,'fontsize', 18);
        title(titles{t}); xlabel('time (s)');
        ylabel('spikes/sec');
        %ylim([0 max(mn(:))+1])
    end
end

%%

amp=peakresp;
useResp = amp(:,1,1)>2 | amp(:,1,2)>2;% |amp(:,2,1)>2 | amp(:,2,2)>2;
data = goodAll & useResp' & inhAll==0 & hasDrift==1 ;%which cells to include
dataInh =goodAll & useResp' & inhAll==1 & hasDrift==1 ;
dataAll =goodAll & useResp' & hasDrift==1 ;


%prepost mean psth, normalized to pre range
for i=1:3
    clear used
    figure
    for t=1:2
        if i==1, set(gcf,'Name','grating lyr5');
            used =data& layerAll==5 & treatment==t;
            mn = squeeze(mean(mnPsth(used,:,:),1))';
            %  range = max(mn(:,1)) - min(mn(:,1)); %set range w/pre data to normalize to
            mnpre = (mn(:,1) - min(mn(:,1))) / range; mnpre = circshift(mnpre,10);
            mnpost = (mn(:,2) - min(mn(:,2)))/range; mnpost = circshift(mnpost,10);
            sem = std(mnPsth(used,:,:),1)/sqrt(sum(used));
            sem=squeeze(sem);sem=sem';
            semrange1 = max(sem(:,1)) - min(sem(:,1));
            semrange2 = max(sem(:,2)) - min(sem(:,2));
            sempre = (sem(:,1) - min(sem(:,1))) / semrange1; sempre = circshift(sempre,10);
            sempost = (sem(:,2) - min(sem(:,2))) / semrange2; sempost = circshift(sempost,10);
            ncells = sum(used);
        elseif i ==2,set(gcf,'Name','grating lyr 2/3');
            used =data& (layerAll==2|layerAll==3) & treatment==t;
            mn = squeeze(mean(mnPsth(used,:,:),1))';
            range = max(mn(:,1))% - min(mn(:,1)); %set range w/pre data to normalize to
            mnpre = (mn(:,1) - min(mn(:,1))) / range; mnpre = circshift(mnpre,10);
            mnpost = (mn(:,2)-min(mn(:,2)))/range; mnpost = circshift(mnpost,10);
            sem = std(mnPsth(used,:,:),1)/sqrt(sum(used));sem=squeeze(sem);sem=sem';
            semrange1 = max(sem(:,1)) - min(sem(:,1)); semrange2 = max(sem(:,2)) - min(sem(:,2));
            sempre = (sem(:,1) - min(sem(:,1))) / semrange1; sempre = circshift(sempre,10);
            sempost = (sem(:,2) - min(sem(:,2))) / semrange2; sempost = circshift(sempost,10);
            ncells = sum(used);
        elseif i ==3, set(gcf,'Name','grating inh');
            used =dataInh& treatment==t;
            mn = squeeze(mean(mnPsth(used,:,:),1))';
            range = max(mn(:,1)) %- min(mn(:,1)); %set range w/pre data to normalize to
            mnpre = (mn(:,1) - min(mn(:,1))) / range; mnpre = circshift(mnpre,10);
            mnpost = (mn(:,2)-min(mn(:,2)))/range; mnpost = circshift(mnpost,10);
            sem = std(mnPsth(used,:,:),1)/sqrt(sum(used));
            sem=squeeze(sem);sem=sem';
            semrange1 = max(sem(:,1)) - min(sem(:,1)); semrange2 = max(sem(:,2)) - min(sem(:,2));
            sempre = (sem(:,1) - min(sem(:,1))) / semrange1; sempre = circshift(sempre,10);
            sempost = (sem(:,2) - min(sem(:,2))) / semrange2; sempost = circshift(sempost,10);
            ncells = sum(used);
        end
        subplot(1,2,t)
        x = ((1:length(mnpre)-5)*dt -dt/2); x2 = ((1:length(mnpost)-5)*dt -dt/2);
        preerr=shadedErrorBar(x,mnpre(1:45,1),sempre(1:45,1)','b',1); hold on;
        posterr=shadedErrorBar(x2,mnpost(1:45,1),sempost(1:45,1)','r',1);
        axis square; set(gca,'fontsize', 18); hold on; xlim([0 2.5])
        title(titles{t}); xlabel('time (s)'); ylabel('spikes/sec');
        text(1.75,max(max((mnpre))+.75), ['n = ' num2str(ncells)],'FontSize',20);
        
    end
end


%for l=1:6
ampmean = squeeze(mean(peakresp,2))
useResp = ampmean(:,1)>0 &ampmean(:,2)>0;
use = goodAll &  inhAll==0 & hasDrift==1& useResp' ;%which cells to include
useN = find(use & (layerAll==2 |layerAll==3) & treatment==t);
figure
%useN= find(useResp' & goodAll==1 & hasDrift & treatment==doi &layerAll==2)
length(useN)
for c=1:length(useN)/2
    % set (gcf,'Name','Lyr 2/3 psth/unit') % set up in loop
    subplot(11,11,c)
    plot(squeeze(mnPsth(useN(c),1,:)));hold on;
    plot(squeeze(mnPsth(useN(c),2,:)),'r')
end


%%


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

%%
useOsi_pre = drift_osi(:,1,1)>.3;
useOsi_post = drift_osi(:,1,2)>.3;
useOsi=(useOsi_pre & useOsi_post)'
useOsi_mn = squeeze(mean(drift_osi,2));
peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont
% peakresp = squeeze(mean(drift_orient,2))-squeeze(drift_spont);%subtract driftspont
useResp = peakresp(:,2,1)>2| peakresp(:,2,2)>2;
ampmean = squeeze(nanmean(peakresp,2));
% useResp = ampmean(:,1)>2 |ampmean(:,2)>2;

data = goodAll & useResp' & ~inhAll & hasDrift==1 &useOsi==1 ;%which cells to include
dataInh =goodAll & useResp' & inhAll & hasDrift==1 & useOsi==1;
dataAll =goodAll & useResp' & hasDrift==1 &useOsi==1;
titles={'Saline', 'DOI'};
for mv= 1:2
    clear R rs_pf
    figure
    if mv==1, set(gcf,'Name', 'stationary');
        useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2; set(gcf,'Name','stationary pref ori');
        
    else
        set(gcf,'Name', 'mv');
        useResp = peakresp(:,2,1)>2 | peakresp(:,2,2)>2; set(gcf,'Name','mv pre ori')
        
    end
    for t=1:2
        data = goodAll & useResp' & ~inhAll & hasDrift==1 &useOsi==1&treatment==t ;%which cells to include
        dataInh =goodAll & useResp' & inhAll & hasDrift==1 & useOsi==1&treatment==t;
        dataAll =goodAll & useResp' & hasDrift==1 &useOsi==1&treatment==t;
        pre=drift_pref_theta(dataAll,mv,1); post= drift_pref_theta(dataAll,mv,2)
        post(pre-post > pi/2) = post(pre-post>pi/2)+pi; pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;
        
        mdl_pf= fitlm(pre,post)
        c=corrcoef(pre,post,'rows','pairwise'); c = c(2,1)
        rs_pf(t) = mdl_pf.Rsquared.Ordinary
        subplot(1,2,t)
         pre=drift_pref_theta(dataInh,mv,1); post= drift_pref_theta(dataInh,mv,2)
        post(pre-post > pi/2) = post(pre-post>pi/2)+pi;
        pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;
       
        plot(pre,post,'ko','Markersize',10,'linewidth',2,'MarkerFaceColor',[.9 .2 .5]);title(titles{t});hold on;
        pre=drift_pref_theta(data,mv,1); post= drift_pref_theta(data,mv,2)
        post(pre-post > pi/2) = post(pre-post>pi/2)+pi;
        pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;
        
        plot(pre,post,'ko','Markersize',10,'linewidth',2,'MarkerFaceColor',[.2 0.5 .8])
        axis square;plot([0 3.25], [0 3.25]);xlim([0 3.25]);ylim([0 3.25]);
        text(.2, 2.8, ['r^2 = ' num2str(rs_pf(t),'%.2f')],'FontSize',18)
        set(gca,'xtick',([pi/4 pi/2 .75*pi pi]))
        set(gca,'xticklabels',({'45','90','135','180'}),'FontSize',20)
        set(gca,'ytick',([0 pi/4 pi/2 .75*pi pi]))
        set(gca,'yticklabels',({'0','45','90','135','180'}),'FontSize',20)
        n_cells(t) = sum(dataAll)-(sum(isnan(pre)) + sum(isnan(post)))
        text(.2, 3, ['n = ' num2str(n_cells(t))],'FontSize',20)
        text(.2, 2.6, ['cc = ' num2str(c,'%.2f')],'FontSize',18)
    end
end


%% pref theta for mean of stat/running

data = goodAll & useResp' & ~inhAll & hasDrift==1 & useOsi==1 ;%which cells to include
dataInh =goodAll & useResp' & inhAll & hasDrift==1 & useOsi==1;
dataAll =goodAll & useResp' & hasDrift==1 & useOsi==1;

prefAlltheta = squeeze(mean(drift_pref_theta,2))
titles={'Saline', 'DOI'};
figure
for t=1:2
    pre=prefAlltheta(dataAll&treatment==t,1); post= prefAlltheta(dataAll&treatment==t,2)
    post(pre-post > pi/2) = post(pre-post>pi/2)+pi; pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;
    c=corrcoef(pre,post,'rows','pairwise'); c = c(2,1)
    subplot(1,2,t)
        hold on; axis square;plot([0 3.25], [0 3.25]);xlim([0 3.25]);ylim([0 3.25]);
    pre=prefAlltheta(dataInh&treatment==t,1); post= prefAlltheta(dataInh&treatment==t,2)
    post(pre-post > pi/2) = post(pre-post>pi/2)+pi;
    pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;
    plot(pre,post,'r.','Markersize',18);title(titles{t});
    pre=prefAlltheta(data&treatment==t,1); post= prefAlltheta(data&treatment==t,2)
    post(pre-post > pi/2) = post(pre-post>pi/2)+pi;
    pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;
    plot(pre,post,'b.','Markersize',18)
    n_cells(t) = sum(dataAll&treatment==t)-(sum(isnan(pre)) + sum(isnan(post)))


    text(.2, 2.8, ['r^2 = ' num2str(rs_pf(t),'%.2f')],'FontSize',18)
    text(.2, 3, ['n = ' num2str(n_cells(t))],'FontSize',20)
    text(.2, 2.6, ['cc = ' num2str(c,'%.2f')],'FontSize',18)
    
    %set(gca,'xtick',0:pi/2:pi,'xticklabel',1:180,'Fontsize',14)
%     set(gca,'xtick',([pi/4 pi/2 .75*pi pi]))
    % set(gca,'xticklabels',({'\pi/4','\pi/2','3/4\pi','\pi'}),'FontSize',20)
        %set(gca,'yticklabels',({'0','\pi/4','\pi/2','3/4\pi','\pi'}),'FontSize',20)

     set(gca,'xtick',([pi/4 pi/2 .75*pi pi]))
    set(gca,'xticklabels',({'45','90','135','180'}),'FontSize',20)
    set(gca,'ytick',([0 pi/4 pi/2 .75*pi pi]))
    set(gca,'yticklabels',({'0','45','90','135','180'}),'FontSize',20)
end

%%
figure
for t=1:2
    useN=dataAll & treatment==t & hasDrift==1
    miDrift= (drift_pref_theta(useN,1,2))-(drift_pref_theta(useN,1,1))...
        ./drift_pref_theta(useN,1,2)+ (drift_pref_theta(useN,1,1));
    h= hist(miDrift,-1:.1:1);
    Mbins=-1:.1:1
    % hist(drift_pref_theta(dataAll & treatment==2,1,1))
    subplot(1,2,t)
    bar(Mbins,h/sum(useN),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .7]); xlim([-1.5 1.5]);axis square; title('layer 2/3')
    
end

%% CV OSI
%%% for paper
drift_osiCV_abs=abs(drift_osiCV);


useOsi_pre = drift_osiCV_abs(:,1,1)>.3;
useOsi_post = drift_osiCV_abs(:,1,2)>.3;
useOsi=(useOsi_pre & useOsi_post)'
useOsi_mn = squeeze(mean(drift_osiCV_abs,2));
peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont

data = goodAll & useResp' & ~inhAll & hasDrift==1 &useOsi==1 ;%which cells to include
dataInh =goodAll & useResp' & inhAll & hasDrift==1 & useOsi==1;
dataAll =goodAll & useResp' & hasDrift==1 &useOsi==1;
titles={'Saline', 'DOI'};
figure; hold on;
for mv= 1
    clear R rs_pf
   % figure
    if mv==1, set(gcf,'Name', 'stationary');
        useResp = peakresp(:,1,1)>1 | peakresp(:,1,2)>1; set(gcf,'Name','stationary pref ori');
        
    else
        set(gcf,'Name', 'mv');
        useResp = peakresp(:,2,1)>1 | peakresp(:,2,2)>1; set(gcf,'Name','mv pre ori')
        
    end
    for t=1
        data = goodAll & useResp' & ~inhAll & hasDrift==1 &useOsi==1&treatment==t ;%which cells to include
        dataInh =goodAll & useResp' & inhAll & hasDrift==1 & useOsi==1&treatment==t;
        dataAll =goodAll & useResp' & hasDrift==1 &useOsi==1&treatment==t;
        pre=drift_pref_thetaCV(dataAll,mv,1); post= drift_pref_thetaCV(dataAll,mv,2)
        post(pre-post > pi/2) = post(pre-post>pi/2)+pi; pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;
        
        mdl_pf= fitlm(pre,post)
        c=corrcoef(pre,post,'rows','pairwise'); c = c(2,1)
        rs_pf(t) = mdl_pf.Rsquared.Ordinary
        %subplot(1,2,t)
%          pre=drift_pref_thetaCV(dataInh,mv,1); post= drift_pref_thetaCV(dataInh,mv,2)
%         post(pre-post > pi/2) = post(pre-post>pi/2)+pi;
%         pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;
       
        plot(pre,post,'bo','Markersize',10,'linewidth',2)%,'MarkerFaceColor',[.9 .2 .5]);
        title(titles{t});hold on;
        pre=drift_pref_thetaCV(data,mv,1); post= drift_pref_thetaCV(data,mv,2)
        post(pre-post > pi/2) = post(pre-post>pi/2)+pi;
        pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;
        
        plot(pre,post,'bo','Markersize',10,'linewidth',2)%,'MarkerFaceColor',[.2 0.5 .8])
        [h stat(t)] = ttest(pre,post);
        axis square;plot([0 3.25], [0 3.25]);xlim([0 3.25]);ylim([0 3.25]);
        text(.2, 2.8, ['r^2 = ' num2str(rs_pf(t),'%.2f')],'FontSize',18)
        text(.2, 2.4, ['p = ' num2str(stat(t))],'FontSize',18)
        
       
        
        set(gca,'xtick',([pi/4 pi/2 .75*pi pi]))
        set(gca,'xticklabels',({'45','90','135','180'}),'FontSize',20)
        set(gca,'ytick',([0 pi/4 pi/2 .75*pi pi]))
        set(gca,'yticklabels',({'0','45','90','135','180'}),'FontSize',20)
        n_cells(t) = sum(dataAll)-(sum(isnan(pre)) + sum(isnan(post)))
        text(.2, 3, ['n = ' num2str(n_cells(t))],'FontSize',20)
        text(.2, 2.6, ['cc = ' num2str(c,'%.2f')],'FontSize',18)
    end
end

%% normal osi
useOsi_pre = drift_osi(:,1,1)>.3;
useOsi_post = drift_osi(:,1,2)>.3;
useOsi=(useOsi_pre & useOsi_post)'
useOsi_mn = squeeze(mean(drift_osi,2));
peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont
% peakresp = squeeze(mean(drift_orient,2))-squeeze(drift_spont);%subtract driftspont
useResp = peakresp(:,2,1)>2| peakresp(:,2,2)>2;
ampmean = squeeze(nanmean(peakresp,2));
% useResp = ampmean(:,1)>2 |ampmean(:,2)>2;

data = goodAll & useResp' & ~inhAll & hasDrift==1 &useOsi==1 ;%which cells to include
dataInh =goodAll & useResp' & inhAll & hasDrift==1 & useOsi==1;
dataAll =goodAll & useResp' & hasDrift==1 &useOsi==1;
titles={'Saline', 'DOI'};
for mv= 1:2
    clear R rs_pf
    figure
    if mv==1, set(gcf,'Name', 'stationary');
        useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2; set(gcf,'Name','stationary pref ori');
        
    else
        set(gcf,'Name', 'mv');
        useResp = peakresp(:,2,1)>2 | peakresp(:,2,2)>2; set(gcf,'Name','mv pre ori')
        
    end
    for t=1:2
        data = goodAll & useResp' & ~inhAll & hasDrift==1 &useOsi==1&treatment==t ;%which cells to include
        dataInh =goodAll & useResp' & inhAll & hasDrift==1 & useOsi==1&treatment==t;
        dataAll =goodAll & useResp' & hasDrift==1 &useOsi==1&treatment==t;
        pre=drift_pref_theta(dataAll,mv,1); post= drift_pref_theta(dataAll,mv,2)
        post(pre-post > pi/2) = post(pre-post>pi/2)+pi; pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;
        
        mdl_pf= fitlm(pre,post)
        c=corrcoef(pre,post,'rows','pairwise'); c = c(2,1)
        rs_pf(t) = mdl_pf.Rsquared.Ordinary
        subplot(1,2,t)
        pre=drift_pref_theta(dataInh,mv,1); post= drift_pref_theta(dataInh,mv,2)
        post(pre-post > pi/2) = post(pre-post>pi/2)+pi;
        pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;
       
        plot(pre,post,'ko','Markersize',10,'linewidth',2,'MarkerFaceColor',[.9 .2 .5]);title(titles{t});hold on;
        pre=drift_pref_theta(data,mv,1); post= drift_pref_theta(data,mv,2)
        post(pre-post > pi/2) = post(pre-post>pi/2)+pi;
        pre(post-pre > pi/2) = pre(post-pre>pi/2)+pi;
        
        plot(pre,post,'ko','Markersize',10,'linewidth',2,'MarkerFaceColor',[.2 0.5 .8])
        axis square;plot([0 3.25], [0 3.25]);xlim([0 3.25]);ylim([0 3.25]);
        text(.2, 2.8, ['r^2 = ' num2str(rs_pf(t),'%.2f')],'FontSize',18)
        set(gca,'xtick',([pi/4 pi/2 .75*pi pi]))
        set(gca,'xticklabels',({'45','90','135','180'}),'FontSize',20)
        set(gca,'ytick',([0 pi/4 pi/2 .75*pi pi]))
        set(gca,'yticklabels',({'0','45','90','135','180'}),'FontSize',20)
        n_cells(t) = sum(dataAll)-(sum(isnan(pre)) + sum(isnan(post)))
        text(.2, 3, ['n = ' num2str(n_cells(t))],'FontSize',20)
        text(.2, 2.6, ['cc = ' num2str(c,'%.2f')],'FontSize',18)
    end
end
%%
peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont

useOsi_pre = drift_osi(:,2,1)>0
useOsi_post = drift_osi(:,2,2)>0
useOsi=(useOsi_pre | useOsi_post)'
%useResp = peakresp(:,2,1)>2 & peakresp(:,2,2)>2;
% data = goodAll & useResp' &  useOsi==1 ;%which cells to include
dataAll = goodAll & useResp' & hasDrift==1 & useOsi ==1 ;

for mv=1:2
    figure
    if mv ==1
        useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2; set(gcf,'Name','stationary OSI');
    else
        useResp = peakresp(:,2,1)>2 | peakresp(:,2,2)>2;
        set(gcf,'Name','mv OSI'), end
    for t=1:2
        dataAll = goodAll & useResp' & hasDrift==1 &treatment==t;
        subplot(1,2,t)
        pre = drift_osi(dataAll&treatment==t,mv,1);
        post = drift_osi(dataAll&treatment==t,mv,2);
        mdl_osi= fitlm(pre,post)
        rs_osi(t) = mdl_osi.Rsquared.Ordinary
        cc = corrcoef(pre,post,'rows','pairwise'); cc = cc(2,1);
        plot(drift_osi(dataAll&~inhAll,mv,1),drift_osi(dataAll& ~inhAll,mv,2),'.','Markersize',20);hold on;
        plot(drift_osi(dataAll & inhAll==1,mv,1),drift_osi(dataAll& inhAll==1,mv,2),'r.','Markersize',20);
        axis square;hold on; plot([0 1],[0 1]);
        n_cells(t) = sum(dataAll)-(sum(isnan(pre)) + sum(isnan(post)))
        % text(.1, .8, ['r^2 = ' num2str(rs_osi(t),'%.2f')],'FontSize',18)
        text(.1, .8, ['c.c = ' num2str(cc,'%.2f')],'FontSize',20)
        text(.1, .9, ['n = ' num2str(n_cells(t))],'FontSize',20)
        
        title(titles{t},'FontSize', 22)
    end
end

%% does OSI or preferred theta change as a function of change in FR?
useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2;
data = goodAll & useResp' & ~inhAll & hasDrift==1 ;%which cells to include
dataInh =goodAll & useResp' & inhAll & hasDrift==1 ;
dataAll =goodAll & useResp' & hasDrift==1 ;

mv=1
for t=1:2
    figure
    pre = drift_osi(dataAll&treatment==t,mv,1);
    post = drift_osi(dataAll&treatment==t,mv,2);
    miDrift = (mean(peakresp(dataAll&treatment==t,1,2),3)-mean(peakresp(dataAll&treatment==t,mv,1),3))./...
        (mean(peakresp(dataAll&treatment==t,mv,2),3)+mean(peakresp(dataAll&treatment==t,mv,1),3));
    h = hist(miDrift,-1:.2:1);
    
    subplot(121)
    plot(miDrift,pre,'.','Markersize',20); xlabel('MI'); ylabel('OSI Pre'); axis square
    xlim([-1.5 1.5]);
    subplot(122);
    plot(miDrift,post,'.','Markersize',20); xlabel('MI'); ylabel('OSI Post');axis square
    xlim([-1.5 1.5]);
    
end
%%

mv=1
for t=1:2
    figure
    pre = drift_pref_theta(dataAll&treatment==t,mv,1);
    post = drift_pref_theta(dataAll&treatment==t,mv,2);
    miDrift = (mean(peakresp(dataAll&treatment==t,mv,2),3)-mean(peakresp(dataAll&treatment==t,mv,1),3))./...
        (mean(peakresp(dataAll&treatment==t,mv,2),3)+mean(peakresp(dataAll&treatment==t,mv,1),3));
    h = hist(miDrift,-1:.2:1);
    
    subplot(121)
    plot(miDrift,pre,'.','Markersize',20); xlabel('MI'); ylabel('OSI Pre'); axis square
    xlim([-1.5 1.5]);
    subplot(122);
    plot(miDrift,post,'.','Markersize',20); xlabel('MI'); ylabel('OSI Post');axis square
    xlim([-1.5 1.5]);
end

%% dsi pref theta
clear peakresp useResp
% useDsi_pre = drift_dsi(:,1,1)>.3;
% useDsi_post = drift_dsi(:,1,2)>.3;
% useDsi=(useDsi_pre & useDsi_post)'

useDsi_mn = squeeze(mean(drift_dsi,2));
peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont
% peakresp = squeeze(mean(drift_orient,2))-squeeze(drift_spont);%subtract driftspont
% useResp = peakresp(:,2,1)>2| peakresp(:,2,2)>2;
%allOrient = squeeze(mean(drift_orient,2));


%ampmean = squeeze(nanmean(peakresp,2));
% useResp = ampmean(:,1)>2 |ampmean(:,2)>2;


titles={'Saline', 'DOI'};
for mv= 1
    clear R rs_pf pre post useResp useDsi
    figure
    if mv==1, set(gcf,'Name', 'stationary');
        useResp = peakresp(:,1,1)>2| peakresp(:,1,2)>2;
        useDsi_pre = drift_dsi(:,1,1)>.3; useDsi_post = drift_dsi(:,1,2)>.3;
        useDsi=(useDsi_pre | useDsi_post)'
        
    else
        set(gcf,'Name', 'mv');
        useResp = peakresp(:,2,1)>2| peakresp(:,2,2)>2;
        useDsi_pre = drift_dsi(:,2,1)>.3; useDsi_post = drift_dsi(:,2,2)>.3;
        useDsi=(useDsi_pre | useDsi_post)'
    end
    
    for t=1:2
        data = goodAll  & inhAll==0 & hasDrift==1 &useDsi==1 &treatment==t &useResp';
        dataInh =goodAll  & inhAll==1 & hasDrift==1 & useDsi==1&treatment==t &useResp';
        dataAll =goodAll & hasDrift==1 &useDsi==1&treatment==t &useResp';
        pre=drift_dsi_theta(dataAll,mv,1);
        post= drift_dsi_theta(dataAll,mv,2);
        post(pre-post > pi) = post(pre-post>pi)+2*pi; pre(post-pre > pi) = pre(post-pre>pi)+2*pi;
        mdl_pf= fitlm(pre,post)
        c=corrcoef(pre,post,'rows','pairwise'); c = c(2,1)
        rs_pf(t) = mdl_pf.Rsquared.Ordinary
        subplot(1,2,t)
        pre=drift_dsi_theta(data&treatment==t,mv,1);
        post= drift_dsi_theta(data&treatment==t,mv,2)
        post(pre-post > pi) = post(pre-post>pi)+2*pi;
        pre(post-pre > pi) = pre(post-pre>pi)+2*pi;
        plot(pre,post,'.','Markersize',18)
        hold on; axis square;plot([0 8.6], [0 8.6]); xlim([0 8.6]);ylim([0 8.6]);
        [h p(t)] = ttest(pre,post)
        
        preInh=drift_dsi_theta(dataInh,mv,1);
        postInh= drift_dsi_theta(dataInh,mv,2)
        postInh(preInh-postInh > pi) = postInh(preInh-postInh>pi)+2*pi;
        preInh(postInh-preInh > pi) = preInh(postInh-preInh>pi)+2*pi;
        plot(preInh,postInh,'r.','Markersize',18);title(titles{t});
        n_cells(t) = sum(dataAll)-(sum(isnan(pre)) + sum(isnan(post)))

        text(.2, 7.4, ['r^2 = ' num2str(rs_pf(t),'%.2f')],'FontSize',18)
        text(.2, 7.8, ['n = ' num2str(n_cells(t))],'FontSize',20)
        text(.2, 7.0, ['cc = ' num2str(c,'%.2f')],'FontSize',18)
        text(.2, 6.8, ['p = ' num2str(p(t))],'FontSize',18)
        
        %set(gca,'xtick',0:pi/2:pi,'xticklabel',1:180,'Fontsize',14)
        %         set(gca,'xtick',([pi/4 pi/2 .75*pi pi]))
        %     set(gca,'xticklabels',({'\pi/4','\pi/2','3/4\pi','\pi'}),'FontSize',20)
        %         set(gca,'xticklabels',({'45','90','135','180'}),'FontSize',20)
        %         set(gca,'ytick',([0 pi/4 pi/2 .75*pi pi]))
        % set(gca,'yticklabels',({'0','\pi/4','\pi/2','3/4\pi','\pi'}),'FontSize',20)
        %         set(gca,'yticklabels',({'0','45','90','135','180'}),'FontSize',20)
    end
end


%%
clear useDsi useResp dataAll rs_dsi c
peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont); %subtract driftspont
%peakresp = squeeze(mean(drift_orient,2))-squeeze(drift_spont);%subtract driftspont

% % data = goodAll & useResp' &  usedsi==1 ;%which cells to include
% dataAll = goodAll & useResp' & hasDrift==1 & usedsi ==1 ;

for mv=2
    figure
    if mv ==1
        useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2; set(gcf,'Name','stationary dsi');
        useDsi_pre = drift_dsi(:,1,1)>0; useDsi_post = drift_dsi(:,1,2)>0; useDsi=(useDsi_pre | useDsi_post)'
    else
        useResp = peakresp(:,2,1)>2 | peakresp(:,2,2)>2; set(gcf,'Name','mv dsi')
        useDsi_pre = drift_dsi(:,2,1)>0; useDsi_post = drift_dsi(:,2,2)>0; useDsi=(useDsi_pre | useDsi_post)'
    end
    for t=1:2
        dataAll = goodAll & useResp' & hasDrift==1 &treatment==t & useDsi==1;
        subplot(1,2,t)
        mdl_dsi= fitlm(drift_dsi(dataAll,mv,1),drift_dsi(dataAll,mv,2))
        rs_dsi(t) = mdl_dsi.Rsquared.Ordinary
        c = corrcoef(drift_dsi(dataAll,mv,1),drift_dsi(dataAll,mv,2),'rows','pairwise');
        c = c(2,1);
       [p(t) h]=ranksum(drift_dsi(dataAll&treatment==t,mv,1),drift_dsi(dataAll&treatment==t,mv,2))

        plot(drift_dsi(dataAll & ~inhAll,mv,1),drift_dsi(dataAll & ~inhAll,mv,2),'.','Markersize',20);hold on;
        plot(drift_dsi(dataAll & inhAll==1,mv,1),drift_dsi(dataAll & inhAll==1,mv,2),'r.','Markersize',20);
        axis square;hold on; plot([0 1],[0 1]);
        n_cells(t) = sum(treatment==t & dataAll)
        text(.02, .85, ['r^2 = ' num2str(rs_dsi(t),'%.2f')],'FontSize',18)
        text(.02 ,.75, ['cc = ' num2str(c,'%.2f')],'FontSize',18)
        text(.02 ,.95, ['n = ' num2str(n_cells(t))],'FontSize',18)
        title(titles{t})
    end
end

%% DSI & pref theta w/MI
useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2;
data = goodAll & useResp' & ~inhAll & hasDrift==1 ;%which cells to include
dataInh =goodAll & useResp' & inhAll & hasDrift==1 ;
dataAll =goodAll & useResp' & hasDrift==1 ;
%
% mv=1
% for t=1:2
%     figure
% pre = drift_dsi(dataAll&treatment==t,mv,1);
% post = drift_dsi(dataAll&treatment==t,mv,2);
% 
% miDrift = (mean(peakresp(dataAll&treatment==t,1,2),3)-mean(peakresp(dataAll&treatment==t,mv,1),3))./...
%     (mean(peakresp(dataAll&treatment==t,mv,2),3)+mean(peakresp(dataAll&treatment==t,mv,1),3));
% h = hist(miDrift,-1:.2:1);
% 
% subplot(121)
% plot(miDrift,pre,'.','Markersize',20); xlabel('MI'); ylabel('DSI Pre'); axis square
% xlim([-1.5 1.5]);
% subplot(122);
% plot(miDrift,post,'.','Markersize',20); xlabel('MI'); ylabel('DSI Post');axis square
% xlim([-1.5 1.5]);
% 
% end

mv=1
for t=1:2
    figure
    pre = drift_dsi_theta(dataAll&treatment==t,mv,1);
    post = drift_dsi_theta(dataAll&treatment==t,mv,2);
    
    miDrift = (mean(peakresp(dataAll&treatment==t,1,2),3)-mean(peakresp(dataAll&treatment==t,mv,1),3))./...
        (mean(peakresp(dataAll&treatment==t,mv,2),3)+mean(peakresp(dataAll&treatment==t,mv,1),3));
    h = hist(miDrift,-1:.2:1);
    
    subplot(121)
    plot(miDrift,pre,'.','Markersize',20); axis square; xlabel('MI'); ylabel('DSI Theta Pre');
 %   xlim([-1.5 1.5]);
    subplot(122);
    plot(miDrift,post,'.','Markersize',20);axis square; xlabel('MI'); ylabel('DSI Theta Post');
   % xlim([-1.5 1.5]);
    
    
end

peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont

useOsi_pre = drift_osi(:,2,1)>0
useOsi_post = drift_osi(:,2,2)>0
useOsi=(useOsi_pre | useOsi_post)'
useResp = peakresp(:,2,1)>2 & peakresp(:,2,2)>2;

% data = goodAll & useResp' &  useOsi==1 ;%which cells to include
dataAll = goodAll & useResp' & hasDrift==1 & useOsi ==1 ;

titles = {'Saline','DOI','5HT'};

%%
for mv=1:2
    figure
    if mv ==1
        useResp = peakresp(:,1,1)>0 & peakresp(:,1,2)>0;
        set(gcf,'Name','stationary mean OSI')
    else
        useResp = peakresp(:,2,1)>0 & peakresp(:,2,2)>0;
        set(gcf,'Name','mv mean OSI');
    end
    
    %dataAll = goodAll & useResp' & hasDrift==1 & useOsi ==1 ;
    %figure
    for t=1:2
        clear p
        osi(1,:) = nanmean(abs(drift_osi(dataAll&treatment==t,mv,1)))
        osierr(1,:) = squeeze(nanstd(drift_osi(dataAll & treatment==t,mv,1)/sqrt(sum(drift_osi(dataAll & treatment==t,mv,1),1))));
        osi(2,:) = nanmean(abs(drift_osi(dataAll&treatment==t,mv,2)))
        osierr(2,:) = squeeze(nanstd(drift_osi(dataAll & treatment==t,mv,2)/sqrt(sum(drift_osi(dataAll & treatment==t,mv,2),1))));
        subplot(1,2,t)
        barweb(osi,osierr); ylim([0 1]); axis square
        [p h]=ranksum(drift_osi(dataAll&treatment==t,mv,1),drift_osi(dataAll&treatment==t,mv,2))
        text(.55, .9, ['p = ' num2str(p,'%.2f')],'FontSize',18)
        n_cells(t) = sum(dataAll&treatment==t)
        text(.55, .8, ['n = ' num2str(n_cells(t))],'FontSize',18)
        title(titles{t})
        set(gca,'xtick',[1 2])
        set(gca,'xticklabel',['Pre     Post'])
        set(gca,'xticklabel',['Pre      Post'])
    end
end

%%
%osi, mean of running/stationary

drift_osi_mn=squeeze(mean(drift_osi,2));

useOsi_pre = drift_osi_mn(:,1)>0
useOsi_post = drift_osi_mn(:,2)>0
useOsi=(useOsi_pre | useOsi_post)'
peakresp_mn = squeeze(nanmean(peakresp,2));
useResp = peakresp_mn(:,1)>1 | peakresp_mn(:,2)>1;
% data = goodAll & useResp' &  useOsi==1 ;%which cells to include
dataAll = goodAll & useResp' & hasDrift==1 %& useOsi ==1 ;

titles = {'Saline','DOI','5HT'};
%osi all layers
clear n_cells p
clear osierr osi
figure
for t=1:2
    % dataAll = goodAll & useResp' & hasDrift==1 & (useOsi ==1)' &treatment==t;
    osi(1,:) = nanmean(abs(drift_osi_mn(dataAll&treatment==t,1)))
    osierr(1,:) = squeeze(nanstd(drift_osi_mn(dataAll&treatment==t,1))./sqrt(sum(dataAll&treatment==t)));
    osi(2,:) = nanmean(abs(drift_osi_mn(dataAll&treatment==t,2)))
    osierr(2,:) = squeeze(nanstd(drift_osi_mn(dataAll&treatment==t,2))./sqrt(sum(dataAll&treatment==t)));
    subplot(1,2,t)
    barweb(osi,osierr); ylim([0 1]); axis square
    [p,h]=ranksum(drift_osi_mn(dataAll&treatment==t,1),drift_osi_mn(dataAll&treatment==t,2))
   text(.55, .9, ['p = ' num2str(p)],'FontSize',18)
    n_cells(t) = sum(dataAll&treatment==t)
    text(.55, .8, ['n = ' num2str(n_cells(t))],'FontSize',18)
    title(titles{t})
    set(gca,'xtick',[1 2])
    set(gca,'xticklabel',['Pre     Post'])
end

%%
clear dataAll
useResp = peakresp(:,1,1)>2 | peakresp(:,1,2)>2;
dataAll = goodAll & useResp' & hasDrift==1 & useOsi ==1 &inhAll==0;

%osi specific layers

for t=2
    figure
    if t == 1, set(gcf,'Name','Saline')
    elseif t==2, set(gcf,'Name','DOI')
    else set(gcf,'Name','5HT'), end
    
    osi(1,:) = nanmean(abs(drift_osi(dataAll&treatment==t&(layerAll==2|layerAll==3),1,1)))
    osierr(1,:) = squeeze(nanstd(drift_osi(dataAll & treatment==t&(layerAll==2|layerAll==3),1,1)...
        /sqrt(sum(drift_osi(dataAll & treatment==t&(layerAll==2|layerAll==3),1,1),1))));
    osi(2,:) = nanmean(abs(drift_osi(dataAll&treatment==t&(layerAll==2|layerAll==3),1,2)))
    osierr(2,:) = (squeeze(nanstd(drift_osi(dataAll & treatment==t&(layerAll==2|layerAll==3),1,2))/...
        sqrt(sum(drift_osi(dataAll & treatment==t&(layerAll==2|layerAll==3),1,2),1))));
    subplot(1,3,1)
    barweb(osi,osierr); ylim([0 1]); axis square
    [p h]=ranksum(drift_osi(dataAll&treatment==t&(layerAll==2|layerAll==3),1,1)...
        ,drift_osi(dataAll&treatment==t&(layerAll==2|layerAll==3),1,2))
    text(.9, 1, ['p = ' num2str(p,'%.2f')],'FontSize',18)
    for l=4:5
        clear osierr
        use = dataAll & treatment==t & layerAll==l;
        osi(1,:) = nanmean(abs(drift_osi(use,1,1)))
        osierr(1,:) = squeeze(nanstd(drift_osi(use,1,1))) / (sqrt(sum(use)));
        %  sqrt(sum(drift_osi(use,1,1),1)));
        osi(2,:) = nanmean(abs(drift_osi(use,1,2)))
        osierr(2,:) = squeeze(nanstd(drift_osi(use,1,2))) / (sqrt(sum(use)));
        %  sqrt(sum(drift_osi(use,1,2),1)));
        
        subplot(1,3,l-2)
        barweb(osi,osierr); ylim([0 1]); axis square
        [p h]=ranksum(drift_osi(dataAll&treatment==t&layerAll==l,1,1),...
            drift_osi(dataAll&treatment==t&layerAll==l,1,2))
        text(.9, 1, ['p = ' num2str(p,'%.2f')],'FontSize',18)
    end
end
%%
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

%%
figure
for t=1:2
%clear allSF 

for g = 1:7
    data = useResp' & hasDrift==1 & treatment==t
    allSF(g,:,t)= sum(SFpref(data,:)==g)/sum(data)
subplot(1,2,t)    
plot(allSF(:,:,t)); axis square; hold on;plot(allSF(:,:,t),'.','Markersize',15);
ylim([.05 .24]);
hold off
% [h(t) p(t)] = kruskalwallis(allSF(:,1,t),allSF(:,2,t))
end

end

%[h p] = kstest2(allSF(:,1),allSF(:,2))
%%
figure
for t= 1:2
    sess=find(sessionTreatment==t)

for g = 1:length(sess)

data = useResp' & hasDrift==1 & treatment==t
aniCell=find(data&sessionNum==sess(g))
indi=SFpref(aniCell,:);

for sf=1:7
    prop(sf,:,:,g,t) = sum(indi==sf)./length(indi)
end

end
 subplot(1,2,t)
 plot(squeeze(prop(:,1,1,:,t)));axis square
end

figure
for t=1:2
h=squeeze(nanmean(prop,4)) %(sf,prepost,treatment)
err=nanstd(prop,[],4)/sqrt(length(prop));err=squeeze(err);
%p(t) = kruskalwallis(h(:,:,t))
subplot(1,2,t)
%plot(h(:,:,t));
errorbar(h(:,:,t),err(:,:,t));axis square; ylim([0 1]);
% hold on;plot([0 8],[1/7 1/7],':','color',[0.8 0.8 0.8])
end


%%
peakresp = squeeze(max(drift_orient,[],2))-squeeze(drift_spont);%subtract driftspont
useResp = peakresp(:,1,1)>0 & peakresp(:,1,2)>0;

for t = 1:2
    figure
    clear data barsf mn  prepost ff low medium high
        sess=find(sessionTreatment==t);
        data = useResp' & hasDrift==1 & treatment==t;
    for s=1:length(sess)
clear aniCell 
for prepost=1:2
    %for g = 1:4
         aniCell=(data&sessionNum==sess(s));
         ff(s) = (sum(SFpref(aniCell,prepost)==1))./length(aniCell);
         low(s) = ((SFpref(aniCell,prepost)==2|SFpref(aniCell,prepost)==3))/length(aniCell); 
         medium(s) = ((SFpref(aniCell,prepost)==4|SFpref(aniCell,prepost)==5))/length(aniCell); 
         high(s) = ((SFpref(aniCell,prepost)==6|SFpref(aniCell,prepost)==7))/length(aniCell);
       

%          ff(s) = (sum(SFpref(aniCell,prepost)==1))/sum(aniCell);
%          low(s) = (sum(SFpref(aniCell,prepost)==2|SFpref(aniCell,prepost)==3))/sum(aniCell); 
%          medium(s) = (sum(SFpref(aniCell,prepost)==4|SFpref(aniCell,prepost)==5))/sum(aniCell); 
%          high(s) = (sum(SFpref(aniCell,prepost)==6|SFpref(aniCell,prepost)==7))/sum(aniCell);

          barsf(s,:)= [ff(s) low(s) medium(s) high(s)];
          bartest(s,:,prepost) =  [ff(s) low(s) medium(s) high(s)];
          mn(:,prepost) = nanmean(barsf,1); 
          err(:,prepost) = (nanstd(barsf,1)/sqrt(length(barsf)));  
   % end
            %mn(prepost,:) = nanmean(barsf,1); err(prepost,:) = (nanstd(barsf)/sqrt(length(sess)));  


%         barweb(mn,err);% hold on;
        %plot([0 8],[1/4 1/4],':','color',[0.8 0.8 0.8])

         ylim([0 .5]); axis square
         xlim([.5 4.5]); ylabel('proportion of cells'); xlabel('preferred SF');
         set(gca, 'XTickLabels', {'FF' 'low' 'medium' 'high'},'FontSize',15);
        %text(.9, 1, ['p = ' num2str(p(g),'%.2f')],'FontSize',18)
        
%         if t==1, title('Saline');else title('DOI')
%         end
    end
          % [h p]= ttest(bartest(:,:,1),bartest(:,:,2))
end
    end


%%
SFpref =[prefSF ;prefSFPost]';
ff = SFpref==1
low = SFpref==2|SFpref==3; medium = SFpref==4|SFpref==5; high = SFpref==6|SFpref==7;
sfs= {ff low medium high};
warning off

figure
for t = 1:2
    clear data sfBar
    for g = 1:length(sfs)
        data = useResp' & hasDrift==1 & treatment==t% & (sfs{1},1)';
        sfBar(:,1) = hist(SFpref(data,1),[1 2:2:7]);sfBar(:,2) = hist(SFpref(data,2),[1 2:2:7]);
        
        SFpref(data,1)==sfs{g}(data,1)
        
     %  sferr(g,1) = squeeze(nanstd(SFpref(data,1)==(sfs{g},1)))/sqrt(sum(sfs{g}(data,1)))
       % sferr(g,2) = squeeze(nanstd(SFpref(data,2)==sfs{g}(data,2)))/sqrt(sum(sfs{g}(data,2)))
        
     %   sferr(g,2) = squeeze(nanstd(sfBar(g,2)))/sqrt(sum(sfBar(g,2)))
        
     %   sferr(g,1) = squeeze(nanstd(SFpref(sfs{g}(data,1))))/sqrt(sum(SFpref(sfs{g}(data,1))))
      %  sferr(g,2) = squeeze(nanstd(SFpref(sfs{g}(data,2))))/sqrt(sum(SFpref(sfs{g}(data,2))));
        
        %[p h] = ranksum(SFpref(data,1),SFpref(data,2))
        [p(g) h(g)] = ranksum(sfs{g}(data,1),sfs{g}(data,2))
        subplot(2,1,t)
      %  bar(sfBar/sum(data))
    %    barweb(sfBar/sum(data),sferr);
        ylim([0 .5]); axis square
        xlim([.5 4.5]); ylabel('proportion of cells'); xlabel('preferred SF');
        set(gca, 'XTickLabels', {'FF' 'low' 'medium' 'high'},'FontSize',22);
        %text(.9, 1, ['p = ' num2str(p(g),'%.2f')],'FontSize',18)
        
        if t==1, title('Saline');else title('DOI')
        end
    end
end



%%
%response to each SF prepost
for mv=1:2
    peakresp = squeeze(max(drift_orient,[],2))%-squeeze(drift_spont);%subtract driftspont
    useResp = peakresp(:,mv,1)>2 | peakresp(:,mv,2)>2;
    data = useResp' & hasDrift==1
    for t=1:2
        figure
        if t==1, set(gcf,'Name','Saline'); else set(gcf,'Name','DOI')
        end
        for sf = 2:7
            subplot(2,3,sf-1)
            plot(drift_sf(data&treatment==t,sf,1,1),drift_sf(data&treatment==t,sf,1,2),'o')
            axis square; hold on; plot([0 15],[0 15])
            
        end
    end
end
%%
data = goodAll & useResp' & hasDrift==1 &inhAll==0;
dataInh = goodAll & useResp' & hasDrift==1 &inhAll==1;

titlesInh = {'Saline Inh','DOI Inh','5HT Inh'};
for mv=1:2
    figure
    if mv==1, set(gcf,'Name', 'mean SF stationary')
    else set(gcf,'Name', 'mean SF mv'), end
    for t=1:2
        SFpre= nanmean(drift_sf(data & treatment==t & hasDrift==1,2:7,mv,1),1)
        stderrorPre = nanstd(drift_sf(data & treatment==t & hasDrift==1,2:7,mv,1)) / sqrt(sum(data & treatment==t & hasDrift==1));
        SFpost = nanmean(drift_sf(data & treatment==t & hasDrift==1,2:7,mv,2),1)
        stderrorPost = nanstd(drift_sf(data & treatment==t & hasDrift==1,2:7,mv,2)) / sqrt(sum(data & treatment==t & hasDrift==1));
        subplot(2,2,t)
        errorbar(SFpre,stderrorPre) %'-^','MarkerEdgeColor','r','Color','r','LineWidth', 2);
        hold on;
        errorbar(SFpost,stderrorPost)
        %plot(SFpost,'r','LineWidth',2);axis square
        title(titles{t}); xlabel('SF'); ylabel('sp/sec');xlim([1 6]); set(gca, 'XTickLabels', [2 3 4 5 6 7]);
        SFpreInh= nanmean(drift_sf(dataInh & treatment==t & hasDrift==1,2:7,mv,1),1)
        stderrorPreInh = nanstd(drift_sf(dataInh & treatment==t & hasDrift==1,2:7,mv,1)) / sqrt(sum(dataInh & treatment==t & hasDrift==1));
        
        SFpostInh = nanmean(drift_sf(dataInh & treatment==t & hasDrift==1,2:7,mv,2),1)
        stderrorPostInh = nanstd(drift_sf(dataInh & treatment==t & hasDrift==1,2:7,mv,2)) / sqrt(sum(dataInh & treatment==t & hasDrift==1));
        
        subplot(2,2,t+2)
        errorbar(SFpreInh,stderrorPreInh)%'LineWidth',2);
        hold on; errorbar(SFpostInh,stderrorPostInh)%'r','LineWidth',2);
        axis square; xlim([1 6]); set(gca, 'XTickLabels', [2 3 4 5 6 7]);
        title(titlesInh{t}); xlabel('SF'); ylabel('sp/sec');
    end
end

%%

titlesInh = {'Saline Inh','DOI Inh','5HT Inh'};
for mv=1:2
    figure
    if mv==1, set(gcf,'Name', 'mean LOW SF stationary')
    else set(gcf,'Name', 'mean LOW SF mv'), end
    for t=1:3
        lowSFpre= nanmean(drift_sf(data & treatment==t & hasDrift==1,2:4,mv,1),1)
        lowSFpost = nanmean(drift_sf(data & treatment==t & hasDrift==1,2:4,mv,2),1)
        subplot(2,3,t)
        plot(lowSFpre,'LineWidth',2);hold on; plot(lowSFpost,'r','LineWidth',2);axis square
        title(titles{t}); xlabel('SF'); ylabel('sp/sec');set(gca, 'XTickLabels', [2 3 4]);
        lowSFpreInh= nanmean(drift_sf(dataInh & treatment==t & hasDrift==1,2:4,mv,1),1)
        lowSFpostInh = nanmean(drift_sf(dataInh & treatment==t & hasDrift==1,2:4,mv,2),1)
        subplot(2,3,t+3)
        plot(lowSFpreInh,'LineWidth',2);hold on; plot(lowSFpostInh,'r','LineWidth',2);axis square
        title(titlesInh{t}); xlabel('SF'); ylabel('sp/sec');set(gca, 'XTickLabels', [2 3 4]);
    end
end

%%

titlesInh = {'Saline Inh','DOI Inh','5HT Inh'};
for mv=1:2
    figure
    if mv==1, set(gcf,'Name', 'mean High SF stationary')
    else set(gcf,'Name', 'mean High SF mv'), end
    for t=1:3
        highSFpre= nanmean(drift_sf(data & treatment==t & hasDrift==1,5:7,mv,1),1)
        highSFpost = nanmean(drift_sf(data & treatment==t & hasDrift==1,5:7,mv,2),1)
        subplot(2,3,t)
        plot(highSFpre,'LineWidth',2);hold on; plot(highSFpost,'r','LineWidth',2);axis square
        title(titles{t}); xlabel('SF'); ylabel('sp/sec');set(gca, 'XTickLabels', [5 6 7]);
        highSFpreInh= nanmean(drift_sf(dataInh & treatment==t & hasDrift==1,4:6,mv,1),1)
        highSFpostInh = nanmean(drift_sf(dataInh & treatment==t & hasDrift==1,4:6,mv,2),1)
        subplot(2,3,t+3)
        plot(highSFpreInh,'LineWidth',2);hold on; plot(highSFpostInh,'r','LineWidth',2);axis square
        title(titlesInh{t}); xlabel('SF'); ylabel('sp/sec');set(gca, 'XTickLabels', [5 6 7]);
        
    end
end

%%

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

%%
layerTet = layerSites(:,2:4:64)

salineSessions = find(sessionTreatment==1)
doiSessions = (sessionTreatment==2)
htSessions = sessionTreatment==3

doi4= layerTet(doiSessions,:)==4
saline4=layerTet(salineSessions,:)==4

% figure
% for s=1:length(salineSessions)
% subplot(5,5,s)
% imagesc(squeeze(LFPallCh(salineSessions(s),LFPallCh,:,1,2)))
% end



%%% is averaging across treatments for each layer

layerTet = layerSites(:,2:4:64)
figure
imagesc(squeeze(LFPallChDark(1,:,:,1,2)))
%
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
%
% use = find(htSessions==1)
% for prepost=1:2
%     figure
%     if prepost==1, set(gcf,'Name','5HT LFP PRE');
%     else set(gcf,'Name','5HT LFP POST'); end
%     for s =1:sum(htSessions)
%         subplot(4,4,s)
%         imagesc(squeeze(LFPallCh(use(s),:,:,1,prepost)))
%         xlabel('frequency (Hz)');ylabel('tetrode #')
%     end
% end
%
% clear s use
%use = find(salineSessions==1)
use = find(sessionTreatment==1)

for mv = 1:2
    figure
    if mv==1, set(gcf,'Name', 'moving');
    else set(gcf,'Name', 'stationary')
    end
    for prepost =1:2
        for s =  1:length(use)
            subplot(4,5,s)
            plot(squeeze(LFPfreq(:,prepost)),squeeze(nanmean(LFPallChDrift(use(s), :,:,mv,prepost))));
            hold on
        end
    end
end
%
% squeeze(nanmean(LFPallCh(sessionTreatment==1',layerTet(sessionTreatment==1,:)==4,:,1,1)))
%
% %%evoked LFP all
% % figure
% % for t=1:4
% % subplot(2,2,t)
% % set(gcf, 'Name', 'evoked LFP')
% % plot(squeeze(LFPall(salineSessions==1,:,1,1)),'b');hold on;
% % plot(squeeze(LFPall(salineSessions==1,:,1,2)),'r'); xlabel 'Frequency (Hz)'; ylabel 'normalized power';
% % end
% % %
% layerTet = layerSites(:,2:4:64)
% use = layerTet(:,:)==4
%
%
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
            plot(squeeze(nanmean(LFPallChDrift(use(s), layerTet(s,:)==4,6:13,mv,prepost))));
            hold on
        end
    end
end


% figure
% set(gcf, 'Name', 'all ch darkness LFP')
% for t=1:3
% subplot(1,3,t)
% plot(squeeze(mean(LFPallDark(sessionTreatment==t,:,1,1))),'b');hold on;
% plot(squeeze(mean(LFPallDark(sessionTreatment==t,:,1,2))),'r');% xlabel 'Frequency (Hz)'; ylabel 'normalized power';
% title(titles{t}); axis square
% end
%%
C = {[0 0 .9],[.8 0 0]}; %red = mv=1, light =pre, dark =post
D = {[0 1 0],[0 .5 0]}; %green = mv=2


%use = find(doiSessions==1)
use = find(sessionTreatment==1)

figure
for s=1:length(use)
    for prepost =1:2
        %for tet=1:16
        %subplot(4,4,s)
        plot(squeeze(LFPallChDrift(use(s),:,:,1,prepost))','Color',C{prepost});hold on; %xlim([0 120]);
        %plot(squeeze(LFPallChDark(9,tet,:,2,prepost)),'Color',D{prepost}); xlim([0 100]);
        xlabel 'Frequency (Hz)'; ylabel 'normalized power'; %mv
        % title(sprintf('layer %d',layerTet(tet)));
    end
end
%end

%
% figure
% subplot(1,2,1)
% imagesc(squeeze(LFPallChDark(1,layerTet==3,:,1,1)));
% subplot(1,2,2)
% imagesc(squeeze(LFPallChDark(1,layerTet==3,:,1,2))); xlabel 'Frequency (Hz)'; ylabel 'normalized power';
%
%%

layerTet = layerSites(:,2:4:64);

salineSessions = find(sessionTreatment==1);
doiSessions = (sessionTreatment==2);
htSessions = sessionTreatment==3;

mnTetDark = squeeze(nanmean(LFPallChDark(:,:,1:end-1,:,:),2));
mnTetDrift = squeeze(nanmean(LFPallChDrift(:,:,1:end-1,:,:),2));

% s = size(mnTetDark);
% mnTetDark = reshape(mnTetDark,[s(2:end) 1]);
% z = size(mnTetDrift);
% mnTetDrift = reshape(mnTetDrift,[z(2:end) 1]);

C = {[0 0 .9],[.8 0 0]}; %red = mv=1, light =pre, dark =post
D = {[0 1 0],[0 .5 0]}; %green = mv=2
for mv=1:2
    clear prenorm postnorm normrange mnpre mnpost
    figure; hold on;
    if mv==1
        set(gcf,'Name','stationary');
    else
        set(gcf,'Name','mv');
    end
    for t=1:2
        if t ==1
            use = find(sessionTreatment==1);
            use = use([1:3 6:13 15]);
        else
            use = find(sessionTreatment==2);
            use = use([1:2 4:6 8:12 14:15]);
        end
        
        for s= 1: length(use)
          %normalize response of all cells used, then average
            prerange = squeeze(mnTetDark(use(s),1:90,mv,1))
            normrange(s,:)= max(prerange)%-min(prerange)
            postrange = squeeze(mnTetDark(use(s),1:90,mv,2))
            normrangepost(s,:)= max(postrange)%-min(prerange)
            prenorm(s,:) =(mnTetDark(use(s),1:90,mv,1)-min(mnTetDark(use(s),1:90,mv,1)))... %
                ./normrange(s,1);
            
            postnorm(s,:) =(mnTetDark(use(s),1:90,mv,2)-min(mnTetDark(use(s),1:90,mv,2)))...
                ./normrangepost(s,1);
            
            prenorm(isinf(prenorm))=nan; postnorm(isinf(postnorm)) = nan ;
            
            mnpre= nanmean(prenorm,1);
            mnpost = nanmean(postnorm,1);
            
            nsessions = length(use);
            sempre = nanstd(prenorm)/sqrt(nsessions);
            sempost = nanstd(postnorm)/sqrt(nsessions);
            %             subplot(2,2,t)
            %             plot(prenorm','Color',C{1},'Linewidth',2); hold on;xlim([0 115]);
            %            plot(mnpre,'Linewidth',4,'Color',C{1});  axis square;% ylim([ 0 3.3e+04]);
            %            plot(postnorm','Color',C{2},'Linewidth',2);
            %            plot(mnpost,'Linewidth',4,'Color',C{2});
            
        end
        [h(t),p(t)] = kstest2(mnpre, mnpost)
        subplot(1,2,t)
        preerr=mseb(1:length(mnpre),mnpre,sempre,'k',1); hold on; axis square
        posterr=mseb(1:length(mnpost),mnpost,sempost,'r',1); xlim([0 90]);
        text(.1 ,.85, ['p = ' num2str(p(t))],'FontSize',20);
        ylim([0 1.0]);
        set(gca,'FontSize',24);set(gca, 'linewidth', 2)
        
        %frequencies are doubled
        
        delta = [nanmean(prenorm(:,1:4),2) nanmean(postnorm(:,1:4),2)]
        theta = [nanmean(prenorm(:,5:8),2) nanmean(postnorm(:,5:8),2)]
        alpha = [nanmean(prenorm(:,9:13),2) nanmean(postnorm(:,9:13),2)]
        gamma = [nanmean(prenorm(:,60:70),2) nanmean(postnorm(:,60:70),2)]
        
        %            delta = [nanmean(prenorm(1:8,1)) nanmean(postnorm(1:8,1))]
        %            theta = [nanmean(prenorm(9:16,1)) nanmean(postnorm(9:16,1))]
        %            alpha = [nanmean(prenorm(17:26,1)) nanmean(postnorm(17:26,1))]
        %            gamma = [nanmean(prenorm(120:140,1)) nanmean(postnorm(120:140,1))]
        
%          pDelta = kruskalwallis(delta)
%          pTheta = kruskalwallis(theta)
%          pAlpha = kruskalwallis(alpha)
%          pGamma = kruskalwallis(gamma)
[h pDelta(t)]=ttest(delta(:,1),delta(:,2))
[h pTheta(t)]=ttest(theta(:,1),theta(:,2))
[h pAlpha(t)]=ttest(alpha(:,1),alpha(:,2))
[h pGamma(t)]=ttest(gamma(:,1),gamma(:,2))




    end
end

%%

salineSessions = find(sessionTreatment==1);
doiSessions = (sessionTreatment==2);
htSessions = sessionTreatment==3;

mnTetDark = squeeze(nanmean(LFPallChDark(:,:,1:end-1,:,:),2));
mnTetDrift = squeeze(nanmean(LFPallChDrift(:,:,1:end-1,:,:),2));
C = {[0 0 .9],[.8 0 0]}; %red = mv=1, light =pre, dark =post
D = {[0 1 0],[0 .5 0]}; %green = mv=2
for mv=1
    clear prenorm postnorm normrange mnpre mnpost delta theta alpha beta gamma
    figure; hold on;
    if mv==1
        set(gcf,'Name','stationary');
    else
        set(gcf,'Name','mv');
    end
    for t=1:2
        if t ==1
            use = find(sessionTreatment==1);
            use = use([1:3 6:13 15]);
        else
            use = find(sessionTreatment==2);
            use = use([1:2 4:6 8:12 14:15]);
        end
        
        for s= 1: length(use)  
            %normalize response of all cells used, then average
            prerange = (mnTetDrift(use(s),1:90,mv,1));postrange = (mnTetDrift(use(s),1:90,mv,2))
            normrange(s,:)= max(prerange)%-min(prerange)
            normrangepost(s,:)= max(postrange);
            prenorm(s,:) =(mnTetDrift(use(s),1:90,mv,1)-min(mnTetDrift(use(s),1:90,mv,1)))... 
                ./normrange(s,1);
            
            postnorm(s,:) =(mnTetDrift(use(s),1:90,mv,2)-min(mnTetDrift(use(s),1:90,mv,2)))...
                ./normrangepost(s,1);
            
            prenorm(isinf(prenorm))=nan; postnorm(isinf(postnorm)) = nan ;
            
            mnpre= nanmean(prenorm,1);
            mnpost = nanmean(postnorm,1);
            
            nsessions = length(use);
            sempre = nanstd(prenorm)/sqrt(nsessions);
            sempost = nanstd(postnorm)/sqrt(nsessions);
            %            subplot(2,2,t)
            %            plot(prenorm','Color',C{1},'Linewidth',2); hold on;xlim([0 115]);
            %            plot(mnpre,'Linewidth',4,'Color',C{1});  axis square;% ylim([ 0 3.3e+04]);
            %            plot(postnorm','Color',C{2},'Linewidth',2);
            %            plot(mnpost,'Linewidth',4,'Color',C{2});    
        end
        [h(t),pks(t)] = kstest2(mnpre, mnpost)
%         subplot(1,2,t)
%        preerr=shadedErrorBar(1:length(mnpre),mnpre,sempre,'k',1); hold on; axis square      
%        posterr=shadedErrorBar(1:length(mnpost),mnpost,sempost,'r',1);
       
%         text(.1 ,.85, ['p = ' num2str(pks(t))],'FontSize',20);
%         xlim([0 90]);
%         ylim([0 1.0]);
        set(gca,'FontSize',24);set(gca, 'linewidth',2)
        %frequencies are doubled
        
%                    delta = [nanmean(prenorm(:,1:4),2) nanmean(postnorm(:,1:4),2)]
%                    theta = [nanmean(prenorm(:,5:8),2) nanmean(postnorm(:,5:8),2)]
%                    alpha = [nanmean(prenorm(:,9:13),2) nanmean(postnorm(:,9:13),2)]
%                    gamma = [nanmean(prenorm(:,60:70),2) nanmean(postnorm(:,60:70),2)]
        
                   %delta = [nanmean(prenorm(:,1:7),2) nanmean(postnorm(:,1:7),2)]
                   delta = [nanmean(prenorm(:,2:8),2) nanmean(postnorm(:,2:8),2)] %test lfo

                  % theta = [nanmean(prenorm(:,10:16),2) nanmean(postnorm(:,10:16),2)]
                   theta = [nanmean(prenorm(:,12:18),2) nanmean(postnorm(:,12:18),2)]

                   %alpha = [nanmean(prenorm(:,16:26),2) nanmean(postnorm(:,16:26),2)]
                   alpha = [nanmean(prenorm(:,18:26),2) nanmean(postnorm(:,18:26),2)]
                   
                   beta =  [nanmean(prenorm(:,26:56),2) nanmean(postnorm(:,26:56),2)]
                   gamma = [nanmean(prenorm(:,56:70),2) nanmean(postnorm(:,56:70),2)]
        
                %ttests for diff frequencies, avg across all animals pre vs
                %post
                    [o pdelta(t)]=ttest(delta(:,1),delta(:,2))
                    [o ptheta(t)]=ttest(theta(:,1),theta(:,2))
                    [o palpha(t)]=ttest(alpha(:,1),alpha(:,2))
                    [o pbeta(t)]=ttest(beta(:,1),beta(:,2))
                    [o pgamma(t)]=ttest(gamma(:,1),gamma(:,2))
                   
                  % delta = nanmean(delta); theta = nanmean(theta); alpha = nanmean(alpha);beta = nanmean(beta);gamma=nanmean(gamma);
                  %  freqs = [delta' theta' alpha' gamma'];
                  %  freqs = [delta' theta' alpha' beta' gamma'];
                  %  [h freqp(t)]= kstest2(freqs(1,:),freqs(2,:));

                    
%                     [o pdelta(t)]=ttest(delta(1),delta(2))
%                     [o ptheta(t)]=ttest(theta(1),theta(2))
%                     [o palpha(t)]=ttest(alpha(1),alpha(2))
%                     [o pbeta(t)]=ttest(beta(1),beta(2))
%                     [o pgamma(t)]=ttest(gamma(1),gamma(2))

%%% avg pts within each bin for each animal, to get 12x5 pts, then ttest
% [hD(t,mv) pDelta(t,mv)]=ttest(prenorm(~isnan(prenorm(:,1:8))),postnorm(~isnan(postnorm(:,1:8))))
% [hT(t,mv) pTheta(t,mv)]=ttest(prenorm(~isnan(prenorm(:,12:16))),postnorm(~isnan(postnorm(:,12:16))))
% [hA(t,mv) pAlpha(t,mv)]=ttest(prenorm(~isnan(prenorm(:,17:25))),postnorm(~isnan(postnorm(:,17:25))))
% [hB(t,mv) pBeta(t,mv)]=ttest(prenorm(~isnan(prenorm(:,26:59))),postnorm(~isnan(postnorm(:,26:59))))
% [hG(t,mv) pGamma(t,mv)]=ttest(prenorm(~isnan(prenorm(:,60:70))),postnorm(~isnan(postnorm(:,60:70))))

% subplot(1,2,t)
% hold on;
% plot(1,delta(:,1),'o')
% plot(2,delta(:,2),'o')
% 
% plot(4,theta(:,1),'o')
% plot(5,theta(:,2),'o')
% 
% plot(7,alpha(:,1),'o')
% plot(8,alpha(:,2),'o')
% 
% plot(10,beta(:,1),'o')
% plot(11,beta(:,2),'o')
% 
% plot(13,gamma(:,1),'o')
% plot(14,gamma(:,2),'o')
     end
end

