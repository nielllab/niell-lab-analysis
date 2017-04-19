clear all
close all
dbstop if error

batchDOIephys_filtered; %%% load batch file
%set(groot,'defaultFigureVisible','off') %disable figure plotting
set(groot,'defaultFigureVisible','on')

%%% select the sessions you want based on filters
%use =  find(strcmp({files.notes},'good data'))%useSess = use;
%use =  find( strcmp({files.treatment},'5HT') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}) )

%for specific experiment:
use =  find(strcmp({files.notes},'good data') & strcmp({files.expt},'030417'))
sprintf('%d selected sessions',length(use))

saline=1; doi=2; ht=3; ketanserin=4; ketandoi=5; mglur2=6; mglur2doi=7; lisuride=8;
savePDF=0;
redo = 1;
n=0; ncorr=0; %%% number of units loaded, ncorr= number of correlation pairs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
%     goodD = ones(1,nc);
%     goodD(files(use(i)).predark & files(use(i)).badsites) =0;
%     goodDark(cellrange)=goodD;
    
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
                LFPallDark(i,:,:,prepost) = NaN;

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% plot speed histogram
titles = {'Saline','DOI','5HT';
figure 
for t = 1:3
subplot(1,3,t)
plot(0.5:1:25,squeeze(mean(speedHistWn(sessionTreatment==t,1:25,:),1))); title(titles{t}); xlabel('speed');
axis square; 
end

legend('pre','post')

%%%filter data before plotting%%%
clear max
low_wn = squeeze(min(wn_crf,[],2));
max_wn = squeeze(max(wn_crf,[],2));
amp_wn = max_wn-low_wn
useResp = amp_wn(:,1,1)>0& amp_wn(:,1,2)>0 & amp_wn(:,2,1)>0 & amp_wn(:,2,2)>0;
data_wn = goodAll==1 & useResp';

clear max amp low
peakresp = squeeze(max(drift_orient,[],2));
low = squeeze(min(drift_orient,[],2));
amp = peakresp-low;
useResp = amp(:,1,1)>1 & amp(:,1,2)>1 |amp(:,2,1)>1 & amp(:,2,2)>1;
data = goodAll & useResp' & ~inhAll ;%which cells to include
dataInh =goodAll & useResp' & inhAll ;

%evoked FR from drift
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
        
        mdl_drift = fitlm(allPre,allPost)
        rsquared_drift(t) = mdl_drift.Rsquared.Ordinary   
        plot(mpre,mpost,'+k','Markersize',10,'Linewidth',2)
        plot([0 35], [0 35])
        plot(mean(drift_orient(goodAll==1 & inhAll==1  &treatment==t,:,mv,1),2),mean(drift_orient(goodAll==1 & inhAll==1 &treatment==t,:,mv,2),2),'r.','Markersize',10);
        n_cells(t) = sum(goodAll==1 &treatment==t)
        text(.5, 11.5, ['n = ' num2str(n_cells(t))],'FontSize',18)
         text(.5, 10.75, ['r^2 = ' num2str(rsquared_drift(t))],'FontSize',18)
    end
end
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
    useN=goodAll & treatment==t
    miDrift= (mean(drift_orient(useN &(layerAll==2|layerAll==3),:,1,2),2)-mean(drift_orient(useN&(layerAll==2|layerAll==3),:,1,1),2))./(mean(drift_orient(useN&(layerAll==2|layerAll==3),:,1,2),2)+mean(drift_orient(useN&(layerAll==2|layerAll==3),:,1,1),2));
    h= hist(miDrift,-1:.2:1);
    Mbins=-1:.2:1
    subplot(1,3,1)
    bar(Mbins,h/sum(useN &(layerAll==2|layerAll==3)),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .4]); xlim([-1.5 1.5]);axis square; title('layer 2/3')
    for i=4:5
        miDrift= (mean(drift_orient(useN &(layerAll==i),:,1,2),2)-mean(drift_orient(useN&(layerAll==i),:,1,1),2))./(mean(drift_orient(useN&(layerAll==i),:,1,2),2)+mean(drift_orient(useN&(layerAll==i),:,1,1),2));
        h= hist(miDrift,-1:.2:1);
        subplot(1,3,i-2)
        Mbins=-1:.2:1
        bar(Mbins,h/sum(useN&layerAll==i),'FaceColor',[0 .5 .5],'Linewidth',2);ylim([0 .4]);xlim([-1.5 1.5]); axis square
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

%MI for drift spont
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



clear cycR
for f = 1:20;
    cycR(:,f,:,:) = nanmean(wn_frameR(:,f:20:end,:,:),2);
end
spont = squeeze(mean(cycR(:,[1 2 19 20],:,:),2));
evoked = squeeze(mean(cycR(:,9:11,:,:),2)) - spont;

%wn evoked mv and stationary
titles = {'Saline','DOI','5HT'}
for mv=1:2
if mv==1, set(gcf,'Name', 'stationary');
else set(gcf,'Name', 'mv');end
figure
for t=1:3
    subplot(1,3,t)
    set(gcf,'Name','wn evoke prepost')
    plot(wn_evoked(data_wn&treatment==t&~inhAll,mv,1),wn_evoked(data_wn&treatment==t&~inhAll,mv,2),'.');hold on;axis square;
    xlim([-10 30]); ylim([-10 30]);
    plot([-30 30],[-30 30]);
    plot(wn_evoked(data_wn&treatment==t&inhAll==1,mv,1),wn_evoked(data_wn&treatment==t&inhAll==1,mv,2),'.r');
    title(titles{t})
end
end

titles = {'Saline','DOI','5HT'}
for mv=1:2
if mv==1, set(gcf,'Name', 'stationary');
else set(gcf,'Name', 'mv');end
figure
for t=1:3
    subplot(1,3,t)
    set(gcf,'Name','wn evoke prepost')
    plot(wn_spont(data_wn&treatment==t&~inhAll,mv,1),wn_spont(data_wn&treatment==t&~inhAll,mv,2),'.');hold on;axis square;
    xlim([min(wn_spont(data_wn,mv,1))-.5 max(wn_spont(data_wn,mv,1))+.5]); ylim([min(wn_spont(data_wn,mv,1))-.5 max(wn_spont(data_wn,mv,1))+.5]);
    plot([-30 30],[-30 30]);
    plot(wn_spont(data_wn&treatment==t&inhAll==1,mv,1),wn_spont(data_wn&treatment==t&inhAll==1,mv,2),'.r');
    title(titles{t})
end
end

%evoked calculcated from cyc avg
for mv = 1:2
    figure
    if mv==1, set(gcf,'Name', 'evoked pre vs post stop')
    else set(gcf,'Name', 'evoked pre vs post mv'), end
    for t=1:3
        subplot(1,3,t)
        plot(evoked(find(data& treatment==t),mv,1),evoked(find(data& treatment==t),mv,2),'.'); title(titles{t}); hold on; plot([0 50],[0 50]); axis square
        xlabel('pre'); ylabel('post'); xlim([min(evoked(data_wn,mv,1))-.5 max(evoked(data_wn,mv,1))+.5]); ylim([min(evoked(data_wn,mv,1))-.5 max(evoked(data_wn,mv,1))+.5]);
        plot(evoked(find(dataInh& treatment==t),mv,1),evoked(find(dataInh& treatment==t),mv,2),'.r'); hold on; plot([0 50],[0 50]); %axis([-5 10 -5 10])
    end
end

%spont calculcated from cyc avg
titles = {'saline','doi','ht'};
for mv = 1:2 
figure 
if mv==1, set(gcf,'Name', 'spont pre vs post stop')
else set(gcf,'Name', 'spont pre vs post mv'), end
for t=1:3
    subplot(1,3,t)
    plot(spont(find(data & treatment==t),mv,1),spont(find(data& treatment==t),mv,2),'.'); title(titles{t}); hold on; axis square; plot([0 50],[0 50]);
    xlabel('pre'); ylabel('post');
    plot(spont(find(dataInh&inhAll& treatment==t),mv,1),spont(find(dataInh & inhAll& treatment==t),mv,2),'.r'); hold on; plot([0 50],[0 50]); axis([0 10 0 10])
end
end


% clear data
% titles={'Saline', 'DOI', '5HT'}
% data=(meanRdark)
% figure
% for t=1:3
%     subplot(1,3,t)
%     plot(meanRdark(~inhAll & treatment==t,1), meanRdark(~inhAll &treatment==t,2),'.','Markersize',12);hold on; axis square;
%     set(gca,'FontSize',18)
%    % plot(meanRdark(inhAll  & treatment==t,1), meanRdark(inhAll &treatment==t,2),'r.','Markersize',12);hold on; axis square;
% %     mdl_dark= fitlm(meanRdark(treatment==t,1), meanRdark(treatment==t,2),2)
% %     rsquared_dark(t) = mdl_dark.Rsquared.Ordinary
% %     xlim([0 5]);ylim([0 5]); %xlabel('Pre spikes/sec');ylabel('Post spikes/sec');
%     % mpre=nanmean(mean(meanRdark(goodAll==1 & treatment==t,1),2))
%     % mpost=nanmean(mean(meanRdark(goodAll==1 & treatment==t,2),2))
%     % plot(mpre,mpost,'+k','Markersize',12,'Linewidth',2)
%     plot([0 35],[0 35],'Linewidth',2);
%  %   n_cells(t) = sum(goodAll==1 &treatment==t)
%     %     text(.25, 14, ['n = ' num2str(n_cells(t))],'FontSize',18)
%     %     text(.25, 13, ['r^2 = ' num2str(rsquared_dark(t))],'FontSize',18)
%     %title(titles{t});
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
for t=1:3
    preDriftCorr = driftCorr(corrTreatment==t,1)
    postDriftCorr=driftCorr(corrTreatment==t,2)
    mdl = fitlm(preDriftCorr,postDriftCorr)
    rsquared(t) = mdl.Rsquared.Adjusted
end

titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI','MGluR2', 'MGluR2 + DOI','Lisuride'};
figure
for i = 1:3
    subplot(1,3,i);
    plot(driftCorr(corrTreatment==i,1),driftCorr(corrTreatment==i,2),'.'); hold on; axis equal
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
    xlabel('pre drift corr'); ylabel('post')
    text(-.25, .85, ['r^2 = ' num2str(rsquared(i))])
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
for t=1:3
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
    plot([-0.5 1],[-0.5 1]); axis([-0.5 1 -0.5 1]); title(titles{i});
    xlabel('pre dark corr'); ylabel('post')
     text(-.25, .85, ['r^2 = ' num2str(rsquared(i))])
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

%% nx & ny prepost
figure
titles={'Saline','DOI','5HT'};
useSTA = sta_exp_var(:,1) & sta_exp_var(:,2)>.6
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

useSTA = sta_exp_var(:,1) & sta_exp_var(:,2)>.6
useSTA = useSTA & data_wn'% &(layerAll==4)'
figure
titles={'Saline', 'DOI','5HT'}
for t=1:3
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
for t=1:3
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


useSTA = find((sta_exp_var(:,1) & sta_exp_var(:,2)>.5) &(treatment==doi)' &(layerAll==2|3)');
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
useSTA = find((sta_exp_var(:,1) & sta_exp_var(:,2)>.5) &(treatment==doi)' &(layerAll==2|3)');
for prepost=1:2
    figure
    if prepost==1
        set(gcf,'Name','pre')
    else
        set(gcf,'Name','post');end
    for c=1:length(useSTA)
        subplot(10,7,c)
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
            title(titles{l});

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
 
 
%% to do: seperate to have multiple figures for each treatment
for t = 1:3
    figure
    if t==1, set(gcf,'Name','saline wn CRF'),
    elseif t==2, set(gcf,'Name','doi wn CRF'),
    else set(gcf, 'Name','ht wn CRF'); end
%     elseif t==3, set(gcf,'Name','ht wn CRF'),
%     elseif t==4, set(gcf,'Name','ketanserin wn CRF')
%     elseif t==5, set(gcf,'Name','ketanserin + doi wn CRF')
%     elseif t==6, set(gcf,'Name','MGluR2 wn CRF')
%     elseif t==7, set(gcf, 'Name','MGluR2+ DOI wn CRF')
%     else set(gcf,'Name','Lisuride'),end

    useN = find(goodAll & treatment==t)%& layerAll==4)
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
useN =data
for c= 1:length(useN)
    for prepost =1:2
        [respmax oind] = max(drift_orient(c ,:,1,prepost)');
        [g h]=(max(respmax(:)));
        prefOri(prepost,c)=squeeze(oind(h));
    end
end
% size(prefOri')
prefOri_pre = prefOri(1,:)'
prefOri_post = prefOri(2,:)'

titles ={'Saline','DOI','5HT'};
figure
for t=1:3
    useOri_pre= prefOri_pre(treatment==t)
    useOri_post= prefOri_post(treatment==t)
    %     useOri_post= find(prefOri_post & (treatment==t)')
    subplot(2,3,t)
    hist(useOri_pre);xlabel('pre orientation pref'); ylabel('proportion of cells')
    xlim([0 12]); ylim([0 .03]);
    Mbins= 0:1:12
    hpre= hist(useOri_pre,0:1:12);
    bar(Mbins,hpre/sum(useOri_pre)); xlim([0 13]);ylim([0 .035]);axis square;
    title(titles{t});
    subplot(2,3,t+3)
    hpost=hist(useOri_post,0:1:12);
    bar(Mbins,hpost/sum(useOri_post))
    xlabel('post orientation preference');ylabel('proportion of cells')
    axis square
    xlim([0 13]);ylim([0 .035]);
%     [p h]=ranksum(prefOri_pre(treatment==t),prefOri_post(treatment==t))
%     text(6, .0275, ['p = ' num2str(p)],'FontSize',18)
end



figure
for t=1:3
    clear useN
    useN = data & treatment==t
    for c = 1:sum(useN)
        subplot(1,3,t)
        mn = squeeze(drift_cond_tcourse(useN(c),1,1,prefOri_pre(useN(c)),1,:))
        mn = mn - repmat(mn(1,:),[50 1]);
        mn = circshift(mn,10);
        subplot(1,3,t)
        plot((1:length(mn)-5)*dt -dt/2,mn(1:45,:)) %,'LineWidth',2);axis square; set(gca,'fontsize', 18);
    end
end

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
     
        title(titles{t}); xlabel('time (s)');
        ylabel('spikes/sec');
        %ylim([0 max(mn(:))+1])
    end
end

titles = {'Saline','DOI','5HT','ketanserin', 'ketanserin + DOI', 'MGluR2','MGluR2 + DOI','Lisuride'};
figure
for t=1:3
    mnpre = squeeze(mean(mnPsth(data & layerAll==4 & treatment==t,1,:),1))';
    mnpost = squeeze(mean(mnPsth(data & layerAll==4 &treatment==t,2,:),1))';
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

%osi all layers
figure
for t=1:3
    osi(1,:) = nanmean(abs(drift_osi(data&treatment==t,1,1)))
    osierr(1,:) = squeeze(nanstd(drift_osi(data & treatment==t,1,1)/sqrt(sum(drift_osi(data & treatment==t,1,1),1))));
    osi(2,:) = nanmean(abs(drift_osi(data&treatment==t,1,2)))
    osierr(2,:) = squeeze(nanstd(drift_osi(data & treatment==t,1,2)/sqrt(sum(drift_osi(data & treatment==t,1,2),1))));
    subplot(1,3,t)
    barweb(osi,osierr); ylim([0 1]); axis square
    [p h]=ranksum(drift_osi(data&treatment==t,1,1),drift_osi(data&treatment==t,1,2))
    text(.9, 1, ['p = ' num2str(p)],'FontSize',18)
end

%osi specific layers
figure
for t=1:3
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
for t = 1:3
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
for t = 1:3
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


