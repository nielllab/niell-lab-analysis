clear all
close all

dbstop if error

batchDOIephys_filtered; %%% load batch file

%%% select the sessions you want based on filters
%%% example
use =  find(strcmp({files.notes},'good data')& ~cellfun(@isempty,{files.predark})& ~cellfun(@isempty,{files.postdark}) & ~cellfun(@isempty,{files.prewn}) & strcmp({files.treatment},'DOI')  )

%use =  find( strcmp({files.treatment},'KetanserinDOI') & strcmp({files.notes},'good data') & ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}) )
%use =  find( strcmp({files.treatment},'DOI') &  ~cellfun(@isempty,{files.predark}) & ~cellfun(@isempty,{files.postdark}))

%for specific experiment:
%use =  find(strcmp({files.notes},'bad data')  & ~cellfun(@isempty,{files.predark})& ~cellfun(@isempty,{files.postdark}) & strcmp({files.expt},'083115'))
sprintf('%d selected sessions',length(use))

useSess = use;

saline=1; doi=2; ketanserin=3; ketandoi=4; lisuride=5;

savePDF=0;
redo = 1;
n=0; ncorr=0; %%% number of units loaded, ncorr= number of correlation pairs
for i = 1:length(use)
    
    %%% extract filenames
    afile = [pathname '\' files(use(i)).dir '\' files(use(i)).analysisfile '.mat'];
    clustfile = [pathname '\' files(use(i)).dir '\' files(use(i)).clusterfile '.mat'] ;
    cfile = [{[pathname '\' files(use(i)).dir '\' files(use(i)).predark_camera '.mat']}; {[pathname '\' files(use(i)).dir '\' files(use(i)).postdark_camera '.mat']}]';
    
    %%% get cell type based on waveform
    [inh mid] = getWaveform(clustfile,afile,0);
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
        layerInfo = getLayer(clustfile,afile,files(use(i)).tip1,files(use(i)).tip2,files(use(i)).angle, 1); %(needs histo information, but will give layers for all sites)
        layerAll(cellrange) = layerInfo.units;
        layerSites= layerInfo.sites;
    end
    %
    %     if ~isempty(files(use(i)).blockWn{1}) & ~isempty(files(use(i)).blockWn{2})
    %
    %         for prepost =1:2
    %             eyes = getEyes(clustfile,afile,cfile{:,prepost}, files(use(i)).blockWn{prepost},1);
    %             rad{:,prepost} = eyes.rad
    %             t{prepost} = eyes.t
    %         end
    %     end
    
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
            spd = getSpeed(clustfile,afile,files(use(i)).blockWn{prepost},0);
            speedHistWn(i,:,prepost) = hist(spd.v,0.5:1:100)/length(spd.v);
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
                drift_trial_psth{cellrange(c),prepost} = squeeze(drift.trialPsth(c,:,:));
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
        
        %%% get wn response
        for prepost = 1:2
            wn = getWn_mv(clustfile,afile,files(use(i)).blockWn{prepost},0,300);
            wn_crf(cellrange,:,:,prepost)=wn.crf;
            wn_spont(cellrange,:,prepost)=wn.spont;
            wn_evoked(cellrange,:,prepost)=wn.evoked;
            wn_gain(1,cellrange,prepost)=wn.gain;
            wn_frameR(cellrange,:,:,prepost) = downsamplebin(wn.frameR,2,15)/15;
        end
    end
    
    figure
    subplot(2,2,1); imagesc(squeeze(wn_frameR(cellrange,:,1,1)),[ 0 50]); title('pre stop')
    subplot(2,2,2); imagesc(squeeze(wn_frameR(cellrange,:,1,2)),[ 0 50]); title('post stop')
    subplot(2,2,3); imagesc(squeeze(wn_frameR(cellrange,:,2,1)),[ 0 50]); title('pre move');
    subplot(2,2,4); imagesc(squeeze(wn_frameR(cellrange,:,2,2)),[ 0 50]); title('post move')
    set(gcf,'Name',sprintf('%s %s',files(use(i)).expt, files(use(i)).treatment));
    drawnow
    
    %   getSorting(clustfile,afile,sprintf('%s %s',files(use(i)).expt,files(use(i)).treatment));
    %   drawnow
    %
    
    %% lfp power
    %%%(right now averages over all sites, should use layer info)
    %     for prepost=1:2
    %         lfpMove = getLfpMovement(clustfile,afile,files(use(i)).blockWn{prepost},1);
    %         LFPall(i,:,:,prepost) =squeeze(nanmedian(lfpMove.meanSpect, 1));
    %         size(lfpMove.meanSpect)
    %         if size(lfpMove.meanSpect,1)==16
    %             LFPallCh(i,:,:,:,prepost) = lfpMove.meanSpect;
    %         end
    %     end
    
    %     for prepost=1:2
    %         lfpMoveDark = getLfpMovement(clustfile,afile,files(use(i)).blockDark{prepost},0);
    %         LFPallDark(i,:,:,prepost) =squeeze(nanmedian(lfpMove.meanSpect, 1));
    %         if size(lfpMove.meanSpect,1)==16
    %             LFPallChDark(i,:,:,:,prepost) = lfpMove.meanSpect;
    %         end
    %     end
    %
    %     %% darkness / correlation analysis
    %
    %     %%% need to keep track of n^2 values for correlations
    %     corrRange=ncorr+1:ncorr+nc^2;
    %     corrTreatment(corrRange)=treatment(cellrange(1));
    %
    %
    %     %%% prepost correlation for white noise
    %     dt = 1;
    %     [preCorr postCorr cv2 R eigs] =  prepostDOIdarkness(clustfile,afile,files(use(i)).blockWn,dt,0);
    %     wnCorr(corrRange,1) = preCorr(:); wnCorr(corrRange,2)=postCorr(:);
    %
    %     %%%% prepost correlation in darkness
    %     dt = 1;
    %     [preCorr postCorr cv2 R eigs] = prepostDOIdarkness(clustfile,afile,files(use(i)).blockDark,dt,0);
    %     darkCorr(corrRange,1) = preCorr(:); darkCorr(corrRange,2)=postCorr(:);
    %
    %     cv2Dark(cellrange,:) = cv2;
    %     meanRdark(cellrange,:) = mean(R,2);
    %
    %     %%% keep track of cell type for correlations
    %     corrType1 = zeros(size(preCorr)); corrType2 = corrType1;
    %     for j= 1:length(inh);
    %         corrType1(j,:)=inh(j); corrType2(:,j)=inh(j);
    %     end
    %     corrType1all(corrRange) = corrType1(:) ; corrType2all(corrRange)= corrType2(:);
    %
    n= n+nc;
    
    %     ncorr= ncorr+nc^2;
end

keyboard

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

    clear ori_data
    
    %%% average across reps and both directions of motion
    %%% ori_data(cell,mv,prepost,orientation,t)
    
    for ori = 1:6
        ori_data(:,:,:,ori,:) = squeeze(nanmean(nanmean(drift_cond_tcourse(:,:,:,[ori ori+6],:,:),5),4));
    end
    
    for n= 1:size(ori_data,1);
        [y peak] = max(squeeze(mean(mean(mean(ori_data(n,:,:,:,5:30),5),3),2)));
        ori_data(n,:,:,:,:) = circshift(ori_data(n,:,:,:,:),-peak,4);
    end
    
    ori_data =downsamplebin(ori_data,5,2,1);
    ori_data = circshift(ori_data,5,5);
    for prepost=1:2
        for mv = 1:2
            subplot(2,2,prepost+2*(mv-1));
            plot(squeeze(nanmean(ori_data(find(goodAll & select),mv,prepost,:,:),1))');
            title(titles{prepost+2*(mv-1)}); if c<=2; ylim([0 4]); else ylim([0 10]); end
        end
    end
end



mv = 1
clear ori_data

%%% average over repetitions and opposite directions of motion
%%% oridata(cell,prepost,orientation, time)
for ori = 1:6
    ori_data(:,:,ori,:) = squeeze(nanmean(nanmean(drift_cond_tcourse(:,mv,:,[ori ori+6],:,:),5),4));
end

%%% align by peak orientation
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

dsmall = downsamplebin(d,2,5,1);
dist = pdist(dsmall(:,1:end/2),'correlation');  %%% sort based on correlation coefficient
display('doing cluster')
tic, Z = linkage(dist,'ward'); toc
figure
subplot(3,4,[1 5 9 ])
display('doing dendrogram')
[h t perm] = dendrogram(Z,0,'Orientation','Left','ColorThreshold' ,5);
axis off
subplot(3,4,[2 3 4 6 7 8 10 11 12 ]);
imagesc((d(perm,:)),[0 10]); axis xy   %%% show sorted data

lyr = layerAll(goodAll & clean'); sess = sessionNum(goodAll & clean'); type = celltype(goodAll & clean');
figure
imagesc(lyr(perm)'); title('layers clustered'); colormap jet; axis xy
figure
imagesc(sess(perm)'); title('sessions clustered'); colormap jet; axis xy
figure
imagesc(type(perm)'); title('type clustered'); colormap jet; axis xy

figure
imagesc(layerAll(goodAll & clean')'); title('layer by session'); colormap jet;

figure
imagesc(inhAll(goodAll & clean')'); title('inh by session'); colormap jet;

celltype = layerAll; celltype(inhAll) = 7;
[y order] = sort(celltype(find(clean & goodAll')));

figure
imagesc(resp(order,:),[0 1]); hold on; title('layers');
borders = find(diff(y)); for i = 1:length(borders); plot([1 size(d(:,:),2)],[borders(i) borders(i)],'m','Linewidth',2);end

figure
imagesc(resp,[0 1]); hold on; borders = find(diff(sessionNum(find(goodAll & clean')))); for i = 1:length(borders); plot([1 size(d(:,:),2)],[borders(i) borders(i)],'m','Linewidth',2);end
title('sessions')

figure
imagesc(d(order,:),[0 10]); hold on; title('layers');
borders = find(diff(y)); for i = 1:length(borders); plot([1 size(d(:,:),2)],[borders(i) borders(i)],'m','Linewidth',2);end

figure
imagesc(d(:,:),[0 10]); hold on; borders = find(diff(sessionNum(find(goodAll & clean')))); for i = 1:length(borders); plot([1 size(d(:,:),2)],[borders(i) borders(i)],'m','Linewidth',2);end
title('sessions')



tic
[W,V,whichMarg] = dpca(data, 20, ...
    'combinedParams', combinedParams);
toc

explVar = dpca_explainedVariance(data, W, V, ...
    'combinedParams', combinedParams);


time= (1:25)*0.1; timeEvents = 0.5;
dpca_plot(data, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);


figure
imagesc(V(order,1:15),[-0.5 0.5]); colormap jet; title('encoder by layer');

keyboard



d= squeeze(drift_tcourse(:,:,1,:,1));  %%% stationary first timepoint
d = downsamplebin(d,2,2,1);
d = reshape(d,size(d,1),size(d,2)*size(d,3));
nanratio = mean(isnan(d),2);
figure
plot(nanratio)

for prepost = 1:2
    
    tcourse = squeeze(drift_tcourse(:,:,1,prepost,:));
    tcourse = downsamplebin(tcourse(:,:,2:49),3,4,1);
    tcourse = downsamplebin(tcourse,2,2,1);
    tcourse = permute(tcourse,[1 3 2]);
    trace = reshape(tcourse,size(tcourse,1),size(tcourse,2)*size(tcourse,3));
    
    figure
    imagesc(isnan(squeeze(trace)))
    
    if prepost ==1
        s = nanstd(trace,[],2);
    end
    figure
    hist(s,0.25:0.5:15)
    usespks = find(s>1 & nanratio<0.04);
    figure
    imagesc(isnan(squeeze(trace(usespks,:))));
    
    mn = nanmean(trace,2);
    
    normtrace = (trace - repmat(mn,[1 size(trace,2)]))./repmat(s,[1 size(trace,2)]);
    figure
    imagesc(normtrace(usespks,:),[ -5 5])
    
    if prepost==1
        [coeff score latent] = pca(normtrace(usespks,:)');  score(:,12) = mean(normtrace(usespks,:),1);
        %[coeff score latent] = pca(normtrace(usespks,:)');
        
        figure
        plot(latent(1:20)/sum(latent));
        
    end
    data = normtrace(usespks,:)';
    score = data*coeff;
    
    figure
    for i = 1:5;
        subplot(5,1,i);
        plot(score(:,i));
    end
    
    condtuning = reshape(score',size(score,2),size(tcourse,2),size(tcourse,3));
    tuning = reshape(condtuning,size(condtuning,1),size(condtuning,2),3,12);
    
    figure
    for i = 1:12
        subplot(3,4,i)
        d=circshift(squeeze(nanmean(tuning(i,:,:,:),4)),2); %%% average over either sf or orientation
        plot(d);  range = max(abs(d(:)))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 size(d,1)+0.5])
    end
    
    figure
    for i = 1:12
        subplot(3,4,i)
        d=circshift(squeeze(mean(tuning(i,:,:,:),3)),2);
        plot(d);  range = max(abs(d(:)))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 size(d,1)+0.5])
    end
    
    figure
    for i = 1:12
        subplot(3,4,i)
        d=squeeze(mean(mean(tuning(i,2:7,:,:),3),2))- squeeze(mean(mean(tuning(i,10:12,:,:),3),2));
        plot(d); range = max(abs(d))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 12.5])
    end
    
    figure
    for i = 1:12
        subplot(3,4,i)
        d=squeeze(mean(tuning(i,2:7,:,:),2))- squeeze(mean(tuning(i,10:12,:,:),2));
        plot(d'); range = max(abs(d(:)))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 12.5])
    end
    
    figure
    for i = 1:12
        subplot(3,4,i)
        d=squeeze(mean(mean(tuning(i,1,:,:),3),2))- squeeze(mean(mean(tuning(i,10:12,:,:),3),2));
        plot(d); range = max(abs(d))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 12.5])
    end
    
end

tcourse = squeeze(drift_tcourse(:,:,1,:,:));
tcourse = downsamplebin(tcourse(:,:,:,2:49),4,4,1);
tcourse = downsamplebin(tcourse,2,2,1);
tcourse = permute(tcourse,[1 4 2 3]);
trace = reshape(tcourse,size(tcourse,1),size(tcourse,2)*size(tcourse,3)*size(tcourse,4));

figure
imagesc(isnan(squeeze(trace)))

if prepost ==1
    s = nanstd(trace,[],2);
end
figure
hist(s,0.25:0.5:15)
usespks = find(s>1 & nanratio<0.04);
figure
imagesc(isnan(squeeze(trace(usespks,:))));

mn = nanmean(trace,2);

normtrace = (trace - repmat(mn,[1 size(trace,2)]))./repmat(s,[1 size(trace,2)]);
figure
imagesc(normtrace(usespks,:),[ -5 5])

[coeff scoreAll latent] = pca(normtrace(usespks,:)');
%[coeff score latent] = pca(normtrace(usespks,:)');


figure
for i = 1:5;
    subplot(5,1,i);
    plot(scoreAll(:,i));
end

for prepost = 1:2
    if prepost==1
        score = scoreAll(1:size(scoreAll,1)/2,:);
    else
        score = scoreAll((size(scoreAll,1)/2+1):end,:);
    end
    
    condtuning = reshape(score',size(score,2),size(tcourse,2),size(tcourse,3));
    tuning = reshape(condtuning,size(condtuning,1),size(condtuning,2),3,12);
    
    
    for i = 1:12
        subplot(3,4,i)
        d=circshift(squeeze(nanmean(tuning(i,:,:,:),4)),2);
        plot(d);  range = max(abs(d(:)))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 size(d,1)+0.5])
    end
    
    figure
    for i = 1:12
        subplot(3,4,i)
        d=circshift(squeeze(mean(tuning(i,:,:,:),3)),2);
        plot(d);  range = max(abs(d(:)))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 size(d,1)+0.5])
    end
    
    figure
    for i = 1:12
        subplot(3,4,i)
        d=squeeze(mean(mean(tuning(i,2:7,:,:),3),2))- squeeze(mean(mean(tuning(i,10:12,:,:),3),2));
        plot(d); range = max(abs(d))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 12.5])
    end
    
    figure
    for i = 1:12
        subplot(3,4,i)
        d=squeeze(mean(tuning(i,2:7,:,:),2))- squeeze(mean(tuning(i,10:12,:,:),2));
        plot(d'); range = max(abs(d(:)))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 12.5])
    end
    
    figure
    for i = 1:12
        subplot(3,4,i)
        d=squeeze(mean(mean(tuning(i,1,:,:),3),2))- squeeze(mean(mean(tuning(i,10:12,:,:),3),2));
        plot(d); range = max(abs(d))*1.25; range = max(range,4); ylim([-range range]); xlim([0.5 12.5])
    end
end



%%% plot mean timecourse across all drift stim and cells, by layer / cell-type
%%% this should be updated to select active cells and separate
%%% moving/stationary. Also maybe choose optimal stim for each cell?
dt = 0.05;
for i=1:3
    if i==1
        mn = squeeze(mean(mnPsth(find(goodAll & ~inhAll & layerAll==5),:,:),1))'; t = 'grating lyr 5';
    elseif i==2
        mn = squeeze(mean(mnPsth(find(goodAll & ~inhAll & layerAll<5),:,:),1))'; t = 'grating lyr 2-4';
    elseif i==3
        mn = squeeze(mean(mnPsth(find(goodAll & inhAll),:,:),1))'; t = 'grating inh';
    end
    mn = mn - repmat(mn(1,:),[50 1]);
    mn = circshift(mn,10)
    figure
    plot((1:length(mn)-5)*dt -dt/2,mn(1:45,:)); title(t); xlabel('secs'); ylabel('sp/sec')
    %ylim([0 max(mn(:))+1])
end





figure
subplot(2,2,1); imagesc(squeeze(wn_frameR(find(goodAll),:,1,1)),[ 0 50]); title('pre stop')
subplot(2,2,2); imagesc(squeeze(wn_frameR(find(goodAll),:,1,2)),[ 0 50]); title('post stop')
subplot(2,2,3); imagesc(squeeze(wn_frameR(find(goodAll),:,2,1)),[ 0 50]); title('pre move');
subplot(2,2,4); imagesc(squeeze(wn_frameR(find(goodAll),:,2,2)),[ 0 50]); title('post move')
set(gcf,'Name',sprintf('%s %s',files(use(i)).expt, files(use(i)).treatment));
drawnow


figure
subplot(2,2,1); plot(squeeze(nanmean(wn_frameR(find(goodAll),:,1,1),1))); title('pre stop'); ylim([0 5])
subplot(2,2,2); plot(squeeze(nanmean(wn_frameR(find(goodAll),:,1,2),1))); title('post stop');ylim([0 5])
subplot(2,2,3); plot(squeeze(nanmean(wn_frameR(find(goodAll),:,2,1),1))); title('pre move'); ylim([0 5])
subplot(2,2,4); plot(squeeze(nanmean(wn_frameR(find(goodAll),:,2,2),1))); title('post move'); ylim([0 5])
set(gcf,'Name',sprintf('%s %s',files(use(i)).expt, files(use(i)).treatment));
drawnow

clear cycR
for f = 1:20;
    cycR(:,f,:,:) = nanmean(wn_frameR(:,f:20:end,:,:),2);
end

spont = squeeze(mean(cycR(:,[1 2 19 20],:,:),2));
evoked = squeeze(mean(cycR(:,9:11,:,:),2)) - spont;

clear spontbar evokedbar spont_err evoked_err

figure
for lyr = 2:5
    for i=1:2
        if i==1
            use = find(goodAll & layerAll==lyr & ~inhAll & (squeeze(evoked(:,2,1)>1)' & squeeze(evoked(:,2,2)>1)') ); symb = 'bo';
        else
            use = find(goodAll & layerAll==lyr  & inhAll ); symb = 'ro';
        end
        subplot(4,4,1 + 4*(lyr-2));
        plot(spont(use,1,1),spont(use,1,2),symb); hold on; plot([0 50],[0 50]);  axis square; axis([0 30 0 30])
        if lyr==2,  title({'spont stop','layer 2'}), end; if lyr>2, title(sprintf('layer %d',lyr)); end
        subplot(4,4,2+ 4*(lyr-2));
        plot(spont(use,2,1),spont(use,2,2),symb);  hold on; plot([0 50],[0 50]); axis square; axis([0 30 0 30])
        if lyr==2,  title({'spont move',''}), end
        subplot(4,4,3+ 4*(lyr-2));
        plot(evoked(use,1,1),evoked(use,1,2),symb); hold on; plot([-30 30],[-30 30]);  axis square; axis([-10 20 -10 20])
        if lyr==2,  title({'evoked stop',''}), end
        subplot(4,4,4+ 4*(lyr-2))
        plot(evoked(use,2,1),evoked(use,2,2),symb); hold on; plot([-30 30],[-30 30]);  axis square; axis([-10 20 -10 20])
        if lyr==2,  title({'evoked move',''}), end
    end
    
end



figure
for lyr = 2:5
    use = find(goodAll & layerAll==lyr & ~inhAll & (squeeze(evoked(:,2,1)>1)' & squeeze(evoked(:,2,2)>1)') ); symb = 'bo';
    evokedbar(lyr-1,:,:) = median(evoked(use,:,:)); evoked_err(lyr-1,:,:) = std(evoked(use,:,:))/sqrt(length(use));
end

figure
for lyr = 2:5
    use = find(goodAll & layerAll==lyr & ~inhAll  ); symb = 'bo';
    spontbar(lyr-1,:,:) = median(spont(use,:,:)); spont_err(lyr-1,:,:) = std(spont(use,:,:))/sqrt(length(use));
end



figure
for lyr = 2:5
    for i=1:2
        if i==1
            use = find(goodAll & layerAll==lyr & ~inhAll ); symb = 'bo';
        else
            use = find(goodAll & layerAll==lyr  & inhAll ); symb = 'ro';
        end
        subplot(4,4,1 + 4*(lyr-2));
        plot(spont(use,1,1),spont(use,2,1),symb); hold on; plot([0 50],[0 50]);  axis square; axis([0 30 0 30])
        if lyr==2,  title({'spont pre','layer 2'}), end; if lyr>2, title(sprintf('layer %d',lyr)); end
        subplot(4,4,2+ 4*(lyr-2));
        plot(spont(use,1,2),spont(use,2,2),symb);  hold on; plot([0 50],[0 50]); axis square; axis([0 30 0 30])
        if lyr==2,  title({'spont post',''}), end
        subplot(4,4,3+ 4*(lyr-2));
        plot(evoked(use,1,1),evoked(use,2,1),symb); hold on; plot([-30 30],[-30 30]);  axis square; axis([-10 20 -10 20])
        if lyr==2,  title({'evoked pre',''}), end
        subplot(4,4,4+ 4*(lyr-2))
        plot(evoked(use,1,2),evoked(use,2,2),symb); hold on; plot([-30 30],[-30 30]);  axis square; axis([-10 20 -10 20])
        if lyr==2,  title({'evoked post',''}), end
    end
end

figure
barweb(squeeze(spontbar(:,2,:)), squeeze(spont_err(:,2,:)), [],{'layer 2/3','layer 5'},'spontaneous',[],'sp/sec',[],[],{'pre','post'});
ylim([0 1.25])

figure
barweb(squeeze(evokedbar(:,2,:)), squeeze(evoked_err(:,2,:)), [],{'layer 2/3','layer 5'},'evoked',[],'sp/sec',[],[],{'pre','post'});
ylim([0 7])

lyr23 = [squeeze(spontbar(1,2,:)) squeeze(evokedbar(1,2,:))]; lyr23_err = [squeeze(spont_err(1,2,:)) squeeze(evoked_err(1,2,:))];

lyr5 = [squeeze(spontbar(4,2,:)) squeeze(evokedbar(4,2,:))]; lyr5_err = [squeeze(spont_err(4,2,:)) squeeze(evoked_err(4,2,:))];

figure
barweb(lyr23', lyr23_err', [],{'spont','evoked'},'layer 2/3',[],'sp/sec',[],[],{'pre','post'});
ylim([0 7])

figure
barweb(lyr5', lyr5_err', [],{'spont','evoked'},'layer 5',[],'sp/sec',[],[],{'pre','post'});
ylim([0 4])


for i = 1:size(LFPall,1)
    figure
    plot((1:116)/2,squeeze(LFPall(i,1:116,1,1)),'r--'); hold on
    plot((1:116)/2,squeeze(LFPall(i,1:116,2,1)),'g--');
    plot((1:116)/2,squeeze(LFPall(i,1:116,1,2)),'r');
    plot((1:116)/2,squeeze(LFPall(i,1:116,2,2)),'g');
    title(sprintf('wn session %d %s %s,',i,files(useSess(i)).expt,files(useSess(i)).treatment))
end

% for i = 1:size(LFPall,1)
%     figure
%     plot((1:116)/2,squeeze(LFPallDark(i,1:116,1,1)),'r--'); hold on
%     plot((1:116)/2,squeeze(LFPallDark(i,1:116,2,1)),'g--');
%     plot((1:116)/2,squeeze(LFPallDark(i,1:116,1,2)),'r');
%     plot((1:116)/2,squeeze(LFPallDark(i,1:116,2,2)),'g');
%     title(sprintf('dark session %d %s %s,',i,files(useSess(i)).expt,files(useSess(i)).treatment))
% end

for i = 1:size(LFPallCh,1)
    
    d = LFPallCh(i,:,:,:,:);
    range = prctile(d(:),99);
    if range>0
        figure
        subplot(2,2,1); imagesc(squeeze(LFPallCh(i,:,:,1,1)),[0 range]); title('pre stop')
        subplot(2,2,2); imagesc(squeeze(LFPallCh(i,:,:,1,2)),[0 range]);  title('post stop')
        subplot(2,2,3); imagesc(squeeze(LFPallCh(i,:,:,2,1)),[0 range]); title('pre move');
        subplot(2,2,4); imagesc(squeeze(LFPallCh(i,:,:,2,2)),[0 range]);
        
        title(sprintf('dark session %d %s %s,',i,files(useSess(i)).expt,files(useSess(i)).treatment))
    end
end


for lyr = 3:6
    figure
    if lyr<6
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




pre = [cycR(find(goodAll & ~inhAll),:,1,1) cycR(find(goodAll & ~inhAll),:,2,1)];

post = [cycR(find(goodAll & ~inhAll),:,1,2) cycR(find(goodAll & ~inhAll),:,2,2)];

data = [pre post];
[y ind] = sort(layerAll(goodAll & ~inhAll));
data = data(ind,:);
for i = 1:size(data,1);
    normdata(i,:) = data(i,:)/std(data(i,:));
end

% [coeff score latent] = pca(data','variableweight','variance');
% [coeff score latent] = pca(normdata');

[score A w] = fastica(normdata,'numOfIC',3,'lastEig',4);
score = score';

figure
imagesc(A);
figure
imagesc(w');

figure
for i = 1:min(size(score,2),5);
    subplot(5,1,i);
    plot(score(:,i))
end

figure
plot(score(1:20,1),score(1:20,2),'r--'); hold on; plot(score(21:40,1),score(21:40,2),'g--');
plot(score(41:60,1),score(41:60,2),'r'); plot(score(61:80,1),score(61:80,2),'g');

figure
plot(score(1:20,1),score(1:20,3),'r--'); hold on; plot(score(21:40,1),score(21:40,3),'g--');
plot(score(41:60,1),score(41:60,3),'r'); plot(score(61:80,1),score(61:80,3),'g');

figure
plot(score(1:20,2),score(1:20,3),'r--'); hold on; plot(score(21:40,2),score(21:40,3),'g--');
plot(score(41:60,2),score(41:60,3),'r'); plot(score(61:80,2),score(61:80,3),'g');



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


