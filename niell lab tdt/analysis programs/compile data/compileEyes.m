clear all
close all 

dbstop if error

batchEyetracking_angie; %%% load batch file
%set(groot,'defaultFigureVisible','off') %disable figure plotting
set(groot,'defaultFigureVisible','on')

%%% select the session you want based on filters
%use =  find(strcmp({files.notes},'good data'))
 use =  find(strcmp({files.notes},'good data') & strcmp({files.expt},'011817'))
sprintf('%d selected session',length(use))

movieFile = 'C:\Users\Angie Michaiel\Desktop\movie files\cortex\DetectionStim3contrast10min.mat';
%movieFile = 'C:\Users\Angie Michaiel\Desktop\movie files\cortex\flashStim.mat';
load(movieFile);

saline=1; doi=2; ht=3;

for i = 1:length(use)
    
    cfile = [{[pathname '\' files(use(i)).dir '\' files(use(i)).predetection_camera '.mat']};
        {[pathname '\' files(use(i)).dir '\' files(use(i)).postdetection_camera '.mat']}]';

% session = use(i)
% 
%    sessionNum(use)=i;
%     for j=1:length(session); expt{session(j)} = files(use(i)).expt; end
%     if strcmp(files(use(i)).treatment,'Saline'), treatment(session)=saline, end;
%     if strcmp(files(use(i)).treatment,'DOI'), treatment(session)=doi, end;
%     if strcmp(files(use(i)).treatment,'5HT'), treatment(session)=ht, end;
%     sessionTreatment(i) = treatment(sessionNum(i));

if ~isempty(files(use(i)).blockDetect{1}) & ~isempty(files(use(i)).blockDetect{2})
    
    for prepost =1:2
        eyes = eye_detection_move(cfile{:,prepost}, prepost, files(use(i)).blockDetect{prepost},files(use(i)).Tank_Name, 0);
        rad{i,:,prepost} = eyes.rad
        rInterp{i,:,prepost}=eyes.rInterp
        vInterp{i,:,prepost}=eyes.vInterp
        cameraT{i,:,prepost} = eyes.t
        frameNum{i,:,prepost} = eyes.frameNum;
        frameT{i,:,prepost} = eyes.frameT;
        fInterpR{i,:,prepost} = eyes.fInterpR ;
        fInterpX{i,:,prepost} = eyes.fInterpX ;
        fInterpY{i,:,prepost} = eyes.fInterpY ;
        trials{i,:,:,prepost} = eyes.trials;
        trialSamp{i,:,prepost} = eyes.trialSamp;
    end
end   
% 
% 
% if ~isempty(files(use(i)).blockHigh{1}) & ~isempty(files(use(i)).blockHigh{2})
%     
%     for prepost =1:2
%         eyes = getEyes_angie(cfile{:,prepost}, files(use(i)).blockHigh{prepost},1);
%         radHigh{:,prepost} = eyes.rad
%         rInterpHigh{:,prepost}=eyes.rInterp
%         vInterpHigh{:,prepost}=eyes.vInterp
%         tHigh{:,prepost} = eyes.t
%     end
% end   
% 
%    if ~isempty(files(use(i)).blockLow{1}) & ~isempty(files(use(i)).blockLow{2})
% 
%        %low doesn't work---all empty
%    for prepost =1:2
%         eyes = getEyes_angie(cfile{:,prepost}, files(use(i)).blockLow{prepost},1);
%         radLow{:,prepost} = eyes.rad
%         rInterpLow{:,prepost}=eyes.rInterp
%         vInterpLow{:,prepost}=eyes.vInterp
%         tLow{:,prepost} = eyes.t
%    end
% end
end

figure
plot(frameT{1},fInterpR{1});hold on; plot(frameT{2},fInterpR{2},'r')



% dt = 0.5;
% figure
% for prepost=1:2
% %[corr_vrad lags] = xcorr((~isnan(rInterpHigh{prepost}-nanmean(rInterpHigh{prepost})),(~isnanvInterpHigh{prepost}-nanmean(vInterpHigh{prepost}))),60/dt,'coeff');
% [corr_vrad lags] = cellfun(@xcorr,(~isnan(rInterpHigh{prepost}-(mean(rInterpHigh{prepost})))),(~isnan(vInterpHigh-(mean(vInterpHigh)))),60/dt,'coeff')
% 
% subplot(1,2,prepost)
% plot(lags*dt,corr_vrad);
% % hold on
% % plot([0 0],[0 1],'g-')
% end

trialSampPre = [trialSamp{1}];% trialSampPre = trialSampPre(1:length(contrast))
trialSampPost = [trialSamp{2}]; %trialSampPost = trialSampPost(1:length(contrast))
try
trialContPre = [trialSampPre;contrast(1:length(trialSampPre))]'
catch
trialContPre = [trialSampPre(1:length(contrast));contrast;]'
end
trialContPost = [trialSampPost;contrast(1:length(trialSampPost))]'

use = unique(contrast)
figure
for i=1:length(use)
%     subplot(2,2,1)
%     plot(trialContPre(contrast==0)); hold on;plot(trialContPost(contrast==0))
    subplot(2,2,i)
    plot(trialContPre(contrast(1:length(trialContPre))==use(i))); hold on; plot(trialContPost(contrast(1:length(trialContPre))==use(i)))
end

figure
Labels = {'0', '0.01', '0.04', '1.0'};
for i=1:length(use)
    subplot(1,4,i)
    premean = nanmean(trialContPre(contrast(1:length(trialContPre))==use(i))); postmean = nanmean(trialContPost(contrast(1:length(trialContPre))==use(i)))
    preEr = nanstd(trialContPre(contrast(1:length(trialContPre))==use(i)))/sqrt(length(trialContPre(contrast(1:length(trialContPre))==use(i))))
    postEr = nanstd(trialContPost(contrast(1:length(trialContPre))==use(i)))/sqrt(length(trialContPost(contrast(1:length(trialContPre))==use(i))))
    mn = [premean postmean]; err = [preEr postEr];
    barweb(mn,err);axis square
    set(gca, 'XTick', 1:4, 'XTickLabel', Labels(i));
end
legend('pre', 'post')


trials =squeeze(trials)'
preTrials = [trials{1}]
postTrials = [trials{2}]
figure; subplot(1,2,1)
imagesc(preTrials); subplot(1,2,2); imagesc(postTrials)

%need to normalize
titles = {'0', '0.01', '0.04', '1.0'};
figure
for i= 1:length(use)
subplot(2,4,i);imagesc(preTrials(contrast(1:length(preTrials))==use(i),:));axis square
title(titles{i});
subplot(2,4,i+4)
imagesc(postTrials(contrast==use(i),:));axis square
end

titles = {'0', '0.01', '0.04', '1.0'};
figure
for i= 1:length(use)
subplot(1,4,i);plot(nanmean(preTrials(contrast(1:length(preTrials))==use(i),:)));
xlim([0 80]);ylim([-.6 1]);axis square; hold on
title(titles{i});
subplot(1,4,i)
plot(nanmean(postTrials(contrast==use(i),:)));
xlim([0 80]);ylim([-.6 1]);axis square
end

% figure
% plot(nanmean(preTrials(contrast(1:length(preTrials))==0,:)));
% hold on; plot(nanmean(preTrials(contrast(1:length(preTrials))==1,:)));
% xlim([0 80]);ylim([-.6 1]);axis square; hold on


