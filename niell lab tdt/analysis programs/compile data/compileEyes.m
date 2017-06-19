clear all
close all 

dbstop if error

batchEyetracking_angie; %%% load batch file
%set(groot,'defaultFigureVisible','off') %disable figure plotting
set(groot,'defaultFigureVisible','on')

%%% select the session you want based on filters
%use =  find(strcmp({files.notes},'good data'))
use =  find(strcmp({files.notes},'good data') & strcmp({files.expt},'012017b'))
sprintf('%d selected session',length(use))

movieFile = 'C:\Users\Angie Michaiel\Desktop\movie files\cortex\DetectionStim2contrast_LOW_7_25min.mat';
%movieFile = 'C:\Users\Angie Michaiel\Desktop\movie files\cortex\DetectionStim3contrast10min.mat';
%movieFile = 'C:\Users\Angie Michaiel\Desktop\movie files\cortex\flashStim.mat';
load(movieFile);

saline=1; doi=2; ht=3;

for i = 1:length(use)
    
%cfile = {[pathname '\' files(use(i)).dir '\' files(use(i)).predetection_camera '.mat']};
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
        eyes = eye_detection_move(cfile{:,prepost}, prepost, files(use(i)).blockDetect{prepost},files(use(i)).Tank_Name, 1);
        rad{i,:,prepost} = eyes.rad
        rInterp{i,:,prepost}=eyes.rInterp
        vInterp{i,:,prepost}=eyes.vInterp
        cameraT{i,:,prepost} = eyes.t
        frameNum{i,:,prepost} = eyes.frameNum;
        frameT{i,:,prepost} = eyes.frameT;
        fInterpR{i,:,prepost} = eyes.fInterpR ;
        fInterpX{i,:,prepost} = eyes.fInterpX ;
        fInterpY{i,:,prepost} = eyes.fInterpY ; 
        fInterpV{i,:,prepost} = eyes.fInterpV ;
        trials{i,:,:,prepost} = eyes.trials;
        trialV{i,:,prepost} = eyes.trialV;
        trialSamp{i,:,prepost} = eyes.trialSamp;
    end
end   

end
   
% t = length(vInterp{1})./(length(trials{1}))
% vTrials = interp1(vInterp{1}, trials{1},t)


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
%end




%corr of velocity and radius, check lags
% does corr become decoupled post?

% dt = 0.5;
% figure
% for prepost=1:2
% [corr_vrad lags] = cellfun(@xcorr,(~isnan(rInterpHigh{prepost}-(mean(rInterpHigh{prepost})))),(~isnan(vInterpHigh-(mean(vInterpHigh)))),60/dt,'coeff')
% [corr_vrad lags] = cellfun(@xcorr,(~isnan(rInterp{prepost}-(mean(rInterp{prepost})))),(~isnan(vInterp-(mean(vInterp)))),60/dt,'coeff')
% 


% subplot(1,2,prepost)
% plot(lags*dt,corr_vrad);
% % hold on
% % plot([0 0],[0 1],'g-')
% end

% trialSampPre = [trialSamp{1}];% trialSampPre = trialSampPre(1:length(contrast))
% trialSampPost = [trialSamp{2}]; %trialSampPost = trialSampPost(1:length(contrast))
% try
% trialContPre = [trialSampPre;contrast(1:length(trialSampPre))]'
% catch
% trialContPre = [trialSampPre(1:length(contrast));contrast;]'
% end
% try
% trialContPost = [trialSampPost;contrast(1:length(trialSampPost))]'
% catch
% trialContPost = [trialSampPost(1:length(contrast));contrast]'
% end
% % 
%  use = unique(contrast)
%  useX = unique(xpos)
% figure
% for i=1:length(use)
% %     subplot(2,2,1)
% %     plot(trialContPre(contrast==0)); hold on;plot(trialContPost(contrast==0))
%     subplot(2,2,i)
%     plot(trialContPre(contrast(1:length(trialContPre))==use(i))); hold on; plot(trialContPost(contrast(1:length(trialContPre))==use(i)))
% end
% 
% figure
% Labels = {'0', '0.01', '0.04', '1.0'};
% for i=1:length(use)
%     subplot(1,4,i)
%     premean = nanmean(trialContPre(contrast(1:length(trialContPre))==use(i))); postmean = nanmean(trialContPost(contrast(1:length(trialContPre))==use(i)))
%     preEr = nanstd(trialContPre(contrast(1:length(trialContPre))==use(i)))/sqrt(length(trialContPre(contrast(1:length(trialContPre))==use(i))))
%     postEr = nanstd(trialContPost(contrast(1:length(trialContPre))==use(i)))/sqrt(length(trialContPost(contrast(1:length(trialContPre))==use(i))))
%     mn = [premean postmean]; err = [preEr postEr];
%     barweb(mn,err);axis square
%     set(gca, 'XTick', 1:4, 'XTickLabel', Labels(i));
% end
% legend('pre', 'post')
% 
% % 
% trials =squeeze(trials)'
% preTrials = [trials{1}]; preTrials = preTrials(1:length(contrast))
% postTrials = [trials{2}]
% preV = [trialV{1}]
% postV = [trialV{2}]
% figure; subplot(1,2,1)
% imagesc(preTrials); 
% subplot(1,2,2); 
% imagesc(postTrials)
% figure
% subplot(1,2,1)
% imagesc(preV);
% subplot(1,2,2)
% imagesc(postV)
% 
% figure
% plot(fInterpR{1});hold on;plot(fInterpV{1})
% % 
% % % % %need to normalize
% titles = {'0', '0.01', '0.04'};
% figure
% for i= 1:length(use)
% subplot(1,3,i);
% imagesc(preTrials(contrast==use(i)));
% axis square
% title(titles{i});
% % subplot(2,2,i+2)
% % imagesc(postTrials(contrast==use(i),:));axis square
% end
% % 
% titles = {'x pos L','x pos R'};
% figure
% for i= 1:length(useX)
% subplot(1,2,i);
% plot(nanmean(preTrials(xpos==useX(i)& contrast==.04)));
% % xlim([0 120]); %ylim([15 25]); 
% axis square;hold on;
% % plot([60 60],[14 25],'g');
% title(titles{i})
% end
% % 
% clear use
% use = unique(contrast);
% titles = {'0','0.01','0.04'};
% figure
% for i= 1:length(use)
% subplot(1,3,i);
% plot(mean(preTrials(contrast==use(i),:)));
% % xlim([0 120]);
% % axis square; hold on
% % % ylim([-.6 1])
% % %title(titles{i});
% % subplot(1,3,i)
% % plot([60 60],[16 22],'g');
% % %plot(nanmean(postTrials(contrast==use(i),:)));
% % %line([60 0],'g')
% % ylim([16 22]);
% % axis square
% end
% 
% figure
% c=find(contrast==1)
% for d=1:16
% subplot(4,4,d)
% plot(preTrials(c(d),:));
% hold on
% plot(preV(c(d),:),'r')
% xlim([0 130]);axis square;
% %ylim([-4 7])
% %title(titles{i});
% % subplot(1,4,i)
% % plot(nanmean(postTrials(contrast==use(i),:)));
% ylim([-2 40])
% plot([60 60],[-2 40],'g');
% plot ([75 75],[-2 40],'g');
% axis square
% end
% 
% clear c
% for i =1:2
%     c=find(contrast==1 & xpos==useX(i))
%     figure
%     if i ==1
%         set(gcf,'Name','trials L position')
%     else
%         set(gcf,'Name','trials R position')
%     end
%     for d=1:8
%         subplot(2,4,d);
%         plot(preTrials(c(d),:)); hold on;
%         plot(preV(c(d),:));
%         xlim([0 130]);axis square;
%         ylim([-2 50])
%         plot([60 60],[-2 50],'g');
%         plot ([75 75],[-2 50],'g');
%        % legend('radius','velocity')
%     end
% end
% 
% figure
% %subplot(1,2,1)
% plot(nanmean(preTrials(contrast(1:length(preTrials))==0,:)));
% hold on; plot(nanmean(preTrials(contrast(1:length(preTrials))==1,:)));
% plot(nanmean(preV(contrast(1:length(preTrials))==1,:)));
% plot([60 60],[0 30],'g');  plot ([75 75],[0 30],'g');
% legend('0 contrast', '1 contrast','velocity');
% 


figure
plot(rad{1,1,1}); hold on; plot(rad{1,1,2}); ylim([8 30]);
xlim([0 3000]);set(gca,'xtick',600:600:7.25*600,'xticklabel',1:5,'Fontsize',14)
xlabel ('time (min)');

prerad= rad{1}; postrad = rad{2};
figure
h1=hist(prerad(1:1500),1:1:max(prerad))
h2 = hist(postrad(1:1500),1:1:max(prerad))
plot(h1/1500,'Linewidth',2); hold on; plot(h2/1500,'r','Linewidth',2); xlim([10 25]);ylim([0 .5])



