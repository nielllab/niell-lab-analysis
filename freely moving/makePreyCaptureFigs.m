close all;
clear all;
%    load('ACCAnalyzed_AllAnimals_010820_noDLS.mat')
% load('ACC_AllAnimals_021520_a.mat')
% load('DEINTER_Analyzed_AllAnimals_051820.mat');
%load('ACC_deInter_Analyzed_AllAnimals_011520_a.mat')

%%% NOTE: find and replace entire script deInter==0 for original dataset or ==1 for
%%% deinterlaced - it will load the appropriate dataset and will adjust framerate

deInter=1;
if deInter
    frRate=60;
    %     load('DEINTER_Analyzed_AllAnimals_051820.mat');
%     load('DEINTERLACED_Analyzed_AllAnimals_052020.mat')
%    load('DEINTERLACED_Analyzed_AllAnimals_052320_halfShift.mat')
   load('DEINTERLACED_Analyzed_AllAnimals_060220_halfShift.mat')
else
    frRate=30;
    load('ACC_AllAnimals_021520_a.mat')
end

savePDF=1;
if savePDF
    psfilename = 'C:\analysisPS_B.ps';
    if exist(psfilename,'file')==2;delete(psfilename);end
end
pname = 'T:\PreyCaptureAnalysis\Data\';

%% FIGURE 1: a freely moving eye and head tracking system
load('T:\PreyCaptureAnalysis\Data\ControlAnalysis\compileAllAnimals_CONTROL_090219.mat','mouseSp','appEpoch','thetaHead');

deInter=1;
if deInter
    frRate=60;
else
    frRate=30;
end


for vid=1:length(thetaHead)
    prop(vid)=sum(~isnan(thetaHead{vid}))./length(thetaHead{vid});
end
useData=find(prop>.85); % 85% of points needs to be present for exp to be used

for vid=1:length(useData)
    appTime=appEpoch{useData(vid)};
    %      w=gausswin(15);
    %     runningSmooth=filter(w,1,mouseSp{useData(vid)});
    runningSmooth =medfilt1(mouseSp{useData(vid)},(frRate/2));
    Cntrl_speed(vid,1) =sum(runningSmooth(:,appTime==0)>5)./length(runningSmooth(:,appTime==0));
    
    Cntrl_speed(vid,3) =nanmean(runningSmooth(:,appTime==0));
    
    if sum(appTime==1)>frRate
        Cntrl_speed(vid,2) =sum(runningSmooth(:,appTime==1)>5)./length(runningSmooth(:,appTime==1));
        Cntrl_speed(vid,4) =nanmean(runningSmooth(:,appTime==1));
        
    else
        Cntrl_speed(vid,2)=nan;
        Cntrl_speed(vid,4)=nan;
        
    end
end


clear mouseSp appEpoch useData appTime thetaHead vid runningSmooth
% load('ACCAnalyzed_AllAnimals_010820_noDLS.mat')
% load('ACC_AllAnimals_021520_a.mat')
% load('DEINTER_Analyzed_AllAnimals_051820.mat');
deInter=1;
if deInter
    frRate=60;
    %     load('DEINTER_Analyzed_AllAnimals_051820.mat');
%     load('DEINTERLACED_Analyzed_AllAnimals_052020.mat')
%   load('DEINTERLACED_Analyzed_AllAnimals_052320_halfShift.mat')
       load('DEINTERLACED_Analyzed_AllAnimals_060220_halfShift.mat')

else
    frRate=30;
    load('ACC_AllAnimals_021520_a.mat')
end


savePDF=1;
psfilename = 'C:\analysisPS_B.ps';
for vid=1:length(useData)
    appTime=appEpoch{vid};
    w=gausswin(frRate/2);
    %     runningSmooth=filter(w,1,mouseSp{useData(vid)});
    runningSmooth =medfilt1(mouseSp{useData(vid)},frRate/2);
    
    Cam_speed(vid,1) = sum(runningSmooth(:,appTime==0)>5)./(length(runningSmooth(:,appTime==0)));
    
    
    Cam_speed(vid,3) = nanmean(runningSmooth(:,appTime==0));
    if sum(appTime==1)>frRate
        Cam_speed(vid,2) =sum(runningSmooth(:,appTime==1)>5)./(length(runningSmooth(:,appTime==1)));
        Cam_speed(vid,4) = nanmean(runningSmooth(:,appTime==1));
    else
        Cam_speed(vid,2)=nan;
        Cam_speed(vid,4)=nan;
    end
end

Cntrl_speed(Cntrl_speed(:,4)>60,4)=NaN; % noisy control data, speeds at a few timepoints go above 150 cm/sec - bad DLC tracking
Cntrl_err= nanstd(Cntrl_speed)/sqrt(length(useData));
Cam_err= nanstd(Cam_speed)/sqrt(length(useData));
% figure; barweb([nanmean(Cam_speed(:,1)) nanmean(Cam_speed(:,2)) nanmean(Cntrl_speed(:,1)) nanmean(Cntrl_speed(:,2))],...
%     [nanmean(Cam_err(:,1)) nanmean(Cam_err(:,2)) nanmean(Cntrl_err(:,1)) nanmean(Cntrl_err(:,2))]);
% title('proportion of time moving')
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


[h p]= ttest2(Cntrl_speed(:,3),Cam_speed(:,3))
[h p]= ttest2(Cntrl_speed(:,4),Cam_speed(:,4))
%
% [h,p,ci,stats] = ttest2(Cntrl_speed(:,3),Cam_speed(:,3), 'both','Alpha',.05)
% figure; barweb([nanmean(Cam_speed(:,3)) nanmean(Cam_speed(:,4)) nanmean(Cntrl_speed(:,3)) nanmean(Cntrl_speed(:,4))],...
%     [nanmean(Cam_err(:,3)) nanmean(Cam_err(:,4)) nanmean(Cntrl_err(:,3)) nanmean(Cntrl_err(:,4))]);
% title('mean speed')

figure;plot(2,Cam_speed(:,3),'ko'); hold on
plot(2,nanmean(Cam_speed(:,3)),'b*','Markersize',15);
plot(2,nanmean(Cam_speed(:,3))+(Cam_err(:,3)),'r*','Markersize',10);
plot(2,nanmean(Cam_speed(:,3))-(Cam_err(:,3)),'r*','Markersize',10);
plot(3,Cam_speed(:,4),'go')
plot(3,nanmean(Cam_speed(:,4)),'b*','Markersize',15)
plot(3,nanmean(Cam_speed(:,4))+(Cam_err(:,4)),'r*','Markersize',10);
plot(3,nanmean(Cam_speed(:,4))-(Cam_err(:,4)),'r*','Markersize',10);

plot(5,Cntrl_speed(:,3),'ko')
plot(5,nanmean(Cntrl_speed(:,3)),'b*','Markersize',15)
plot(5,nanmean(Cntrl_speed(:,3))+(Cntrl_err(:,3)),'r*','Markersize',10);
plot(5,nanmean(Cntrl_speed(:,3))-(Cntrl_err(:,3)),'r*','Markersize',10);

plot(6,Cntrl_speed(:,4),'go')
plot(6,nanmean(Cntrl_speed(:,4)),'b*','Markersize',15)
plot(6,nanmean(Cntrl_speed(:,4))+(Cntrl_err(:,4)),'r*','Markersize',10);
plot(6,nanmean(Cntrl_speed(:,4))-(Cntrl_err(:,4)),'r*','Markersize',10);
xlim([1 7]); axis square
legend('cam','cam app','control','control app');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
[h p]= kstest2(Cntrl_speed(:,3),Cntrl_speed(:,4))
[h p]= kstest2(Cam_speed(:,3),Cam_speed(:,4))
[h p]= ttest2(Cam_speed(:,3),Cntrl_speed(:,3))
[h p]= ttest2(Cam_speed(:,4),Cntrl_speed(:,4))

%%%average number of captures per session %avg numbers of crickets for each animal, across sessions, from experiment spreadsheet
exp=[4.6	6.2	5.2	4.333333333	3.166666667	4.166666667	5.2	6.6	3.5	4.166666667];
cntrl=[4.933333333	6.366666667	5.366666667	5.666666667	5.833333333	4.5	5	6.6	5.166666667	5.333333333];
expM=mean(exp); errExp=std(exp)./(length(sqrt(exp)));
mnC=mean(cntrl); errC=std(cntrl)./(length(sqrt(cntrl)));
comp=[expM mnC];err=[errExp,errC];
figure
barweb(comp,err)
[sig p]=ttest2(exp',cntrl');

%%% old ex trace
% Figure 1G: gyro & DLC traces

for vid=105
    
    use=appEpoch{vid};
    if deInter
        range=2*1148:2*1655
        tR=Rtheta{useData(vid)}(:,range)-nanmedian(Rtheta{useData(vid)}(:,range));
        tL=Ltheta{useData(vid)}(:,range)-nanmedian(Ltheta{useData(vid)}(:,range));
        tR=interpNan(tR,10,'linear'); tL=interpNan(tL,10,'linear');
        mnEye=.5*(tR+tL); mnEye=interpNan(mnEye,10,'linear');
        headTh=thetaHead{useData(vid)}(:,range)-nanmedian(thetaHead{useData(vid)}(:,range));
    else
        range =1148:1655
        tR=Rtheta{useData(vid)}(range,:)-nanmedian(Rtheta{useData(vid)}(range,:));
        tL=Ltheta{useData(vid)}(range,:)-nanmedian(Ltheta{useData(vid)}(range,:))
        mnEye=.5*(tR+tL);
        
    end
    headTh=thetaHead{useData(vid)}(:,range)-nanmedian(thetaHead{useData(vid)}(:,range));
    eyeVelocity=interpNan(diff(mnEye), 10,'linear');
    
    figure
    subplot(7,1,1)
    plot(tR); hold on
    plot(tL);
    plot(mnEye-nanmedian(mnEye));xlim([1 length(range)])
    ylim([-32 32]);
    subplot(7,1,2)
    plot(eyeVelocity);
    xlim([1 length(range)]);
    title('eye velocity')
    subplot(7,1,3)
    plot(dist{useData(vid)}(:,range));xlim([1 length(range)]); hold on
    ylim([0 40])
    plot([296 296],[-60 60],'g')
    title('distance')
    subplot(7,1,4)
    plot(rad2deg(az{useData(vid)}(:,range)));xlim([1 length(range)])
    title('azimuth')
    subplot(7,1,5);
    plot(mouseSp{useData(vid)}(:,range))%-nanmedian(mouseSp{useData(vid)}(:,range)));
    title('speed')
    xlim([1 length(range)])
    subplot(7,1,6)
    plot(headTh);
    % plot(cumsum(accelChannels{useData(vid)}(range,6))-nanmedian(cumsum(accelChannels{useData(vid)}(range,6))));
    xlim([1 length(range)])
    subplot(7,1,7)
    plot(accelChannels{useData(vid)}(range,2)-nanmedian(accelChannels{useData(vid)}(range,2)));
    xlim([1 length(range)])
end
% %%%RMS err for dTheta measures
% for vid=1:length(useData)
% 
% comp(vid,:)=nanmean(d_Theta{useData(vid)}-accelChannels{useData(vid)}(:,6));
% end
% nanstd(comp)

% figure
% range=350:830;
% for vid=10%1:10%length(useData)
%     %     subplot(2,5,vid);
%     plot(d_Theta{useData(vid)}(range,:));
%     hold on;
%     plot(accelChannels{useData(vid)}(range,6));
% end
% legend('DLC head theta','gyroscope')
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

% Figure 1H: measure of similarity of gyro & DLC (scatter or hist of
% difference)


figure
clear h
subplot(1,2,1)
plot(gyro3All(1:100:end),dlcDhth(1:100:end),'.'); xlabel('gyro yaw'); ylabel('DLC yaw')
axis equal; xlim([-20 20]); ylim([-20 20]);
bins=-20:2:20
h=hist(gyro3All-dlcDhth,bins)
subplot(1,2,2)
plot(bins,h/length(gyro3All)); hold on
xlabel('gyro - DLC (degrees');
suptitle('Figure 1H: DLC is good measure of head angle')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

% Figure 1I (where is the right place to put this?): head angle is directed
% towards cricket - corr of az and change in head yaw


%% FIGURE 2: Coordination of eyes during free movement

% Figure 2A: scatter plots of R & L eyes yaw vs pitch, approach and non-approach
%%
figure
for vid=1:length(useData)
    nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
    pR=Rphi{useData(vid)}(1:nframe); pL=Lphi{useData(vid)}(1:nframe);
    pR=nanmean(pR)-pR; pL=nanmean(pL)-pL;
    tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
    tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
    clear use
    useN = find(appEpoch{vid}==0); use = find(appEpoch{vid}==1);
    ptsUsed(vid,1,1)=length(1:480:length(useN));
    ptsUsed(vid,1,2)=length(1:480:length(use));
    ptsUsed(vid,2,1)=length(useN);
    ptsUsed(vid,2,2)=length(use);
    subplot(1,2,1)
    plot(tR(useN(1:480:length(useN))),pR(useN(1:480:length(useN))),'.b'); hold on; xlim([-50 50]);ylim([-50 50]); axis square
    plot(tR(use(1:480:length(use))),pR(use(1:480:length(use))),'.g'); hold on; xlim([-50 50]);ylim([-50 50]); axis square
    title('R eye')
    xlabel('yaw (deg)'); ylabel('pitch(deg)');
    subplot(1,2,2)
    plot(tL(useN(1:480:length(useN))),pL(useN(1:480:length(useN))),'.b'); hold on; xlim([-50 50]);ylim([-50 50]); axis square
    plot(tL(use(1:480:length(use))),pL(use(1:480:length(use))),'.g'); hold on; xlim([-50 50]);ylim([-50 50]); axis square
    xlabel('yaw (deg)'); ylabel('pitch(deg)');
    title('L eye')
end

suptitle('Figure 2A: R & L eye yaw & pitch');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

nansum(ptsUsed(:,1,1))/nansum(ptsUsed(:,2,1))
nansum(ptsUsed(:,1,2))/nansum(ptsUsed(:,2,2))

% figure
% for vid=1:length(useData)
%     nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
%     pR=Rphi{useData(vid)}(1:nframe); pL=Lphi{useData(vid)}(1:nframe);
%     pR=nanmean(pR)-pR; pL=nanmean(pL)-pL;
%     tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
%     tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
%         clear use
%         useN = appEpoch{vid}==0;use = appEpoch{vid}==1;
%         subplot(1,2,1)
%         polarscatter(tR(useN(1:100:length(useN))),pR(useN(1:100:length(useN))),'.b'); hold on;% xlim([-100 100]);ylim([-100 100]); axis square
%         polarscatter(tR(use(1:50:length(use))),pR(use(1:50:length(use))),'.g'); hold on; rlim([0 40]);%ylim([-100 100]); axis square
%         title('R eye')
%      %   xlabel('yaw (deg)'); ylabel('pitch(deg)');
%         subplot(1,2,2)
%         polarscatter(tL(useN(1:100:length(useN))),pL(useN(1:100:length(useN))),'.b'); hold on; rlim([0 40])% xlim([-100 100]);ylim([-100 100]); axis square
%         polarscatter(tL(use(1:50:length(use))),pL(use(1:50:length(use))),'.g'); hold on;% xlim([-100 100]);ylim([-100 100]); axis square
%       %  xlabel('yaw (deg)'); ylabel('pitch(deg)');
%         title('L eye')
% end
% suptitle('Figure 2A: R & L eye yaw & pitch');
%%
% Figure 2B: histograms of yaw and pitch from panel above
clear tNcounts tAcounts allThN allThA
for vid=1:length(useData)
    nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
    pR=Rphi{useData(vid)}(1:nframe); pL=Lphi{useData(vid)}(1:nframe);
    pR=nanmean(pR)-pR; pL=nanmean(pL)-pL;
    tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
    tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
    useN= appEpoch{vid}==0; use= appEpoch{vid}==1;
    clear r l
    if sum(useN)>4
        nBins=-90:5:90; nBinsP=-50:5:50
        r= (hist(tR(useN),nBins))/sum(useN);l= (hist(tL(useN),nBins))/sum(useN);
        rP= (hist(pR(useN),nBinsP))/sum(useN);lP= (hist(pL(useN),nBinsP))/sum(useN);
    else
        r=NaN;l=NaN;
        rP=NaN;lP=NaN;
    end
    tNcounts(vid,:,1)=r;
    tNcounts(vid,:,2)=l;
    tNcountsP(vid,:,1)=rP;
    tNcountsP(vid,:,2)=lP;
    
end

for vid=1:length(useData)
    nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
    tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
    tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
    pR=Rphi{useData(vid)}(1:nframe); pL=Lphi{useData(vid)}(1:nframe);
    pR=nanmean(pR)-pR; pL=nanmean(pL)-pL;
    use = (appEpoch{vid});
    %     plot(tR(use(1:20:end)),tL(use(1:20:end)),'go'); hold on; axis square; %ylim([-frRate frRate]); xlim([-frRate frRate]);
    clear r l
    if sum(use)>4
        nBins=-90:5:90; nBinsP=-50:5:50
        r= (hist(tR(use),nBins))/sum(use);l= (hist(tL(use),nBins))/sum(use);
        rP= (hist(pR(use),nBinsP))/sum(use);lP= (hist(pL(use),nBinsP))/sum(use);
        
    else
        r=NaN;l=NaN;
        rP=NaN;lP=NaN;
    end
    tAcounts(vid,:,1)=r;
    tAcounts(vid,:,2)=l;
    tAcountsP(vid,:,1)=rP;
    tAcountsP(vid,:,2)=lP;
end


allThN(:,1)=nansum(tNcounts(:,:,1),1)./(length(useData))
allThN(:,2)=nansum(tNcounts(:,:,2),1)./(length(useData))

allThA(:,1)=nansum(tAcounts(:,:,1),1)./(length(useData))
allThA(:,2)=nansum(tAcounts(:,:,2),1)./(length(useData))

mns=squeeze(nanmean(tNcounts,2));
mnsA=squeeze(nanmean(tAcounts,2));

[h p]=ttest2(mns(:,1),mnsA(:,1))%R eye
[h p]=ttest2(mns(:,2),mnsA(:,2))%L eye

figure;subplot(1,2,1); plot(nBins,allThN(:,1),'b'); hold on; plot(nBins,allThA(:,1),'g'); title('Right Eye Theta'); axis square
ylim([0 .27]);xlim([-50 50]); xlabel('theta (degrees)'); ylabel('proportion of time');
subplot(1,2,2);plot(nBins,allThN(:,2),'b'); hold on; plot(nBins,allThA(:,2),'g');title('Left Eye Theta');axis square
ylim([0 .27]);xlim([-50 50]); xlabel('theta (degrees)'); ylabel('proportion of time');
legend('non approach','approach')
suptitle('Figure 2B: eyes are more centered in yaw during app')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


allThNP(:,1)=nansum(tNcountsP(:,:,1),1)./(length(useData))
allThNP(:,2)=nansum(tNcountsP(:,:,2),1)./(length(useData))

allThAP(:,1)=nansum(tAcountsP(:,:,1),1)./(length(useData))
allThAP(:,2)=nansum(tAcountsP(:,:,2),1)./(length(useData))

mnsPhi=squeeze(nanmean(tNcountsP,2));
mnsAPhi=squeeze(nanmean(tAcountsP,2));

[h p]=ttest2(mnsPhi(:,1),mnsAPhi(:,1))%R eye
[h p]=ttest2(mnsPhi(:,2),mnsAPhi(:,2))%L eye


figure;subplot(1,2,1); plot(nBinsP,allThNP(:,1),'b'); hold on; plot(nBinsP,allThAP(:,1),'g'); title('Right Eye Phi'); axis square
ylim([0 .35]);
subplot(1,2,2);plot(nBinsP,allThNP(:,2),'b'); hold on; plot(nBinsP,allThAP(:,2),'g');title('Left Eye Phi');axis square
ylim([0 .35]);
xlabel('phi (degrees)'); ylabel('proportion of time');
legend('non approach','approach')
suptitle('Figure 2B: eyes are more centered in phi during app')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%
% Figure 2C: overlaid trace of two eyes converging and diverging
clear all
load('ACC_deInter_Analyzed_AllAnimals_011520_a.mat') %original preprint
% used first deinterlaced dataset for example - video 40
% load('DEINTERLACED_Analyzed_AllAnimals_052020.mat')
% load('DEINTERLACED_Analyzed_AllAnimals_052320_halfShift.mat')


deInter=1;
if deInter
    frRate=60;
else
    frRate=30;
end
savePDF=1;
psfilename = 'C:\analysisPS_B.ps';
% 
% 
% range=295:1495; %20 seconds
% figure
% appTimes=find(appEpoch{24}(:,range)==1)
% rT=Rtheta{useData(24)}(:,range);
% lT=Ltheta{useData(24)}(:,range)
% plot(rT-nanmean(rT)); hold on;
% plot(lT-nanmean(lT)); hold on;
% xlim([1 length(range)])
% % ylim([-32 32]);
% legend('r theta','l theta');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% 
% figure
% w = gausswin(10);
% y = filter(w,1,mouseSp{useData(24)}(:,range))
% plot(y/8); hold on;
% xlim([1 length(range)])
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


% figure
% appTimes=find(appEpoch{40}(1:5400)==1)
% plot(Rtheta{useData(40)}(1:5400)); hold on; %3 minute segment
% plot(Ltheta{useData(40)}(1:5400)); hold on;
% % plot(appTimes,20,'og')
% legend('r theta','l theta');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
%
% range=3000:4800
% figure
% plot(Rtheta{useData(40)}(range,:)); hold on; % 1 minute segment
% plot(Ltheta{useData(40)}(range,:)); hold on;
% legend('r theta','l theta');
% % plot(accelChannels{useData(40)}(range,1)); %roll
% % plot(Rtheta{useData(40)}(range,:)-Ltheta{useData(40)}(range,:))
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% range=3600:4400
% figure
% plot(Rtheta{useData(40)}(range,:)); hold on; % 30 second segment
% plot(Ltheta{useData(40)}(range,:)); hold on;
% legend('r theta','l theta');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %
% figure
% w = gausswin(10);
% y = filter(w,1,mouseSp{useData(40)}(:,range))
% plot(y); hold on; % 30 second segment
%
%
% close all
% range=4000:4600
% % range=2000:2300
% for vid = 40%(useData)
% 
%    if length(Rtheta{useData(vid)})>=100
% % tr=medfilt1(Rtheta{useData(vid)}(range,:),10);
% % pr=medfilt1(Rphi{useData(vid)}(range,:),10);
% %
% % tl=medfilt1(Ltheta{useData(vid)}(range,:),10);
% % pl=medfilt1(Lphi{useData(vid)}(range,:),10);
% 
% tr=Rtheta{useData(vid)}(range,:); tr=tr-nanmean(tr);
% pr=Rphi{useData(vid)}(range,:); pr=pr-nanmean(pr);
% 
% tl=(Ltheta{useData(vid)}(range,:));tl=tl-nanmean(tl);
% pl=(Lphi{useData(vid)}(range,:)); pl=pl-nanmean(pl);
% 
% filtwin=gausswin(10);
% tr=filter(filtwin,5,tr);
% pr=filter(filtwin,5,pr);
% tl=filter(filtwin,5,tl);
% pl=filter(filtwin,5,pl);
% 
% figure(1)
% subplot(1,2,1)
% plot(tr-nanmean(tr),pr-nanmean(pr),'k'); axis square; hold on; xlabel('eye theta'); ylabel('eye phi'); axis equal
% title('right eye'); xlim([-23 23]); ylim([-23 23]); colormap jet; colorbar
% figure(2)
% subplot(1,2,2)
% plot(tl-nanmean(tl),pl-nanmean(pl),'k'); axis square; hold on; xlabel('eye theta'); ylabel('eye phi'); axis equal
% title('left eye'); xlim([-23 23]); ylim([-23 23]); colorbar
% 
% for i =1:length(tl)
%     figure(1)
% subplot(1,2,1)
% plot(tr(i)-nanmean(tr),pr(i)-nanmean(pr),'.','Markersize',15,'Color', cmapVar(i,1,length(tl),jet)); axis square; hold on
% subplot(1,2,2)
% figure(2)
%  plot(tl(i)-nanmean(tl),pl(i)-nanmean(pl),'.','Markersize',15,'Color', cmapVar(i,1,length(tl),jet)); axis square; hold on
% 
% end
% %
% %
% %
% % % % % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% % % % %
%     end
% end
%%

clear all
% load('ACCAnalyzed_AllAnimals_010820_noDLS.mat')
% %  load('ACC_AllAnimals_021520_a.mat')
% load('DEINTERLACED_Analyzed_AllAnimals_052320_halfShift.mat')
   load('DEINTERLACED_Analyzed_AllAnimals_060220_halfShift.mat')


deInter=1;
if deInter
    frRate=60;
else
    frRate=30;
end
savePDF=1;
psfilename = 'C:\analysisPS_B.ps';

% Figure 2D: example of convergence during approach
% 
%  appRange=60:140;
% appRange=97-45:151+30;
% appTimes=97-45:106;
% % appRange=474-100:502+100
% 
% tR=Rtheta{useData(67)}(appRange,:); tR=tR-nanmean(tR);
% tL=Ltheta{useData(67)}(appRange,:);tL=tL-nanmean(tL);
% tilt=accelChannels{useData(67)}(appRange,2);  tilt=tilt-nanmean(tilt); tilt=medfilt1(tilt,8);
% roll=accelChannels{useData(67)}(appRange,1);  roll=roll-nanmean(roll); roll=medfilt1(roll,5);
% pR=Rphi{useData(67)}(appRange,:); pR=pR-nanmean(pR);
% pL=Lphi{useData(67)}(appRange,:);pL=pL-nanmean(pL);
% pVg=pR-pL; pVg=pVg-nanmean(pVg);
% hT=thetaHead{useData(67)}(1,appRange); hT=hT-nanmean(hT);
% gaze=(.5*(tR+tL))+hT'; %gaze=gaze-nanmean(gaze);
% %
% figure
% subplot(2,1,1)
%  appTimes=find(appEpoch{67}(appRange)==1)
% % appTimes=97:151;
% plot(tR); hold on;
% plot(tL); hold on;
% plot([appTimes(1) appTimes(1)],[-frRate frRate],'g'); plot([appTimes(end) appTimes(end)],[-frRate frRate],'r');
% ylim([-60 60]);
% ylabel('degrees');xlabel('time');
% legend('right eye','left eye');
% 
% subplot(2,1,2)
% plot(roll);hold on;
% plot(tilt);
% plot([appTimes(1) appTimes(1)],[-frRate frRate],'g'); plot([appTimes(end) appTimes(end)],[-frRate frRate],'r');
% hold on;
% % plot(tR(appRange)-tL(appRange));
% ylim([-60 60]); title('acc tilt');
% legend('roll','tilt');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%
% figure
% subplot(2,1,1)
% % appTimes=find(appEpoch{67}(appRange)==1)
% % appTimes=97:151;
% plot(tR-tL); hold on;
% plot(tilt); hold on;
% plot([appTimes(1) appTimes(1)],[-frRate frRate],'g'); plot([appTimes(end) appTimes(end)],[-frRate frRate],'r');
% ylim([-60 60]);
% ylabel('degrees');xlabel('time');
% legend('vergence','tilt');
% 
% subplot(2,1,2)
% plot(pVg);hold on;
% plot(roll);
% plot([appTimes(1) appTimes(1)],[-frRate frRate],'g'); plot([appTimes(end) appTimes(end)],[-frRate frRate],'r');
% hold on;
% % plot(tR(appRange)-tL(appRange));
% ylim([-60 60]); title('acc tilt');
% legend('phi diff','roll');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% 
% figure
% % appTimes=find(appEpoch{67}(appRange)==1)
% % appTimes=97:151;
% plot(gaze,'g'); hold on;
% plot(hT,'k'); hold on;
% % plot(tR-tL,'--b');
% plot([appTimes(1) appTimes(1)],[-frRate frRate],'g'); plot([appTimes(end) appTimes(end)],[-frRate frRate],'r');
% ylabel('degrees');xlabel('time');
% legend('gaze','head');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% 
% % Figure 2D: example of convergence during approach
% appRange=3676:3855;
% tR=Rtheta{useData(35)}; tR=tR-nanmean(tR);
% tL=Ltheta{useData(35)};tL=tL-nanmean(tL);
% tilt=accelChannels{useData(35)}(:,2);  tilt=tilt-nanmean(tilt); tilt=medfilt1(tilt,8);
% roll=accelChannels{useData(35)}(:,1);  roll=roll-nanmean(roll); roll=medfilt1(roll,5);
% 
% figure
% subplot(2,1,1)
% appTimes=find(appEpoch{35}(appRange)==1)
% plot(tR(appRange)); hold on; %3 minute segment
% plot(tL(appRange)); hold on;
% plot([appTimes(1) appTimes(1)],[-frRate frRate],'g'); plot([appTimes(end) appTimes(end)],[-frRate frRate],'r');
% ylim([-60 60]);
% ylabel('degrees');xlabel('time');
% legend('right eye','left eye');
% 
% subplot(2,1,2)
% plot(roll(appRange));hold on;
% plot(tilt(appRange));
% plot([appTimes(1) appTimes(1)],[-frRate frRate],'g'); plot([appTimes(end) appTimes(end)],[-frRate frRate],'r');
% hold on;
% % plot(tR(appRange)-tL(appRange));
% ylim([-60 60]); title('acc tilt');
% legend('roll','tilt');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%

clear all

% %  load('ACC_AllAnimals_021520_a.mat')
% load('DEINTER_Analyzed_AllAnimals_051820.mat');

deInter=1;
if deInter
    frRate=60;
    %     load('DEINTER_Analyzed_AllAnimals_051820.mat');
%load('DEINTERLACED_Analyzed_AllAnimals_052020.mat')
%    load('DEINTERLACED_Analyzed_AllAnimals_052320_halfShift.mat')
      load('DEINTERLACED_Analyzed_AllAnimals_060220_halfShift.mat')

  
else
    frRate=30;
    load('ACC_AllAnimals_021520_a.mat')
end

%%% subtract off median of gyro3, which is a small bias in measurement
gyroBias=nanmedian(gyro3All)
gyro3All = gyro3All-gyroBias;


savePDF=1;
psfilename = 'C:\analysisPS_B.ps';


% Figure 2D: scatter plot of R vs L eye yaw to show conv/divergence
clear numPts congMv convMv divMv Rmv Lmv
figure
for vid=1:length(useData)
    clear nframe
    nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
    nframe=min(nframe, length(appEpoch{vid}));
    tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
    tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
    useN = appEpoch{vid}(1:nframe)==0; use = appEpoch{vid}(1:nframe)==1;
    plot(tR(useN(1:20:end)),tL(useN(1:20:end)),'bo'); hold on; axis square; xlim([-50 50]);ylim([-50 50]);
    plot(tR(use(1:20:end)),tL(use(1:20:end)),'go');
    clear tR tL
    tR=dRtheta{useData(vid)}(1:nframe); tL=dLtheta{useData(vid)}(1:nframe);
    tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
    useN = appEpoch{vid}(1:nframe)==0; use = appEpoch{vid}(1:nframe)==1;
    xlabel('R eye theta'); ylabel('L eye theta');
    Rmv=nansum(tR(useN)>=0 & tL(useN)>=0);
    Lmv=nansum(tR(useN)<=0 & tL(useN)<=0);
    
    congMv(vid,1)=(Rmv+Lmv);
    convMv(vid,1)=(nansum(tR(useN)<=0 & tL(useN)>=0));
    divMv(vid,1)=(nansum(tR(useN)>=0 & tL(useN)<=0));
    
    clear Rmv Lmv
    Rmv=nansum(tR(use)>=0 & tL(use)>=0);
    Lmv=nansum(tR(use)<=0 & tL(use)<=0);
    congMv(vid,2)=(Rmv+Lmv);
    convMv(vid,2)=(nansum(tR(use)<=0 & tL(use)>=0));
    divMv(vid,2)=(nansum(tR(use)>=0 & tL(use)<=0));
    nPts(vid,1) = nansum((useN) & sum(~isnan(tR))' & sum(~isnan(tL))');
    nPts(vid,2)= nansum((use) & sum(~isnan(tR))' & sum(~isnan(tL))');
end
title('Figure 2D: R and L eye theta position');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%%% another example
% appRange=6843-60:6880+60
% tR=Rtheta{useData(40)}; tR=tR-nanmean(tR);
% tL=Ltheta{useData(40)};tL=tL-nanmean(tL);
% clear tilt
% tilt=accelChannels{useData(40)}(:,2); tilt=tilt-nanmean(tilt);
%
%
% figure
% subplot(2,1,1)
% appTimes=find(appEpoch{40}(appRange)==1)
% plot(tR(appRange)); hold on; %3 minute segment
% plot(tL(appRange)); hold on;
% plot([appTimes(1) appTimes(1)],[-frRate frRate],'g'); plot([appTimes(end) appTimes(end)],[-frRate frRate],'r');
% % ylim([-40 40]);
%
% subplot(2,1,2)
% plot(tilt(appRange));hold on;
% plot([appTimes(1) appTimes(1)],[-frRate frRate],'g'); plot([appTimes(end) appTimes(end)],[-frRate frRate],'r');
% hold on;
% plot(tR(appRange)-tL(appRange));
% ylim([-40 40]); title('acc tilt');

% Figure 2E: quantification of 2D, bar plots
congErr=nanstd(congMv)./sqrt(length(congMv))
convErr=nanstd(convMv)./sqrt(length(convMv))
divErr=nanstd(divMv)./sqrt(length(divMv))

eyeTypes = [nansum(congMv,1)./(nansum(nPts,1)); nansum(convMv,1)./(nansum(nPts,1)); nansum(divMv,1)./(nansum(nPts,1))];
eyeErr=[congErr./(nansum(nPts,1));convErr./(nansum(nPts,1));divErr./(nansum(nPts,1))];
figure;barweb(eyeTypes,eyeErr);
title('Figure 2E: proportions of congruent and incongruent eye movements')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
[h p]=ttest2(congMv(:,1),congMv(:,2))
[h p]=ttest2(convMv(:,1),convMv(:,2))
[h p]=ttest2(divMv(:,1),divMv(:,2))

% Figure 2F: hist of difference in eye yaw (aka vergence)
clear vergCounts vgDiff nBins vgHist

for c=0:1
    for vid=1:length(useData)
        nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
        nframe=min(nframe, length(appEpoch{vid}))
        tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
        tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
        vg=tR-tL;
        %         vgErr=nanstd(vg)./(sqrt(length(vg)));
        useN= appEpoch{vid}(1:nframe)==c;
        nBins=-90:5:90
        if sum(useN)>4
            if deInter
                vgHist = (hist(vg(useN),nBins))/sum(useN&~isnan(vg));
            else
                vgHist = (hist(vg(useN),nBins))/sum(useN&~isnan(vg)');
            end
        else
            vgHist=NaN;
        end
        vgDiff(vid,:,c+1)=vgHist;
    end
end

vergCounts(:,1)=nansum(vgDiff(:,:,1),1)./(length(useData))
vergCounts(:,2)=nansum(vgDiff(:,:,2),1)./(sum(~isnan(useData))-sum(isnan(vgDiff(:,1,2))))
vergErr=nanstd(vgDiff)./(sqrt(length(vgDiff))); vergErr=squeeze(vergErr);

figure; shadedErrorBar(nBins,vergCounts(:,1),vergErr(:,1),'b',1); hold on
shadedErrorBar(nBins,vergCounts(:,2),vergErr(:,2),'g',1); hold on

title('Figure 2F: vergence (R - L eye)'); axis square
xlabel('vergence (deg)'); ylabel('proportion of time');
ylim([0 .27])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
trMn=squeeze(nanmean(vgDiff,2))
[h, p]=ttest2(trMn(:,1),trMn(:,2))

% Figure 2G: correlation of eye thetas & eye phis
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
corrRAAll=[]
for vid=1:length(useData)
    nframe = min(length(dRtheta{useData(vid)}),length(dLtheta{useData(vid)}));
    dtR=dRtheta{useData(vid)}(1:nframe); dtL=dLtheta{useData(vid)}(1:nframe);
    nonapp=appEpoch{vid}==0;
    use = (nonapp==1)';%& ~isnan(dtR(1:nframe));
    if sum(use)>3
        [corrR lagsR]= nanxcorr(dtR(use),dtL(use),frRate,'coeff');
        uselagsR=(lagsR>=-frRate& lagsR<=frRate);
        
        clear nframe
        nframe = min(length(dRphi{useData(vid)}),length(dLphi{useData(vid)}));
        dpR=dRphi{useData(vid)}(1:nframe); dpL=dLphi{useData(vid)}(1:nframe);
        [corrL lagsL]= nanxcorr(dpR(use),dpL(use),frRate,'coeff');
        uselagsL=(lagsL>=-frRate & lagsL<=frRate);
        
    else
    end
    use=appEpoch{vid}==1;
    if sum(use)>3 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dtR(use),dtL(use),frRate,'coeff');
        uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
        
        nframe = min(length(dRphi{useData(vid)}),length(dLphi{useData(vid)}));
        dpR=dRphi{useData(vid)}(1:nframe); dpL=dLphi{useData(vid)}(1:nframe);
        [corrLA lagsLA]= nanxcorr(dpR(use),dpL(use),frRate,'coeff');
        uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
    else
    end
    
    if sum(uselagsR)==2*frRate+1 & sum(uselagsL)==2*frRate+1 &sum(uselagsRA)==2*frRate+1& sum(uselagsLA)==2*frRate+1
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        %         errVidR(vid)=errRvid;errVidL(vid)=errLvid;
    else
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.7 .7]);
 xlim([41 81]); 
ylim([-.2 .4])
axis square
L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');

legend(L,{'dTheta non-app','dPhi non-app','dTh Approach','dPhi app'}); title('Figure 2G: mean between eye corr, dtheta & dphi');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%% FIGURE 3: during approaches, eyes are more centered (smaller vergence), due to less head movements in pitch and roll

% Figure 3A: overlaid traces of pitch & eye vergence

% Figure 3B: scatter & corr of pitch (tilt) and vergence
close all
figure; clear pts
for vid=1:length(useData)
    clear tiltVid
    nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
    tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
    tR=nanmedian(tR)-tR; tL=nanmedian(tL)-tL;
    
    tiltVid=accelChannels{useData(vid)}(:,2);
    tiltVid = medfilt1(tiltVid,8);
    
    tiltVid=tiltVid-nanmedian(tiltVid);
    vg=tR-tL;
    clear use
    for c=0:1
        use = find(appEpoch{vid}==c);
        if c==0
            plot(tiltVid(use(1:300:length(use))),vg(use(1:300:length(use))),'b.'); hold on;
            xlim([-60 60]);ylim([-60 60]);
            pts(vid,1,1)=length(use(1:300:end))
            
            axis square
        else
            plot(tiltVid(use(1:300:length(use))),vg(use(1:300:length(use))),'g.');
            pts(vid,2,1)=length(use(1:300:end))
            
        end
        %     pts(vid,1,1)=length(use(1:70:end))
        %     pts(vid,2,1)=length(use(1:50:end))
        pts(vid,c+1,2)=length(tiltVid(use));
    end
end
v=squeeze(nansum(pts,1))

close all
for c=0:1
    clear use
    use = find(appAll==c);
    figure(1)
    plot(tiltAll(use(1:60:end)),vergDlc(use(1:60:end)),'.'); axis equal; hold on;
    pts(c+1,1)=length(use(1:60:end));
    pts(c+1,2)=length(use(1:60:end))/length(tiltAll);
    xlabel('acc tilt'); ylabel('eye vergence');
    ylim([-90 90]); xlim([-90 90]);
    R=corrcoef(tiltAll(use), vergDlc(use),'Rows','pairwise')
    text(-80,(80-c*5), ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
    title('Figure 3B: acc tilt and dlc eye vergence');
    if c==1
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(2)
    [corr lags]=(nanxcorr(tiltAll(use),vergDlc(use),frRate,'coeff'));
    plot(lags,corr); hold on; ylabel('correlation coeff');
    % axis square; ylim([-.2 .8]); xlim([-15 15])
    title('Figure 3C: corr of acc tilt & eye vergence');
    if c==1
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
    
    figure(3)
    clear nBins h
    nBins=-90:5:90
    h=hist(tiltAll(use),nBins);
    plot(nBins, h/length(tiltAll(use))); hold on
    title('Figure 3D: acc pitch (tilt) is less during approach)')
    pitchHist(:,c+1)=h;
    if c==1
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
end

clear tiltHist
for c=0:1
    for vid=1:length(useData)
        clear useN tilt tiltHist
        tilt=accelChannels{useData(vid)}(:,1);
        useN= appEpoch{vid}==c;
        nframe=min(length(tilt), length(useN));
        useN=useN(1:nframe);tilt=tilt(1:nframe);
        tilt=tilt-nanmean(tilt);
        nBins=-90:5:90
        if sum(useN)>1
            tiltHist = (hist(tilt(useN),nBins))./sum((~isnan(tilt(useN))))% & find(useN)');
        else
            tiltHist=NaN;
        end
        tiltFull(vid,:,c+1)=tiltHist;
    end
end

tiltCounts(:,1)=nansum(tiltFull(:,:,1),1)./(length(useData));
tiltCounts(:,2)=nansum(tiltFull(:,:,2),1)./((length(useData))-sum(isnan(tiltFull(:,1,2))))
tiltErr=nanstd(tiltFull)./(sqrt(length(tiltFull))); tiltErr=squeeze(tiltErr);

figure; shadedErrorBar(nBins,tiltCounts(:,1),tiltErr(:,1),'b',1); hold on
shadedErrorBar(nBins,tiltCounts(:,2),tiltErr(:,2),'g',1); hold on


title('Figure 3B: pitch'); axis square
xlabel('pitch (deg)'); ylabel('proportion of time');
ylim([0 .27])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
tiltMns=squeeze(nanmean(tiltFull,2));
[h p]=ttest2(tiltMns(:,1),tiltMns(:,2))


% Figure 3C: example trace, less change in pitch during approach

% Figure 3D: hist or corr, less change in pitch during approach
% corr & hist above in loop

% Figure 3E: overlaid traces of roll & eye phi

% Figure 3F: scatter & corr of roll & eye phi
close all
for c=0:1
    clear use
    use = find(appAll==c);
    figure(4)
    plot(rollAll(use),dlcPhi(use),'.'); axis equal; hold on; lsline
    xlabel('acc roll'); ylabel('eye phi');
    ylim([-90 90]); xlim([-90 90]);
    R=corrcoef(rollAll(use), dlcPhi(use),'Rows','pairwise');
    text(-80,(80-c*5), ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
    title('Figure 3F: acc roll & diff of eye phi)');
    if c==1
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
        
    end
    
    figure(5)
    [corr lags]=(nanxcorr(rollAll(use),dlcPhi(use),frRate,'coeff'));
    plot(lags,corr);axis square; hold on
    title('Figure 3G: corr of acc roll & diff between eye phi'); ylabel('corr coeff')
    ylim([-1 1]);
    if c==1
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
    figure(6)
    clear nBins h
    nBins=-90:5:90;
    h=hist(rollAll(use),nBins);
    plot(nBins, h/length(tiltAll(use))); hold on; title('acc roll')
    title('Figure 3H: roll during app is less')
    if c==1
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
        
    end
end

% Figure 3G: example trace, less change in roll during approach

% Figure 3H: hist or corr, less change in roll during approach
% hist and corr above in loop

%% FIGURE 4: Coordination of eyes and head

% Figure 4A: hist of head yaw app and non approach (not different)



figure
for c=0:1
    clear nBins h use
    use = find(appAll==c);
    nBins=-30:2:30
    h=hist(gyro3All(use),nBins);
    headHist(:,c+1)=h;
    plot(nBins, h/length(gyro3All(use))); hold on; title('head yaw from gyro')
    title('Figure 4A: head angle is not diff between non-app and app')
    if c==1
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
end

for c=0:1
    for vid=1:length(useData)
        use=appEpoch{vid}==c;
        headMn(vid,c+1)=nanmean(accelChannels{useData(vid)}(use,6));
    end
end

[h p]= ttest2(headMn(:,1),headMn(:,2))

% Figure 4B: hist of mean eye yaw, app and non-app
d_mnEyeAll=[];mnPhiAll=[];allDLChead=[];mnEyeAll=[]; mouseSpAll=[];
for vid=1:length(useData)
    nframe=min(length(dRtheta{useData(vid)}),length(dLtheta{useData(vid)}));
    nframe=min(nframe, length(d_Theta{useData(vid)}));
    mnEye =.5*(Rtheta{useData(vid)}(1:nframe)+Ltheta{useData(vid)}(1:nframe));
    mnEye=mnEye-nanmean(mnEye);
    
    mnEyeD =.5*(dRtheta{useData(vid)}(1:nframe)+dLtheta{useData(vid)}(1:nframe));
    mnEyeD=mnEyeD-nanmean(mnEyeD);
    mnPhi =.5*(dRphi{useData(vid)}(1:nframe)+dLphi{useData(vid)}(1:nframe));
    mnPhi=mnPhi-nanmean(mnPhi);
    
    dHead=d_Theta{useData(vid)}(1:nframe); dHead=dHead-nanmean(dHead);
    allDLChead = [allDLChead dHead'];
    
    speed =mouseSp{useData(vid)}(1:nframe);
    %speed =speed-nanmean(speed);  commented out 3/11 since this makes speeds inaccurate
    mouseSpAll=[mouseSpAll speed];
    
    if deInter
        mnEyeAll=[mnEyeAll mnEye];
        d_mnEyeAll=[d_mnEyeAll mnEyeD];
        mnPhiAll=[mnPhiAll mnPhi];
    else
        mnEyeAll=[mnEyeAll mnEye'];
        d_mnEyeAll=[d_mnEyeAll mnEyeD'];
        mnPhiAll=[mnPhiAll mnPhi'];
    end
    
end

figure
for c=0:1
    hold on
    clear nBins h use hp
    use = find(appAll==c);
    nBins= -60:2:60
    h=hist(mnEyeAll(use),nBins);
    binoc(:,c+1)=h;
    
    
    subplot(1,2,1)
    plot(nBins,h/length(mnEyeAll(use))); hold on; title('mn Eye Yaw'); axis square
    plot([-20,-20],[.15,0],'k--');
    plot([20,20],[.15,0],'k--');
    
    nBins=-30:2:30
    hp=hist(d_mnEyeAll(use),nBins);
    subplot(1,2,2)
    plot(nBins, hp/length(d_mnEyeAll(use))); axis square; hold on; title('mn d Eye Theta')
    deltaBinoc(:,c+1)=hp;
    if c==1
        suptitle('Figure 4B: mn Eye yaw falls within binocular zone');
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
end


%  [h pB]=kstest2(binoc(:,1),binoc(:,2))
[h pdB]=ttest2(deltaBinoc(:,1),deltaBinoc(:,2))

for vid=1:length(useData)
    tR=Rtheta{useData(vid)}; tL=Ltheta{useData(vid)};
    nframe = min(length(tR), length(appEpoch{vid}))
    for c=0:1
        use=appEpoch{vid}(1:nframe)==c;
        tR=Rtheta{useData(vid)}; tL=Ltheta{useData(vid)};
        eyePos=.5*(tR+tL); eyePos=eyePos-nanmedian(eyePos);
        use=appEpoch{vid}==c;
        expMn(vid,c+1)=nanmean(eyePos(use));
    end
end

[h p]=ttest2(expMn(:,1),expMn(:,2))
% Figure 4C: ex trace when dTh is 0, eyes are also 0, then hist

figure
for c=0:1
    hold on
    clear nBins h use hp
    use = (appAll==c);
    %still=(allDLChead<.5 & allDLChead>-.5);
    still=(gyro3All<.25 & gyro3All>-.25);
    nBins= -10:1:10;
    h=hist(d_mnEyeAll(use&still),nBins);
    stillHist(:,c+1)=h;
    subplot(1,2,1)
    plot(nBins,h/sum(use&still)); hold on; title('mn Eye Yaw'); axis square
    plot([-3,-3],[.6,0],'k--');
    plot([3,3],[.6,0],'k--');
    xlim([-10 10]);
    
    prop3(c+1,1) =1- sum(h(nBins<-1 | nBins>1))/sum(h);
    hp=hist(mnPhiAll(use& still),nBins);
    subplot(1,2,2)
    plot(nBins, hp/sum(use&still)); axis square; hold on; title('mn Eye Phi')
    prop3(c+1,2)=1- sum(hp(nBins<-5 | nBins>5))/sum(hp);
    
    plot([-3,-3],[.6,0],'k--');
    plot([3,3],[.6,0],'k--');
    if c==1
        suptitle('Figure 4C: when head angle is not changing, eyes are also still');
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
end

[h still]=ttest2(stillHist(:,1),stillHist(:,2))
%
figure
clear nBins h use hp
allRunFilt=medfilt1(mouseSpAll,15);
stationary = allRunFilt<1;
% stationary = mouseSpAll<1;

nBins= -40:1:40;
hs=hist(d_mnEyeAll(stationary),nBins);
plot(nBins,hs/sum(stationary)); hold on; title('mn Eye Yaw, stationary'); axis square

hmv=hist(d_mnEyeAll(stationary==0),nBins);
plot(nBins,hmv/sum(stationary==0)); hold on; title('mn Eye Yaw, moving'); axis square
% h=hist(d_mnEyeAll(appAll==1),nBins);
% plot(nBins,h/sum(appAll==1)); hold on; title('mn Eye Yaw, approach'); axis square
legend('stationary','moving');
ylim([0 .9]); xlim([-10 10]);
title('Figure 4C: when mouse is stationary, so are eyes');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


for vid=1:length(useData)
    speed=mouseSp{useData(vid)}(:,1:end-1);
    runFilt = medfilt1(speed,(frRate/2));
    
    running =runFilt>=1;
    if deInter
        mnEyes(vid,1) = nanmedian((.5*(dRtheta{useData(vid)}(:,running==0)+ dLtheta{useData(vid)}(:,running==0))));
        mnEyes(vid,2) = nanmedian((.5*(dRtheta{useData(vid)}(:,running==1)+ dLtheta{useData(vid)}(:,running==1))));
    else
        mnEyes(vid,1) = nanmedian((.5*(dRtheta{useData(vid)}(running==0,:)+ dLtheta{useData(vid)}(running==0,:))));
        mnEyes(vid,2) = nanmedian((.5*(dRtheta{useData(vid)}(running==1,:)+ dLtheta{useData(vid)}(running==1,:))));
        
    end
end

[h stat]=kstest2(mnEyes(:,1),mnEyes(:,2))


% Figure 4D: scatter of change in head yaw and change in eye yaw, shows
% mostly congruent but not all
figure; clear samp
for c=1:2
    clear use
    use = find(appAll==c-1);
    samp(c,:)=randsample(use,200);
    plot(gyro3All(samp(c,:)), d_mnEyeAll(samp(c,:)),'.')
    axis([-35 35 -35 35]);
    axis square
    hold on
end
xlabel('gyro yaw'); ylabel('d eye yaw');

title('Figure 4D: head & eye thetas, not all compensatory')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


fullDist=[]
for vid=1:length(useData)
    clear distance
    distance=dist{useData(vid)}
    fullDist=[fullDist distance];
end

%%% distance to cricket histograms
clear crickHist crickCounts crickErr crickMns crickFull 
for c=0:1  %%% approach vs non-approach
    for vid=1:length(useData)
        clear useN crick crickHist
        crick=dist{useData(vid)}
        useN= appEpoch{vid}==c;
        nframe=min(length(crick), length(useN));
        useN=useN(1:nframe);crick=crick(1:nframe);
        nBins=0:4:50
        if sum(useN)>1
            crickHist = (hist(crick(useN),nBins))./sum((~isnan(crick(useN))))% & find(useN)');
        else
            crickHist=NaN;
        end
        crickFull(vid,:,c+1)=crickHist;
    end
end

crickCounts(:,1)=nansum(crickFull(:,:,1),1)./(length(useData))
crickCounts(:,2)=nansum(crickFull(:,:,2),1)./((length(useData))-sum(isnan(crickFull(:,1,2))))
crickErr=nanstd(crickFull)./(sqrt(length(crickFull))); crickErr=squeeze(crickErr);
figure; shadedErrorBar(nBins,crickCounts(:,1),crickErr(:,1),'b',1); hold on
shadedErrorBar(nBins,crickCounts(:,2),crickErr(:,2),'g',1); hold on
title('distance to cricket');
ylabel ('proportion of time'); xlabel('cm');

crickMns=squeeze(nanmean(crickFull,2));
[h p]=ttest2(crickMns(:,1),crickMns(:,2))
%%
clear ptsUsed
figure
for c=0:1
    clear use
    use = find(appAll==c);
    ptsUsed(c+1,1)=length(use(1:100:end));
    ptsUsed(c+1,2)=length(use(1:100:end))/length(gyro3All(use));
    %     samp(c,:)=randsample(use,7375);
    %     plot(gyro3All(samp(c,:)), d_mnEyeAll(samp(c,:)),'.')
    plot(gyro3All(use(1:100:end)), d_mnEyeAll(use(1:100:end)),'.'); hold on
    %     axis([-35 35 -35 35]);
    axis([-15 15 -15 15]);
    axis square
    hold on
end
xlabel('gyro yaw'); ylabel('d eye yaw');
% Figure 4E: corr of change in head and eye yaw, congruence
% (non-compensation) at short time scales

%%
clear corrYaw corrAllYaw err errA lagsYaw
for c=0:1
    for vid=1:length(useData)
        nframe = min(length(accelChannels{useData(vid)}(:,6)),length(dRtheta{useData(vid)}));
        nframe = min(nframe, length(dLtheta{useData(vid)}));
        nframe=min(nframe, length(appEpoch{vid}));
        use=appEpoch{vid}(1:nframe)==c;
        dth=d_Theta{useData(vid)}(1:nframe);
        g3=accelChannels{useData(vid)}(1:nframe,6); dtR=dRtheta{useData(vid)}(1:nframe); dtL=dLtheta{useData(vid)}(1:nframe);
        mnEye=.5*(dtR+dtL); g3=g3-nanmean(g3); mnEye=mnEye-nanmean(mnEye);
        if sum(use)>5 & sum(~isnan(g3(use)))>10
            [corrYaw lagsYaw]= nanxcorr(g3(use),mnEye(use),frRate,'coeff');
            uselags=(lagsYaw>=-frRate& lags<=frRate);
        end
        if sum(uselags)==2*frRate+1
            corrAllYaw(vid,:,c+1)=corrYaw(uselags);
        else
        end
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
err= nanstd(corrAllYaw(:,:,1))/(sqrt(length(corrAllYaw)));
errA= nanstd(corrAllYaw(:,:,2))/(sqrt(length(corrAllYaw)));
shadedErrorBar(1:size(corrAllYaw,2),nanmean(corrAllYaw(:,:,1),1),err,'-b',1); hold on
shadedErrorBar(1:size(corrAllYaw,2),nanmean(corrAllYaw(:,:,2),1),errA,'-g',1);

plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]);
ylim([-.55 .3]);
% xlim([21 41]); %15 and 46 ==500 ms
% xlim([8.5 53.5])
xlim([frRate-((frRate/2)+1) frRate+((frRate/2)+1)])


xlabel('time'); ylabel('correlation coeff');
axis square
title('Figure 4E: eye movements are mostly compensatory');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%% FIGURE 5: Non-compensatory eye movements

%%% clustering of eye movement!!!

%%% cluster on dHead vs dEye
%%%old way
clear mvmts


mvmts=[gyro3All(appAll==1); d_mnEyeAll(appAll==1); mnEyeAll(appAll==1)];
%mvmts = [gyro3All(appAll==1) + d_mnEyeAll(appAll==1)];
gm = fitgmdist(mvmts',3,'Replicates',10);
idx = cluster(gm,mvmts');

figure
gscatter(gyro3All(appAll==1),d_mnEyeAll(appAll==1),idx~=2); axis equal
title('cluster on dHead vs dEye vs Eye 3 clust')
xlim([-25 25]); ylim([-25 25]); %hold on; plot([-25 25], [25 -25],'r')

%%% cluster on dHead vs dEye vs Eye
clear mvmts
mvmts=[gyro3All(appAll==1); d_mnEyeAll(appAll==1); mnEyeAll(appAll==1)];
%mvmts = [gyro3All(appAll==1) + d_mnEyeAll(appAll==1)];
gm = fitgmdist(mvmts',3,'Replicates',10);
idx = cluster(gm,mvmts');

X=mvmts;
figure
gscatter(gyro3All(appAll==1),d_mnEyeAll(appAll==1),idx~=2); axis equal
title('cluster on dHead vs dEye vs Eye 3 clust')
xlim([-25 25]); ylim([-25 25]); %hold on; plot([-25 25], [25 -25],'r')





%%% cluster on dHead vs dEye

%
dGz= gyro3All(appAll==1) + d_mnEyeAll(appAll==1);

clear mvmts
mvmts=[gyro3All(appAll==1); d_mnEyeAll(appAll==1)] %; mnEyeAll(appAll==1)];
%mvmts = [gyro3All(appAll==1) + d_mnEyeAll(appAll==1)];
gm = fitgmdist(mvmts',3,'Replicates',10);
idx = cluster(gm,mvmts');

X=mvmts;
figure
use=find(appAll==1) % subset op pts - used in paper
gscatter(gyro3All(use(1:5:end)),d_mnEyeAll(use(1:5:end)),idx(1:5:end)); axis equal
title('cluster on dHead vs dEye 3 clust')
xlim([-25 25]); ylim([-25 25]); %hold on; plot([-25 25], [25 -25],'r')




clear mvmts
mvmts=[gyro3All(appAll==1); d_mnEyeAll(appAll==1)] %; mnEyeAll(appAll==1)];
%mvmts = [gyro3All(appAll==1) + d_mnEyeAll(appAll==1)];
gm = fitgmdist(mvmts',3,'Replicates',10);
idx = cluster(gm,mvmts');

bins = -16:0.5:16;
clustGz = gyro3All(appAll==1)+ d_mnEyeAll(appAll==1);
clustGz = clustGz-nanmedian(clustGz);
nclustPts = sum(~isnan(clustGz));

bins = -16:1:16;
clear gazeHist
gazeHist = hist(clustGz,bins)/nclustPts;
figure
plot(frRate*bins,gazeHist,'b');
xlim([-15 15]*frRate)
xlabel('gaze velocity deg/sec')


clear gazeHist
gazeHist(:,1) = hist(clustGz(idx~=2),bins)/nclustPts;
gazeHist(:,2) = hist(clustGz(idx==2),bins)/nclustPts;
%gazeHist(gazeHist==0) = NaN;
gazeHist(abs(bins)>9) = NaN;
gazeHist(abs(bins)<3,2)=NaN;
%gazeHist = log10(gazeHist); gazeHist(gazeHist<-3)=-3;
figure
plot(frRate*bins,gazeHist(:,1),'k');
hold on
plot(frRate*bins,gazeHist(:,2),'r','LineWidth',2);
xlim([-16 16]*frRate); %ylim([0 0.1])
xlabel('gaze velocity deg/sec')

bins = 0.25:0.5:15;
gazeHist = hist(abs(clustGz),bins)/nclustPts;
figure
plot(bins,gazeHist)



X=mvmts;
figure
gscatter(gyro3All(appAll==1),d_mnEyeAll(appAll==1),idx~=2); axis equal
title('cluster on dHead vs dEye 3 clust')
xlim([-25 25]); ylim([-25 25]); %hold on; plot([-25 25], [25 -25],'r')

X=mvmts;
figure
gscatter(gyro3All(appAll==1),dGz,idx~=1); axis equal
title('cluster on dHead vs dEye 3 clust')
xlim([-25 25]); ylim([-25 25]); %hold on; plot([-25 25], [25 -25],'r')
%
% %%% cluster on dGz
% clear mvmts
% mvmts = dGz;
% gm = fitgmdist(mvmts',2,'Replicates',10);
% idx = cluster(gm,mvmts');
% X=mvmts;
% figure
% gscatter(gyro3All(appAll==1),d_mnEyeAll(appAll==1),idx); axis equal
% title('cluster on dGz')
% xlim([-25 25]); ylim([-25 25]); %hold on; plot([-25 25], [25 -25],'r')
% %hold on; plot([-25 25], [21 -29],'r'); plot([-25 25], [29 -21],'r')


%%% histogram of gaze
nBins=-20:0.25:20
dGzHist=hist(dGz,nBins);
figure;
plot(nBins, dGzHist);

% %%%% try to directly fit from dGz. looks like it ends up being a bit restrictive
clear err
sigmas = 0.1:0.1:10; A = max(dGzHist);
for s= 1:length(sigmas);
    g = A*exp(-(nBins.^2)/(0.5*sigmas(s)^2));
    err(s) = sum((g-dGzHist).^2);
end
figure
plot(sigmas,err)
[e min_s] = min(err);
sig = sigmas(min_s)*1.5; fit = A*exp(-(nBins.^2)/(0.5*sig^2));
figure;
plot(nBins, dGzHist); hold on; plot(nBins,fit);
pfit = fit./(dGzHist+0.0001); pfit(pfit>1)=1;  %%% add 0.0001 to prevent 0/0
figure
plot(nBins,pfit)
prob = interp1(nBins,pfit,dGz);
clust1 = rand(size(prob))<prob;

%%% super simple! threshold on gaze velocity
dGzV = dGz*frRate;
gzThresh = 180;  %% was 180
clust1 = abs(dGzV)>gzThresh;
figure
useTime=find(appAll==1) % subset op pts - used in paper
gscatter(gyro3All(useTime(1:10:end)),d_mnEyeAll(useTime(1:10:end)),clust1(1:10:end)); axis equal
% useGyro=gyro3All(appAll==1); useEyes=d_mnEyeAll(appAll==1);
title('clutering based on gaussian fit to dGz')
xlim([-15 15]); ylim([-15 15]);
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end



gzVall = (gyro3All + d_mnEyeAll)*frRate;
figure
[gzHist bins] = hist(gzVall,-600:5:600);
semilogy(bins,gzHist/nansum(gzHist));
xlabel('gaze velocity deg/sec'); hold on
sPts = bins>=gzThresh;
semilogy(bins(sPts),gzHist(sPts)/nansum(gzHist),'r');
sPts = bins<=-gzThresh;
semilogy(bins(sPts),gzHist(sPts)/nansum(gzHist),'r');
plot([gzThresh gzThresh],[0.0001 0.1],'r:')
plot(-[gzThresh gzThresh],[0.0001 0.1],'r:')


if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


figure
gscatter(gyro3All(appAll==1)-gyroBias,d_mnEyeAll(appAll==1),idx); axis equal; hold on
title('all approach pts only')
xlim([-25 25]); ylim([-25 25]); hold on; plot([-25 25],[25 -25],'g')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


% use ~clust1
% fullData=[gyro3All;d_mnEyeAll]%;mnEyeAll];
dGzFull= gyro3All + d_mnEyeAll;
fullData=[gyro3All; d_mnEyeAll];
% fullData=[gyro3All; dGzFull];
idxAll=cluster(gm, fullData');
%
% figure
% gscatter(gyro3All-nanmedian(gyro3All),dGzFull,idxAll); axis equal; hold on

[sacc, clust]= max([nanmean(abs(dGzFull(idxAll==1))),nanmean(abs(dGzFull(idxAll==2))), nanmean(abs(dGzFull(idxAll==3)))])

% figure;
% [corr lags]=nanxcorr(gyro3All(idxAll==clust),d_mnEyeAll(idxAll==clust),frRate,'coeff')
% plot(lags,corr)
%
% propNC(:,1)=(length(find(~appAll'&idxAll==clust)))/(sum(appAll==0));
% propNC(:,2)=(length(find(appAll'&idxAll==clust)))/(sum(appAll==1));
% propC(:,1)=(length(find(~appAll'&idxAll~=clust)))/(sum(appAll==0));
% propC(:,2)=(length(find(appAll'&idxAll~=clust)))/(sum(appAll==1));
%
% figure; bar([propC propNC]);

% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

% mnEye=.5*(Ltheta{useData(67)}+Rtheta{useData(67)}); mnEye=mnEye-nanmean(mnEye);
% g3=accelChannels{useData(67)}(:,6); hT=thetaHead{useData(67)};
% newData=[g3(1:end-1) diff(mnEye) mnEye(1:end-1);];
% idx2=cluster(gm,newData)
% figure
% gscatter(g3(1:end-1),diff(mnEye),idx2); axis equal
%
% figure
% plot(mnEye+hT'); hold on;
% plot(find(idx2==1),mnEye(idx2==1),'og')


%%
maxlag = frRate;
thbins = -60:5:60;
skip = 1; %%% only shows figures at this interval
nthresh  = 60;
headAll=[];
eyeAll=[]; g3All= []; azAll = []; spAll  = []; appAll = []; distAll = [];


ns = 0; saccHeadAll = []; saccEyeAll = []; saccAppAll = []; saccVidAll= []; saccAzAll = []; saccThAll = []; saccEyeRawAll=[]; timetoAppAll =[];
% figure

%%% non-app example is vid 114, approach example is vid 48
allS = 0; clear dheadStable eyeStable dgzStable %%% store out stable periods
for i = 1:length(appEpoch)
    vid = useData(i);
    %%% get approaches
    app = appEpoch{i};
    nonapp=~appEpoch{i};
    %%% get left eye positions
    lth = Ltheta{vid} - nanmedian(Ltheta{vid});
    dlth = dLtheta{vid};
    nl(i) = sum(~isnan(lth(app))); %%% # good eye approach points
    lthHist(:,1,i) = hist(lth(app),thbins)/nl(i);
    lthHist(:,2,i) = hist(lth(~app),thbins)/sum(~isnan(lth(~app)));
    if nl(i)<nthresh
        lthHist(:,:,i) = NaN;
    end
    %%% get right eye positions
    rth = Rtheta{vid} - nanmedian(Rtheta{vid});
    drth = dRtheta{vid};
    nr(i) =sum(~isnan(rth(app))); %%% # good eye approach points
    rthHist(:,1,i) = hist(rth(app),thbins)/nr(i);
    rthHist(:,2,i) = hist(rth(~app),thbins)/sum(~isnan(rth(~app)));
    if nr(i)<nthresh
        rthHist(:,:,i) = NaN;
    end
    
    
    %%% get eye phi
    rphi = Rphi{vid}-nanmean(Rphi{vid});
    drphi = dRphi{vid};
    nrp(i) =sum(~isnan(rphi(app))); %%% # good eye approach points
    rphiHist(:,1,i) = hist(rphi(app),thbins)/nrp(i);
    rphiHist(:,2,i) = hist(rphi(~app),thbins)/sum(~isnan(rphi(~app)));
    if nr(i)<nthresh
        rphiHist(:,:,i) = NaN;
    end
    
    % get left eye phi
    lphi = Lphi{vid}-nanmean(Lphi{vid});
    dlphi = dLphi{vid};
    nlp(i) =sum(~isnan(lphi(app))); %%% # good eye approach points
    lphiHist(:,1,i) = hist(lphi(app),thbins)/nrp(i);
    lphiHist(:,2,i) = hist(lphi(~app),thbins)/sum(~isnan(lphi(~app)));
    if nr(i)<nthresh
        lphiHist(:,:,i) = NaN;
    end
    
    
    %%% get head positions
    hth = thetaHead{vid};
    dth = d_Theta{vid};
    azdeg = az{vid}*180/pi;
    crickD = dist{vid};
    
    %%% get accelerometers
    
    %     figure
    %     hist(azdeg(app==1))
    %   if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    
    if exist('accelData','var')
        tilt = accelChannels{vid}(:,2);
        roll = accelChannels{vid}(:,1);
        acc_dth = accelChannels{vid}(:,6);
        acc_dth=acc_dth-gyroBias;
        acc_hth = nancumsum([hth(1) acc_dth(1:end-1)'],[],2);
    else
        display('no acc')
    end
    
    
    %%% get head positions
    %     hth = thetaHead{vid};
    %     dth = d_Theta{vid};
    %     azdeg = az{vid}*180/pi;
    
    hth = acc_hth;
    dth = acc_dth;
    azdeg = az{vid}*180/pi;
    
    %%% azimuth vs eye histograms
    az_hist(:,1,i) = hist(-azdeg(app),thbins)/sum(~isnan(azdeg(app)));
    n= sum(~isnan(azdeg(app)) & ~isnan(lth(app)'));
    azthL_hist(:,1,i) = hist(-azdeg(app)-lth(app)',thbins)/n;
    n= sum(~isnan(azdeg(app)) & ~isnan(rth(app)'));
    azthR_hist(:,1,i) = hist(-azdeg(app)-rth(app)',thbins)/n;
    
    %%% alignment of eyes during approaches
    %%% vergence is cool! it gets very tight around 0 during approaches
    vergence = rth-lth;
    n= sum(~isnan(vergence(app)));
    vergeHist(:,1,i) = hist(vergence(app),thbins)/n;
    vergeHist(:,2,i) = hist(vergence(~app),thbins)/sum(~isnan(vergence(~app)));
    
    %%% mean eye theta is most important for stabilization
    mnEyeTh = 0.5*(rth+lth);
    n= sum(~isnan(azdeg(app)) & ~isnan(mnEyeTh(app)'));
    azthRL_hist(:,1,i) = hist(-azdeg(app)-mnEyeTh(app)',thbins)/n;
    
    %%% gaze is the sum of head position + mean eye position
    %%% key variable!!!
    if deInter
        gaze = hth + mnEyeTh;
    else
        gaze = hth + mnEyeTh';
    end
    
    %%% find the longest approach, to use as example image
    appPts = find(app);
    newApp = [1 find(diff(appPts)>1)+1];
    endApp = [find(diff(appPts)>1)-1 length(appPts)];
    
    dur = endApp - newApp;
    [mx longest] = max(dur);
    mainApp = appPts(newApp(longest)  : endApp(longest)); %%% approach time only
    %%% add 2 secs on either side
    try
        appStart = max(appPts(newApp(longest))-(2*frRate),1); appOffset = appPts(newApp(longest))-appStart;
        appEnd = min(appPts(endApp(longest))+(2*frRate),length(app)); endOffset = appPts(endApp(longest))-appStart;
        appRange = appStart  : appEnd;
    catch
        appRange = appPts(newApp(longest)  : endApp(longest));
        appOffset =0;
    end
    
    
    %%% get rid of large jumps
    
    % hthnonan = hth;
    %hthnonan(abs(diff(hth))>90)=NaN;
    hthnonan = acc_hth;
    
    % hthApp = hth(appRange)-nanmedian(hth(mainApp));
    hthApp = acc_hth(appRange);
    hthApp = hthApp - nanmedian(acc_hth(mainApp));
    %hthApp = mod(hthApp + 180,360)-180;
    
    
    gzApp = hthApp' +0.5*( rth(appRange) +lth(appRange))';
    mnEyeApp =  0.5*(rth(appRange) +lth(appRange))';
    
    dhthnonan = acc_dth;
    dhthnonan(abs(diff(dth))>90)=NaN;
    
    dhthApp = acc_dth(appRange)-nanmedian(acc_dth(mainApp));
    dhthApp = mod(dhthApp + 180,360)-180;
    
    %%% calculate change in position at different lags, as measure of stability
    
    %%% draw figures
    
    if round(i/skip)==i/skip
        
        %          figure
        
        %         subplot(6,1,1);
        % %         plot(hthnonan,'k'); hold on;
        %         plot(rth,'r');hold on; plot(lth,'b'); legend('right th','left th');
        %         plot(find(app),ones(sum(app),1)*90,'g.'); %ylim([-180 180])
        roll = accelChannels{vid}(:,1); roll=roll-nanmean(roll); roll=medfilt1(roll,5);
        tilt=accelChannels{vid}(:,2); tilt=tilt-nanmean(tilt); tilt=medfilt1(tilt,8);
        g3=accelChannels{vid}(:,6); g3=g3-gyroBias;
        
        %dGaze=diff(gaze);
        dGaze = g3(1:end-1) + 0.5*(drth+dlth)';
        
        
        %         figure
        %         plot(g3(1:end-1), 0.5*(drth+dlth),'.'); axis([-20 20 -20 20])
        
        % nBins=-20:0.25:20;
        
        %%% assign clusters based on previous GMM; clust for analysis is set
        %%% above
        if deInter
        newData=[g3(1:end-1)' ; 0.5*(drth+dlth)];
        else
          newData=[g3(1:end-1) ; 0.5*(drth+dlth)];
        end
    %    idx2=cluster(gm,newData');
        
        %%% cluster based on gaussian fit for compensatory points
        %%% probability distribution from fit is define by pfit (from above)
        %%% see where new data lands on this (using interp1) then use
        %%% random probability to assign each one
        
%         prob = interp1(nBins,pfit,dGaze);
%         clust1 = rand(size(prob))<prob;
%         clust1(isnan(dGaze))=1;
%         idx2=~clust1;
%         clust = 1;
        
    %%% threshold on gaze velocity
     dGazeV = dGaze*frRate;
        idx2 = abs(dGazeV)>=gzThresh;
        clust = 1;
        
        
        propCVid(i,2,1) = length(find(~app' & idx2==clust))./(length(newData(~app))); %non-comp.
        propCVid(i,2,2) = length(find(app' & idx2==clust))./(length(newData(app)));
        propCVid(i,1,1)=  length(find(~app' & idx2~=clust))./(length(newData(~app))); %compensatory, non-approach
        propCVid(i,1,2)= length(find(app' & idx2~=clust))./(length(newData(app))); %compensatory, approach
        
        %         subplot(6,1,2)
        %         plot(roll); hold on;
        %         plot(tilt);
        %         plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        %         legend('roll','tilt');
        
        %         subplot(6,1,3);
        
        
        %%% get saccade starts
        sacc = idx2==clust;
        saccStart = [0; diff(sacc)>0];
        saccs = find(saccStart);
        clear saccHead saccEye saccApp saccVid saccAz saccTh saccEyeRaw timetoApp
        
        %%% grab traces around saccades;
        if deInter
            buffer = 20;
        else
            buffer = 10;
        end
        ns =0;
        
        
        azAll = [azAll; azdeg'];
        eyeAll = [eyeAll; mnEyeTh'];
        g3All = [g3All; g3];
        appAll = [appAll; [app 0]' ];
        spAll = [spAll; mouseSp{vid}'];
        distAll = [distAll; dist{vid}'];
        
        if length(mouseSp{vid}) ~= length(g3)
            keyboard
        end
        
        
        for i = 1:length(saccs);
            
            if saccs(i)>buffer & saccs(i)<length(drth)-buffer
                ns = ns+1;
                saccHead(:,ns) = cumsum(g3(saccs(i) + (-buffer:buffer)));
                saccEye(:,ns) = mnEyeTh(saccs(i) + (-buffer:buffer));
                saccAz(:,ns) = azdeg(saccs(i) + (-buffer:buffer));
                saccTh(:,ns) = hth(saccs(i) + (-buffer:buffer));
                saccEyeRaw(:,ns) = saccEye(:,ns); %%% copy that is not inverted
                if dGaze(saccs(i))<0
                    saccHead(:,ns) = - saccHead(:,ns);
                    saccEye(:,ns) = - saccEye(:,ns);
                end
                saccApp(ns) = app(saccs(i));
                saccVid(ns) = vid;
                distToApp = saccs(i)-find(app);
                [d nearApp] = min(abs(distToApp));
                if isempty(nearApp)
                    timetoApp(ns)=NaN;
                else
                    timetoApp(ns) = distToApp(nearApp);
                end
            end
        end
        if exist('saccEye','var')
            saccEyeAll = [saccEyeAll saccEye];
            saccHeadAll = [saccHeadAll saccHead];
            saccAppAll = [saccAppAll saccApp];
            saccVidAll = [saccVidAll saccVid];
            saccAzAll = [saccAzAll saccAz];
            saccThAll = [saccThAll saccTh];
            saccEyeRawAll = [saccEyeRawAll saccEyeRaw];
            timetoAppAll = [timetoAppAll timetoApp];
        end
        
        %%% calculate stabilization;
        saccEnds = find(diff(sacc)<0)+1;
        if length(saccEnds)>0 & length(saccs>0) &  (saccEnds(1) < saccs(1)) %%% sometimes trace starts with a saccade end, messes things up
            saccEnds = saccEnds(2:end);
        end
        
        for s = 1:length(saccs)-1;
            stable = saccEnds(s):saccs(s+1)-1;
            allS = allS+1;
            dheadStable{allS} = g3(stable);
            dgzStable{allS} = dGaze(stable);
            eyeStable{allS} = mnEyeTh(stable);
            %             if ~isnan(sum(dgzStable{allS})) & std(cumsum(dgzStable{allS}))>25
            %                 keyboard
            %             end
        end
        
        
        
        propCVid(i,2,1) = length(find(~app' & idx2==clust))./(length(newData(~app))); %non-comp.
        propCVid(i,2,2) = length(find(app' & idx2==clust))./(length(newData(app)));
        propCVid(i,1,1)=  length(find(~app' & idx2~=clust))./(length(newData(~app))); %compensatory, non-approach
        propCVid(i,1,2)= length(find(app' & idx2~=clust))./(length(newData(app))); %compensatory, approach
        
        %         subplot(6,1,2)
        %         plot(roll); hold on;
        %         plot(tilt);
        %         plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        %         legend('roll','tilt');
        
        %         subplot(6,1,3);
        
        
        figure('Renderer', 'painters', 'Position', [100 -100 600 900])
        subplot(5,1,1)
        plot(0.5*( rth(appRange) +lth(appRange))','k','LineWidth',2); hold on;
        plot(rth(appRange)','r'); hold on; plot(lth(appRange)','b');
        plot(find(idx2(appRange)==clust),0.5*(rth(appRange(idx2(appRange)==clust)) +lth(appRange(idx2(appRange)==clust)))','og')
        ylim([-30 30]); % xlim([0 max(length(appRange),1)]);
        
        plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        legend('mean','right','left');
        title(sprintf('vid %d',vid));
        
        subplot(5,1,2)
        plot(diff(rth(appRange)),'r');hold on; plot(diff(lth(appRange)),'b');
        plot(diff(0.5*( rth(appRange) +lth(appRange)))','k','LineWidth',2); hold on;
        
        
        subplot(5,1,3)
        
        plot(azdeg(appRange)); hold on
        %         plot(0.5*( drth(appRange) +dlth(appRange))','k','LineWidth',2); hold on;
        %         plot(drth(appRange)','r'); hold on; plot(dlth(appRange)','b');
        %         plot(find(idx2(appRange)==clust),0.5*(drth(appRange(idx2(appRange)==clust)) +dlth(appRange(idx2(appRange)==clust)))','og')
        %         ylim([-20 20]); % xlim([0 max(length(appRange),1)]);
        plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        title('az to cricket')
        xlim([0 max(length(appRange),1)]);
        
        gzApp = hthApp' +0.5*( rth(appRange) +lth(appRange))';
        % subplot(6,1,6);
        subplot(5,1,4)
        hold on
        plot(hthApp,'Color',[0 0.75 0],'LineWidth',2);
        %plot(find(idx2(appRange)==clust),hthApp(idx2(appRange)==clust),'bo');
        
        plot((gzApp'+30),'k','LineWidth',2);
        sacc = find(idx2(appRange)==clust);
        sacc = sacc(sacc<length(gzApp)-1); %%% make sure we don't run off the end of the data
        for s = 1:length(sacc)
            plot(sacc(s):sacc(s)+1,gzApp(sacc(s):sacc(s)+1) + 30,'r','LineWidth',2)
        end
        
        if appOffset>0
            plot([appOffset appOffset],[-60 60],'g');end
        plot([endOffset endOffset],[-60 60],'r');
        if ~isempty(min(hthApp)) & sum(~isnan(hthApp))>10
            ylim([min(hthApp)-20 max(gzApp)+30]);
        else  ylim([-60 60]);end
        % xlim([0 max(length(appRange),1)]);
        xlabel('frames'); ylabel('deg');
        % legend('gaze','head')
        
        %         subplot(6,1,4)
        %         plot(roll(appRange)); hold on;
        %         plot(tilt(appRange));
        %         plot(find(idx2(appRange)==clust),roll(appRange(idx2(appRange)==clust)),'og')
        %         plot(find(idx2(appRange)==clust),tilt(appRange(idx2(appRange)==clust)),'og')
        %         plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        %         plot([1  max(length(appRange),1)],[0 0],'--')
        %         xlim([0 max(length(appRange),1)]);
        %         if ~isempty(appRange)
        %         ylim([((min(min(roll(appRange),tilt(appRange))))-5) ((max(max(roll(appRange),tilt(appRange))))+5)]);
        %         end
        %         legend('roll','tilt');
        %
        %         subplot(6,1,5)
        %         vgPhi=rphi(appRange)-lphi(appRange);
        %         plot(vgPhi); hold on
        %         plot(roll(appRange),'k'); hold on
        %         plot(find(idx2(appRange)==clust),vgPhi((idx2(appRange)==clust)),'go')
        %         plot(find(idx2(appRange)==clust),roll(appRange(idx2(appRange)==clust)),'go')
        %
        %         plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        %         xlim([0 max(length(appRange),1)]);
        %         legend('R phi- L phi','acc roll')
        %             if ~isempty(appRange)
        %         ylim([((min(min(roll(appRange),vgPhi)))-5) ((max(max(roll(appRange),vgPhi)))+5)]);
        %         end
        
        subplot(5,1,5)
        gzApp = hthApp' +0.5*( rth(appRange) +lth(appRange))';
        plot((gzApp'),'k','LineWidth',2); hold on;
        sacc = find(idx2(appRange)==clust);
        sacc = sacc(sacc<length(gzApp)-1); %%% make sure we don't run off the end of the data
        for s = 1:length(sacc)
            plot(sacc(s):sacc(s)+1,gzApp(sacc(s):sacc(s)+1),'r','LineWidth',2)
        end
        
        
        % hold on
        % plot(dhthApp','Color',[0 0.75 0],'LineWidth',2); hold on
        % plot(dhthApp' +0.5*( drth(appRange) +dlth(appRange))','k','LineWidth',.75);
        %         plot(find(idx2(appRange)==clust),dhthApp(idx2(appRange)==clust),'bo')
        % ylim([-50 50])
        % title('d head, dgaze')
        
        if appOffset>0
            plot([appOffset appOffset],[-60 60],'g');end
        plot([endOffset endOffset],[-60 60],'r');
        if ~isempty(min(hthApp))
            ylim([min(gzApp)-5 max(gzApp)+5]);
            %     ylim([-50 50])
        else  ylim([-50 50]);end
        
        %         subplot(6,1,4)
        %         plot(roll(appRange)); hold on;
        %         plot(tilt(appRange));
        %         plot(find(idx2(appRange)==clust),roll(appRange(idx2(appRange)==clust)),'og')
        %         plot(find(idx2(appRange)==clust),tilt(appRange(idx2(appRange)==clust)),'og')
        %         plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        %         plot([1  max(length(appRange),1)],[0 0],'--')
        %         xlim([0 max(length(appRange),1)]);
        %         if ~isempty(appRange)
        %         ylim([((min(min(roll(appRange),tilt(appRange))))-5) ((max(max(roll(appRange),tilt(appRange))))+5)]);
        %         end
        %         legend('roll','tilt');
        %
        %         subplot(6,1,5)
        %         vgPhi=rphi(appRange)-lphi(appRange);
        %         plot(vgPhi); hold on
        %         plot(roll(appRange),'k'); hold on
        %         plot(find(idx2(appRange)==clust),vgPhi((idx2(appRange)==clust)),'go')
        %         plot(find(idx2(appRange)==clust),roll(appRange(idx2(appRange)==clust)),'go')
        %
        %         plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        %         xlim([0 max(length(appRange),1)]);
        %         legend('R phi- L phi','acc roll')
        %             if ~isempty(appRange)
        %         ylim([((min(min(roll(appRange),vgPhi)))-5) ((max(max(roll(appRange),vgPhi)))+5)]);
        %         end
        
        %         subplot(4,1,4)
        %         hold on
        %         plot(dhthApp','Color',[0 0.75 0],'LineWidth',2); hold on
        %         plot(dhthApp' +0.5*( drth(appRange) +dlth(appRange))','k','LineWidth',.75);
        %         plot(find(idx2(appRange)==clust),dhthApp(idx2(appRange)==clust),'bo')
        %         ylim([-50 50])
        %         title('d head, dgaze')
        %
        %         if appOffset>0
        %             plot([appOffset appOffset],[-60 60],'g');end
        %         plot([endOffset endOffset],[-60 60],'r');
        %         if ~isempty(min(hthApp))
        %             %ylim([min(hthApp)-20 max(hthApp)+20]);
        %             ylim([-50 50])
        %         else  ylim([-50 50]);end
        
        
        %        if exist('gm','var')
        %            X = [dth(appRange)'; [diff(mnEyeApp) 0]; mnEyeApp]';
        %         idx = cluster(gm,X);
        %         sacc = find(idx==3);
        %         for i = 1:length(sacc)-1;
        %             plot(sacc(i):sacc(i)+1,gzApp(sacc(i):sacc(i)+1),'r')
        %         end
    end
    
    
    %         appGaze =(hthApp +0.5*( rth(appRange) +lth(appRange))')
    %         hold on;plot(find(resetPt(appRange)),appGaze(resetPt(appRange)),'ob');
    %         plot(find(stablePt(appRange)),appGaze(stablePt(appRange)),'og')
    %         plot(find(headTurnPt(appRange)),appGaze(headTurnPt(appRange)),'or')
    %
    %         clear head eyeTrace useT
    %         useTAll=(find(idx2(appRange)==clust));
    %         useT=useTAll-1;
    %         %head=g3(find(idx2(appRange)==clust));
    %         head =(.5*(rth+lth));head = head(find(idx2(appRange)==clust));
    %         eyeTrace=(.5*(drth+dlth));
    %
    %         if ~isempty(useT) && useT(1)==0|useT(1)<0
    %         useT=useT(2:end);
    %         head=head(2:end);
    %         end
    %         plot(head,eyeTrace(useT),'bo'); hold on;
    drawnow
    
    
    if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    
    %           clear corr lags
    % %           if ~isempty(useT)
    % %               [corr lags]=nanxcorr(head, eyeTrace(useT),frRate,'coeff');
    % %               corrVid(i,:)=corr; lagsVid(i,:)=lags;
    % %           end
    %     headAll=[headAll head'];
    %     eyeAll=[eyeAll eyeTrace(useT)'];
    
    
    % notnan = ~isnan(sum(saccEyeAll,1)) & ~isnan(sum(saccHeadAll,1));
    % saccHeadAll = saccHeadAll(:,notnan);
    % saccEyeAll = saccEyeAll(:,notnan);
    % saccAppAll = saccAppAll(notnan);
    % saccVidAll = saccVidAll(notnan);
    
    % useTrace = ceil(rand(200,1)*size(saccHeadAll,2));
    % figure
    % plot(-buffer:buffer,saccHeadAll(:,useTrace))
    % hold on
    % plot(-buffer:buffer,nanmean(saccHeadAll,2),'g','LineWidth',2)
    % plot(-buffer:buffer,nanmean(saccEyeAll,2),'r','LineWidth',2);
    % ylabel('head position')
    %     if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    
    
    
    %        if exist('gm','var')
    %            X = [dth(appRange)'; [diff(mnEyeApp) 0]; mnEyeApp]';
    %         idx = cluster(gm,X);
    %         sacc = find(idx==3);
    %         for i = 1:length(sacc)-1;
    %             plot(sacc(i):sacc(i)+1,gzApp(sacc(i):sacc(i)+1),'r')
    %         end
    
    
    
    %         appGaze =(hthApp +0.5*( rth(appRange) +lth(appRange))')
    %         hold on;plot(find(resetPt(appRange)),appGaze(resetPt(appRange)),'ob');
    %         plot(find(stablePt(appRange)),appGaze(stablePt(appRange)),'og')
    %         plot(find(headTurnPt(appRange)),appGaze(headTurnPt(appRange)),'or')
    %
    %         clear head eyeTrace useT
    %         useTAll=(find(idx2(appRange)==clust));
    %         useT=useTAll-1;
    %         %head=g3(find(idx2(appRange)==clust));
    %         head =(.5*(rth+lth));head = head(find(idx2(appRange)==clust));
    %         eyeTrace=(.5*(drth+dlth));
    %
    %         if ~isempty(useT) && useT(1)==0|useT(1)<0
    %         useT=useT(2:end);
    %         head=head(2:end);
    %         end
    %         plot(head,eyeTrace(useT),'bo'); hold on;
    drawnow
    
    figure
    plot(hthApp-11);
    hold on
    plot(-0.5*(rth(appRange)+lth(appRange)));
    legend('head','-eye');
    
 figure
    plot(-60:60,nanxcorr(hthApp,0.5*(rth(appRange)+lth(appRange)),60,'coeff'))
    xlim([-30 30]);
    ylim([-1 0.5]);
   hold on; plot([0 0], [-1 0.5],'r')
    
   x = 1:length(Rtheta{vid})-1;
   figure
   plot(x,accelChannels{vid}(1:end-1,6)');
   hold on
   plot(x,-diff(0.5*(Rtheta{vid} + Ltheta{vid})))
   legend('acc3','diff mnEye')


      figure
      plot(accelChannels{vid}(2:end,6)',diff(0.5*(Rtheta{vid} + Ltheta{vid})),'.')
      axis equal
   
    figure
    plot(-60:60,nanxcorr(accelChannels{vid}(1:end-1,6)',diff(0.5*(Rtheta{vid} + Ltheta{vid})),60,'coeff'))
    xlim([-30 30]);
    ylim([-1 0.5]);
   hold on; plot([0 0], [-1 0.5],'r')
   
end

%% analyze eye/head position around saccades
%%% defining 3 clusters

%%% plot all eye traces (select every 100th one)
figure
plot(-buffer:buffer,saccEyeAll(:,1:100:end))
hold on
plot(-buffer:buffer,nanmean(saccEyeAll,2),'g','LineWidth',2);
ylabel('mean Eye position'); ylim([-20 20])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%% can either select all saccs for analysis, or just approaches.
%%% right now, all saccs
useS = saccAppAll>=0;
eyeData = saccEyeAll(:,useS);
headData = saccHeadAll(:,useS);


%%% do kmeans on mean-subtracted traces
clear eyeDataCent
for i = 1:size(eyeData,2);
    eyeDataCent(:,i) = eyeData(:,i)-nanmean(eyeData(:,i),1);
end
idx = kmeans(eyeDataCent(5:20,:)',4);
%idx(idx==4)=1;  %%% group together clust 4 and clust 1, since both are "tracking", large and small

%%% tried used gmm, doesn't work well
% gm = fitgmdist(eyeDataCent(5:end,:)',4,'Replicates',10);
% idx = cluster(gm,eyeDataCent(5:end,:)');


%%% plotting 4 clusters
eyeFig = figure;
headFig = figure;
for i = 1:4
    
    figure(eyeFig);
    use = (idx==i);
    subplot(2,2,i);
    % plot(saccEyeAll(:,use));
    hold on; plot(nanmean(eyeData(:,use),2),'g', 'Linewidth',2);
    ylim([-10 10]); xlim([5 20]);title(sprintf('n = %d',sum(use)))
    ylabel('eye theta')
    
    if i ==4
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(headFig)
    
    subplot(2,2,i);
    % plot(saccHeadAll(:,use));
    hold on; plot(nanmean(headData(:,use),2)-nanmean(headData(5,use),2) ,'g', 'Linewidth',2);
    ylim([0 40]);xlim([5 20])%ylim([-20 20]); title(sprintf('n = %d',sum(use)))
    ylabel('head theta')
    
    if i ==4
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
end


%%% find frequencies of each saccade type during app, non-app
for i = 1:4
    sfreq(i,1) = sum(idx==i & saccAppAll'==0);
    sfreq(i,2) = sum(idx==i & saccAppAll'==1);
end

%% calculate stability and frequency of fixations

%%% remove stabilizations with NaNs (since these may cause to miss a saccade)
for i = 1:length(eyeStable)
    nanCount(i) = sum(isnan(eyeStable{i})) + sum(isnan(dheadStable{i}));
end
sum(nanCount==0)/sum(nanCount>=0)
good = find(nanCount==0);
for i = 1:length(good)
    dheadStableGood{i} = dheadStable{good(i)};
    eyeStableGood{i} = eyeStable{good(i)};
    dgzStableGood{i} = dgzStable{good(i)};
end

%%% calculate stats for each fixation (duration and stdev of positions)
clear stabDur headStd gazeStd eyeStd headMn gazeMn
for i = 1:length(dheadStableGood);
    stabDurGood(i) = length(dheadStableGood{i});
    if length(dheadStableGood{i})>1   %%% can't do stats on one point
        headStd(i) = std(cumsum(dheadStableGood{i}));
        gazeStd(i) = std(cumsum(dgzStableGood{i}));
        eyeStd(i) = std(eyeStableGood{i});
        headMn(i)= nanmean((dheadStableGood{i}));
        gazeMn(i) = nanmean((dgzStableGood{i}));
    else
        headStd(i) = NaN;
        gazeStd(i) = NaN;
        eyeStd(i) = NaN;
        headMn(i)= NaN;
        gazeMn(i)= NaN;
    end
end

%%% plot histogram of durations
hbins = 1:2:150;
figure
h =hist(stabDurGood,hbins);  %%% add error bars here!
plot(hbins/frRate,h/sum(h)); xlabel('secs'); ylabel('fraction');
title(sprintf('fixation duration median = %0.2f +/- %0.2f sec',nanmedian(stabDurGood)/frRate,(1/frRate)*std(stabDurGood)/sqrt(length(stabDurGood)))); xlim([0 2]);

%%% plot histogram of durations (with error bars)
hbins = 1:2:150;
figure
shadedErrorBar(hbins/30,h/sum(h),sqrt(h)/sum(h)); xlabel('secs'); ylabel('fraction');
hold on; plot([nanmedian(stabDurGood)/frRate nanmedian(stabDurGood)/frRate], [0 .14])
title(sprintf('fixation duration median = %0.2f +/- %0.2f sec',nanmedian(stabDurGood)/frRate,(1/frRate)*std(stabDurGood)/sqrt(length(stabDurGood)))); xlim([0 2]);



%%% histogram of stability for head and gaze
figure
hbins = 0.25:0.5:20;
hold on
h = hist(headStd,hbins);
plot(hbins,h/sum(h));
h = hist(gazeStd,hbins);
plot(hbins,h/sum(h));
legend('head','gaze');
xlabel('RMS stabilization (deg)'); ylabel('fraction'); xlim([0 15])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%% histogram of stability for head and gaze (with error bars)
figure
hbins = 0.25:0.5:20;
hold on
h = hist(headStd,hbins);  %%% add error bars here!
shadedErrorBar(hbins,h/sum(h),sqrt(h)/sum(h),'b');
h = hist(gazeStd,hbins);
shadedErrorBar(hbins,h/sum(h),sqrt(h)/sum(h),'k');
legend('head','gaze');
xlabel('RMS stabilization (deg)'); ylabel('fraction'); xlim([0 15])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

[p,h,stats]=ranksum(headStd(~isnan(headStd)),gazeStd(~isnan(gazeStd)),'method','approximate')

[h p]=ttest2(headStd(~isnan(headStd)),gazeStd(~isnan(gazeStd)))

%%% calculate average stability of head and gaze
stability(1) = nanmedian(headStd);
stability(2) = nanmedian(gazeStd);
stabErr(1) = nanstd(headStd)/sqrt(sum(~isnan(headStd)));
stabErr(2) = nanstd(gazeStd)/sqrt(sum(~isnan(gazeStd)));
figure
bar(1:2, stability);
hold on
errorbar(1:2,stability, stabErr,'.');
set(gca,'Xtick',1:2); set(gca,'XtickLabel',{'head','gaze'}); ylabel('RMS stability (deg)'); title('stability during fixations')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%% look at the longest stabilization (out of interest only, maybe demo
%%% that no saccades during head stability
[m longStable] = max(stabDurGood);
figure
plot(cumsum(dheadStableGood{longStable}- median(dheadStableGood{longStable}))); %%% subtract median since the accelerometer drifts
hold on
plot(eyeStableGood{longStable})
plot(cumsum(dgzStableGood{longStable}- median(dheadStableGood{longStable})))
title(sprintf('head %0.1f  gaze %0.1f',headStd(longStable),gazeStd(longStable)));
legend('head','eyes','gaze')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

% %%% plot a random selection of stabilizations
% for stable = 1:100:length(stabDur)
% figure
% plot(cumsum(dheadStableGood{stable}));
% hold on
% %plot(eyeStableGood{longStable})
% plot(cumsum(dgzStableGood{stable}))
% title(sprintf('head %0.1f  gaze %0.1f',headStd(stable),gazeStd(stable)));
% end

apps = find(saccAppAll==1);

s = saccAzAll(:,apps);
azOffset=nanmedian(s(:))

if nanmean(azAll)<=0
    azAll = azAll-azOffset;  %%% don't do it twice
end
eyeAzAll = azAll+eyeAll;

hbins = -190:2:190;
figure
h = hist(azAll,hbins)
plot(hbins,h/sum(h));
hold on
h = hist(eyeAzAll,hbins);
plot(hbins,h/sum(h));
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


use = abs(azAll)<90 & abs(g3All)<1;
nanmean(abs(azAll(use)))
nanmean(abs(eyeAzAll(use)))

figure
plot(azAll(1:10:end),eyeAll(1:10:end),'.')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


figure
plot(azAll,spAll,'.')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


figure
hist(spAll,1:500);
xlim([0 50])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

spAll(spAll>50)= NaN;

range = 1:10:length(azAll);
spSmooth = medfilt1(spAll,5);
azSmooth = medfilt1(azAll,5);
azData = azSmooth(range); spData = spSmooth(range); appData = appAll(range);
distData = distAll(range);


figure
subplot(2,2,1)
plot(azData(spData<10),distData(spData<10),'k.'); xlabel('azimuth'); ylabel('distance'); title('stationary')
hold on
% plot(azData(spData<10 & appData ==1),distData(spData<10 & appData ==1),'g.'); xlabel('azimuth'); ylabel('distance');

subplot(2,2,2)
plot(azData(spData>10),distData(spData>10),'k.'); xlabel('azimuth'); ylabel('distance'); title('moving')
hold on
%plot(azData(spData>10 & appData ==1),distData(spData>10 & appData ==1),'g.'); xlabel('azimuth'); ylabel('distance');


subplot(2,2,3)
plot(azData(spData<10),distData(spData<10),'b.'); xlabel('azimuth'); ylabel('distance'); title('stationary')
hold on
plot(azData(spData<10 & appData ==1),distData(spData<10 & appData ==1),'g.'); xlabel('azimuth'); ylabel('distance');

subplot(2,2,4)
plot(azData(spData>10),distData(spData>10),'b.'); xlabel('azimuth'); ylabel('distance'); title('moving')
hold on
plot(azData(spData>10 & appData ==1),distData(spData>10 & appData ==1),'g.'); xlabel('azimuth'); ylabel('distance');


if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


figure
subplot(1,2,1)
plot(azData,spData,'.')
axis([-150 150 0 40])
subplot(1,2,2)
plot(azData,spData,'.')
axis([-150 150 0 40])
hold on
plot(azData(appData==1),spData(appData==1),'g.')

newApp = abs(azData)<30 & spData>10;
newApp = medfilt1(double(newApp),5);

figure
subplot(2,1,1)
plot(azData(spData<10),distData(spData<10),'b.'); xlabel('azimuth'); ylabel('distance'); title('stationary')
hold on
plot(azData(spData<10 & newApp ==1),distData(spData<10 & newApp ==1),'g.'); xlabel('azimuth'); ylabel('distance');

subplot(2,1,2)
plot(azData(spData>10),distData(spData>10),'b.'); xlabel('azimuth'); ylabel('distance'); title('moving')
hold on
plot(azData(spData>10 & newApp ==1),distData(spData>10 & newApp ==1),'g.'); xlabel('azimuth'); ylabel('distance');

figure
plot(newApp(1:1000))



%% analyze eye and head relative to cricket for saccades during approaches
apps = find(saccAppAll==1);

%%% there is ~6deg asymmetry in head position relative to cricket. this is
%%% probably based on the definition of theta from "model" head points.
s = saccAzAll(:,apps);
azOffset=nanmedian(s(:))
saccAzAllRaw = saccAzAll; %%% save raw value
saccAzAll = saccAzAll - azOffset;

%%% get eyes relative to cricket
eyeAz = saccAzAll + saccEyeRawAll;

%%% plot azimuth for all approaches
figure
plot(abs(saccAzAll(:,apps)));
hold on
plot(nanmean(abs(saccAzAll(:,apps)),2),'g','Linewidth',2); title('head azimuth for all approaches')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

figure
plot(abs(eyeAz(:,apps)));
hold on
plot(nanmean(abs(eyeAz(:,apps)),2),'g','Linewidth',2)
title('gaze azimuth for all approaches')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

figure
plot(nanmedian(abs(saccAzAll(:,apps)),2),'b','Linewidth',2)
hold on
plot(nanmedian(abs(eyeAz(:,apps)),2),'g','Linewidth',2)
legend('head azimuth','gaze azimuth')
if deInter
    xlim([15 30]);
end
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
%%% head and eye relative to cricket (azimuth) on all saccades
figure
if deInter
    trange=3:34
else
    trange=3:17;
end
hold on
data = nanmedian(abs(eyeAz(:,apps)),2);
err = nanstd(abs(eyeAz(:,apps)),[],2) ./ sqrt(sum(~isnan(eyeAz(:,apps)),2))
data = data(trange); err = err(trange);
t = (0:length(data)-1)/frRate;
shadedErrorBar(t,data,err,'k')
data = nanmedian(abs(saccAzAll(:,apps)),2);
err = nanstd(abs(saccAzAll(:,apps)),[],2) ./ sqrt(sum(~isnan(saccAzAll(:,apps)),2))
data = data(trange); err = err(trange);
t = (0:length(data)-1)/frRate;
shadedErrorBar(t,data,err,'b')
legend('head azimuth','gaze azimuth');
ylim([7.5 27.5]);
xlim([t(1)-1/(2*frRate) t(end)+1/(2*frRate)]);
plot([t(19) t(19)],[7 28])
axis square
xlabel('secs');
ylabel('azimuth to cricket (deg)')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
%%% plot for all 4 clusters
figure
for i =1:4
    subplot(2,2,i)
    use = saccAppAll'==1 & idx ==i;
    plot(nanmedian(abs(saccAzAll(:,use)),2),'r','Linewidth',2)
    hold on
    plot(nanmedian(abs(eyeAz(:,use)),2),'g','Linewidth',2)
    ylim([0 40])
    title(sprintf('clust %d',i))
    legend('head','eye');
    ylabel('azimuth')
    %%% add mean of all clusters
end

if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%% where does eye end up after saccade?
if deInter
    preSaccT = 21; postSaccT=23;
else
    preSaccT = 11; postSaccT=12;  %%% timepoints to compare pre/post saccade
end
bins = -(2*frRate):10:(2*frRate);

%%% get histograms of eye&head for pre/post saccade
clear eyeAzHist
eyeAzHist(:,1) = hist(eyeAz(preSaccT,apps),bins);
eyeAzHist(:,2) = hist(eyeAz(postSaccT,apps),bins);

AzHist(:,1) = hist(saccAzAll(preSaccT,apps),bins);
AzHist(:,2) = hist(saccAzAll(postSaccT,apps),bins);

%%% plot pre
figure
plot(bins,eyeAzHist(:,1))
hold on
plot(bins,AzHist(:,1))
xlabel('azimuth'); legend('gaze','head')
title('pre saccade')
xlim([-60 60])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%% plot post
figure
plot(bins,eyeAzHist(:,2))
hold on
plot(bins,AzHist(:,2))
xlabel('azimuth'); legend('gaze','head')
title('post saccade')
xlim([-60 60])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%% plot pre (with error bars)
figure
shadedErrorBar(bins,eyeAzHist(:,1)/sum(eyeAzHist(:,1)), sqrt(eyeAzHist(:,1))/sum(eyeAzHist(:,1)),'b')
hold on
shadedErrorBar(bins,AzHist(:,1)/sum(AzHist(:,1)), sqrt(AzHist(:,1))/sum(AzHist(:,1)),'r')
xlabel('azimuth (deg)'); legend('gaze','head')
title('pre saccade')
xlim([-60 60])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%% plot post (with error bars)
figure
shadedErrorBar(bins,eyeAzHist(:,2)/sum(eyeAzHist(:,2)), sqrt(eyeAzHist(:,2))/sum(eyeAzHist(:,2)),'b')
hold on
shadedErrorBar(bins,AzHist(:,2)/sum(AzHist(:,2)), sqrt(AzHist(:,2))/sum(AzHist(:,2)),'r')
xlabel('azimuth (deg)'); legend('gaze','head')
title('post saccade')
xlim([-60 60])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end




%%% get median offset of eye&head for pre/post saccade
eyeOffset(1) = nanmedian(abs(eyeAz(preSaccT,apps)));
eyeOffset(2) = nanmedian(abs(eyeAz(postSaccT,apps)))
eyeOffsetErr(1) = nanstd(abs(eyeAz(preSaccT,apps)))/sqrt(sum(~isnan(eyeAz(preSaccT,apps))));
eyeOffsetErr(2) = nanstd(abs(eyeAz(postSaccT,apps)))/sqrt(sum(~isnan(eyeAz(postSaccT,apps))))

headOffset(1) = nanmedian(abs(saccAzAll(preSaccT,apps)));
headOffset(2) = nanmedian(abs(saccAzAll(postSaccT,apps)));
headOffsetErr(1) = nanstd(abs(saccAzAll(preSaccT,apps)))/sqrt(sum(~isnan(saccAzAll(preSaccT,apps))));
headOffsetErr(2) = nanstd(abs(saccAzAll(postSaccT,apps)))/sqrt(sum(~isnan(saccAzAll(postSaccT,apps))))

eyez=nanmedian(abs(eyeAz(preSaccT,apps)),1);
headz=nanmedian(abs(saccAzAll(preSaccT,apps)),1);
[h p]=ttest2(eyez,headz)

eyez=nanmedian(abs(eyeAz(postSaccT,apps)),1);
headz=nanmedian(abs(saccAzAll(postSaccT,apps)),1);
[h p]=ttest2(eyez,headz)

%%% plot errorbar
figure
barweb([headOffset; eyeOffset]',[headOffsetErr; eyeOffsetErr]')
ylabel('azimuth to cricket (deg)');
set(gca,'Xticklabel',{'pre saccade','post saccade'})
legend('head','gaze')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
for vid=1:length(useData)
    clear hT azC use useN corrR corrRA corrL corrLA
    nframe = min(length(az{useData(vid)}),length(accelChannels{useData(vid)}(:,6)));
    
    nframe = min(nframe,length(appEpoch{vid}));
    useN=appEpoch{vid}(1:nframe)==0; use=appEpoch{vid}(1:nframe)==1;
    
    azC=-(az{useData(vid)});
    % hT=d_Theta{useData(vid)}; hT=hT-nanmedian(hT);
    
    hT=(accelChannels{useData(vid)}(:,6))-gyroBias;
    if sum(useN)>3 & sum(~isnan(hT(useN)))>20 & sum(~isnan(azC(useN)))>20
        [corrR lagsR]= nanxcorr(azC(useN),hT(useN),frRate,'coeff');
        hold on;
    else
        corrR=NaN;
    end
    
    %     subplot(4,5,vid)
    %     plot(rad2deg(azC(use))); hold on; plot(hT(use));
    if sum(use)>4 & sum(~isnan(hT(use)))>20
        [corrRA lagsRA]= nanxcorr(azC(use),hT(use),frRate,'coeff');
    else
        corrRA=NaN;
    end
    corrRAll(vid,:)=corrR;
    corrRAAll(vid,:)=corrRA;
end

figure%('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); xlim([1 2*frRate+1]); axis square; ylim([-.2 .5]);

clear L
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'g-');
legend(L,{'az to cricket','dtheta head'});

title('Figure 4H: Corr of az and dHead Theta, head is directed to cricket')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%


clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll mnEye
for vid=1:length(useData)
    clear hT azC use useN corrR corrRA corrL corrLA
    nframe = min(length(az{useData(vid)}),length(accelChannels{useData(vid)}(:,6)));
    
    nframe = min(nframe,length(appEpoch{vid}));
    useN=appEpoch{vid}(1:nframe)==0; use=appEpoch{vid}(1:nframe)==1;
    
    azC=-(az{useData(vid)});
    hT=(accelChannels{useData(vid)}(:,6))-gyroBias;
    
    mnEye=.5*(Rtheta{useData(vid)}+Ltheta{useData(vid)}); mnEye=mnEye-nanmedian(mnEye);
    hT=(accelChannels{useData(vid)}(:,6))-gyroBias;
    dGazeVid= diff(hT+mnEye);
    
    if sum(useN)>3 & sum(~isnan(dGazeVid(useN)))>20 & sum(~isnan(azC(useN)))>20
        [corrR lagsR]= nanxcorr(azC(useN),dGazeVid(useN),frRate,'coeff');
        hold on;
    else
        corrR=NaN;
    end
    
    %     subplot(4,5,vid)
    %     plot(rad2deg(azC(use))); hold on; plot(hT(use));
    if sum(use)>4 & sum(~isnan(dGazeVid(use)))>20
        [corrRA lagsRA]= nanxcorr(azC(use),dGazeVid(use),frRate,'coeff');
    else
        corrRA=NaN;
    end
    corrRAll(vid,:)=corrR;
    corrRAAll(vid,:)=corrRA;
end

figure%('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); xlim([1 2*frRate+1]); axis square; ylim([-.2 .5]);

clear L
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'g-');
legend(L,{'az to cricket','dtheta head'});

title('Figure 4I: Corr of az and dGaze, gaze is less directed to cricket')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%% old ways of trying to cluster sacades

% plot(diff(mean(saccHeadAll,2)))
% plot(mean(saccEyeAll,2))
% plot(diff(mean(saccEyeAll,2)))
% vhead = diff(saccHeadAll,1);
%
% maxpre = max(vhead(1:9,:),[],1);
% meanpre = mean(vhead(1:9,:),1);
% meanpost = mean(saccHeadAll(15:20,:),1);
%
%     stationary_sm = abs(meanpre)<2 &  abs(maxpre)<2  &meanpost< 20;
%     stationary_lg = abs(meanpre)<2 &  abs(maxpre)<2 & meanpost>20;
%
%     stationary = abs(meanpre)<1 &  abs(maxpre)<2;
%     moving = meanpre>1 & abs(maxpre)>5 ;
%     if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% figure
% plot(nanmean(saccEyeAll(:,stationary_sm),2));
% hold on
% plot(nanmean(saccEyeAll(:,stationary_lg),2));
% plot(nanmean(saccEyeAll(:,moving),2));
%     if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% figure
% plot(nanmean(saccHeadAll(:,stationary_sm),2))
% hold on
% plot(nanmean(saccHeadAll(:,stationary_lg),2))
% plot(nanmean(saccHeadAll(:,moving),2))
%     if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
%
% figure
% plot(saccEyeAll(:,stationary_sm))
%     if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% figure
% plot(saccEyeAll(:,stationary_lg))
%     if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% figure
% plot(saccEyeAll(:,moving))
%     if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% figure
% plot(nanmean(saccHeadAll(:,stationary),2))
% hold on
% plot(nanmean(saccHeadAll(:,moving),2))
%     if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
%
% v = unique(saccVidAll)
% %for i = 1:length(v)
%     for i = 1:10
%
%     stationary = abs(meanpre)<2 &  abs(maxpre)<2 & saccVidAll==v(i);
% moving = meanpre>1 & abs(maxpre)>5 & saccVidAll==v(i);
%
% figure
% plot(saccEyeAll(:,stationary));
% hold on
% plot(saccEyeAll(:,moving));
%     if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% end
%
% figure
% plot(nanmean(saccEyeAll,2))
%     if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end



% ylabel('delta eye theta n-1');
% xlabel('eye position')
% xlim([-30 30]);ylim([-30 30]);
%
% errC=nanstd(propCVid)./(sqrt(length(useData)));
% errNC=nanstd(propNCVid)./(sqrt(length(useData)));
% figure; barweb([nanmean(propCVid) nanmean(propNCVid)],[errC errNC]);
% ylim([0 1]);
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% figure
% for vid=1:15
%     subplot(3,5,vid)
%     plot(lagsVid(vid,:),corrVid(vid,:)); axis square
%
% end
% Figure 5B: example of target selection saccades

% Figure 5C: quantification of target selection sacc
%
% clusters(1,:)=sum(idx==1|idx==3)/length(idx);
% clusters(2,:)=sum(idx==2)/length(idx);
%
% clustersErr(1,:) = squeeze(nanstd(idx(idx==1|idx==3)))/sqrt(length(idx));
% clustersErr(2,:) = squeeze(nanstd(idx(idx==2)))/sqrt(length(idx));
% barweb(clusters,clustersErr)
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

% Figure 5D: example of resetting saccades

% Figure 5E: quantification of resetting sacc



if savePDF
    pSname='T:\PreyCaptureAnalysis\Data\';
    
    filen=sprintf('%s','PaperFigs_deinter_052620_halfShift_a_gazeVelocity','.pdf')
    %     filen=sprintf('%s','PaperFigs_011519_c','.pdf')
    
    pdfilename=fullfile(pSname,filen);
    dos(['ps2pdf ' psfilename ' ' pdfilename]);
    delete(psfilename);
else
    pSname='T:\PreyCaptureAnalysis\Data\';
end









