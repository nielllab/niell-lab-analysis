close all; clear all;
%    load('ACCAnalyzed_AllAnimals_010820_noDLS.mat')
load('ACC_AllAnimals_021520_a.mat')
 %load('ACC_deInter_Analyzed_AllAnimals_011520_a.mat')

savePDF=1;
if savePDF
    psfilename = 'C:\analysisPS_B.ps';
    if exist(psfilename,'file')==2;delete(psfilename);end
end
pname = 'T:\PreyCaptureAnalysis\Data\';

%% FIGURE 1: a freely moving eye and head tracking system
load('T:\PreyCaptureAnalysis\Data\ControlAnalysis\compileAllAnimals_CONTROL_090219.mat','mouseSp','appEpoch','thetaHead');

for vid=1:length(thetaHead)
    prop(vid)=sum(~isnan(thetaHead{vid}))./length(thetaHead{vid})
end
useData=find(prop>.85); % 90% of points needs to be present for exp to be used

for vid=1:length(useData)
    appTime=appEpoch{useData(vid)};
    Cntrl_speed(vid,1) =sum(mouseSp{useData(vid)}(:,appTime==0)>5)./(length(mouseSp{useData(vid)}(:,appTime==0)));
    Cntrl_speed(vid,3) =nanmean(mouseSp{useData(vid)}(:,appTime==0));

    if sum(appTime==1)>30
    Cntrl_speed(vid,2) =sum(mouseSp{useData(vid)}(:,appTime==1)>5)./(length(mouseSp{useData(vid)}(:,appTime==1)));
    Cntrl_speed(vid,4) =nanmean(mouseSp{useData(vid)}(:,appTime==1));

    else
        Cntrl_speed(vid,2)=nan;
        Cntrl_speed(vid,4)=nan;

    end
end
Cntrl_speed(Cntrl_speed(:,4)>60,4)=NaN; % noisy control data, speeds at a few timepoints go above 150 cm/sec
Cntrl_err= nanstd(Cntrl_speed)/sqrt(length(useData));

clear mouseSp appEpoch useData appTime thetaHead vid
% load('ACCAnalyzed_AllAnimals_010820_noDLS.mat')
load('ACC_AllAnimals_021520_a.mat')
savePDF=1;
psfilename = 'C:\analysisPS_B.ps';
for vid=1:length(useData)
    appTime=appEpoch{vid};
    Cam_speed(vid,1) = sum(mouseSp{useData(vid)}(:,appTime==0)>5)./(length(mouseSp{useData(vid)}(:,appTime==0)));
    Cam_speed(vid,3) = nanmean(mouseSp{useData(vid)}(:,appTime==0));
    if sum(appTime==1)>30
        Cam_speed(vid,2) =sum(mouseSp{useData(vid)}(:,appTime==1)>5)./(length(mouseSp{useData(vid)}(:,appTime==1)));
        Cam_speed(vid,4) = nanmean(mouseSp{useData(vid)}(:,appTime==1));
    else
        Cam_speed(vid,2)=nan;
        Cam_speed(vid,4)=nan;
    end
end
Cam_err= nanstd(Cam_speed)/sqrt(length(useData));
% figure; barweb([nanmean(Cam_speed(:,1)) nanmean(Cam_speed(:,2)) nanmean(Cntrl_speed(:,1)) nanmean(Cntrl_speed(:,2))],...
%     [nanmean(Cam_err(:,1)) nanmean(Cam_err(:,2)) nanmean(Cntrl_err(:,1)) nanmean(Cntrl_err(:,2))]);
% title('proportion of time moving')
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

[h p]= kstest2(Cntrl_speed(:,1),Cntrl_speed(:,2))
[h p]= kstest2(Cam_speed(:,1),Cam_speed(:,2))
[h p]= kstest2(Cntrl_speed(:,1),Cam_speed(:,1))
[h p]= kstest2(Cntrl_speed(:,2),Cam_speed(:,2))
% 
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
[h , ~]= kstest2(Cam_speed(:,3),Cam_speed(:,4))
[h p]= kstest2(Cam_speed(:,3),Cntrl_speed(:,3))
[h p]= kstest2(Cam_speed(:,4),Cntrl_speed(:,4))

%%%average number of captures per session
exp=[4.6	6.2	5.2	4.333333333	3.166666667	4.166666667	5.2	6.6	3.5	4.166666667]; %avg numbers of crickets for each animal, across sessions, from experiment spreadsheet
cntrl=[4.933333333	6.366666667	5.366666667	5.666666667	5.833333333	4.5	5	6.6	5.166666667	5.333333333];
expM=mean(exp); errExp=std(exp)./(length(sqrt(exp)));
mnC=mean(cntrl); errC=std(cntrl)./(length(sqrt(cntrl)));
comp=[expM mnC];err=[errExp,errC];
figure
barweb(comp,err)
[sig p]=ttest2(exp',cntrl');

% Figure 1G: gyro & DLC traces
figure
range=350:830;
for vid=10%1:10%length(useData)
%     subplot(2,5,vid);
       plot(d_Theta{useData(vid)}(range,:)); 
    hold on;
    plot(accelChannels{useData(vid)}(range,6));
end
legend('DLC head theta','gyroscope')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

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

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
for vid=1:length(useData)
    clear hT azC use useN corrR corrRA corrL corrLA
    nframe = min(length(az{useData(vid)}),length(d_Theta{useData(vid)}));
    
    nframe = min(nframe,length(appEpoch{vid}));
    useN=appEpoch{vid}(1:nframe)==0; use=appEpoch{vid}(1:nframe)==1;

    azC=-(az{useData(vid)}); hT=d_Theta{useData(vid)}; hT=hT-nanmedian(hT);
    
    if sum(useN)>3 & sum(~isnan(hT(useN)))>20 & sum(~isnan(azC(useN)))>20
        [corrR lagsR]= nanxcorr(azC(useN),hT(useN),30,'coeff');
        hold on;
    else
            corrR=NaN;
    end
    
%     subplot(4,5,vid)
%     plot(rad2deg(azC(use))); hold on; plot(hT(use));
   if sum(use)>4 & sum(~isnan(hT(use)))>20
        [corrRA lagsRA]= nanxcorr(azC(use),hT(use),30,'coeff');
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
plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); xlim([1 61]); axis square; ylim([-.2 .5]);

clear L
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'g-');
legend(L,{'az to cricket','dtheta head'});

title('Figure 1I: Corr of az and dHead Theta, head is directed to cricket')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

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
        useN = appEpoch{vid}==0;use = appEpoch{vid}==1;
        subplot(1,2,1)
        plot(tR(useN(1:60:length(useN))),pR(useN(1:60:length(useN))),'.b'); hold on; xlim([-60 60]);ylim([-60 60]); axis square
        plot(tR(use(1:60:length(use))),pR(use(1:60:length(use))),'.g'); hold on; xlim([-60 60]);ylim([-60 60]); axis square
        title('R eye')
        xlabel('yaw (deg)'); ylabel('pitch(deg)');
        subplot(1,2,2)
        plot(tL(useN(1:60:length(useN))),pL(useN(1:60:length(useN))),'.b'); hold on; xlim([-60 60]);ylim([-60 60]); axis square
        plot(tL(use(1:60:length(use))),pL(use(1:60:length(use))),'.g'); hold on; xlim([-60 60]);ylim([-60 60]); axis square
        xlabel('yaw (deg)'); ylabel('pitch(deg)');
        title('L eye')
end
suptitle('Figure 2A: R & L eye yaw & pitch');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end



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
%     plot(tR(use(1:20:end)),tL(use(1:20:end)),'go'); hold on; axis square; %ylim([-30 30]); xlim([-30 30]);
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

figure;subplot(1,2,1); plot(nBins,allThN(:,1),'b'); hold on; plot(nBins,allThA(:,1),'g'); title('Right Eye Theta'); axis square
ylim([0 .25]); xlabel('theta (degrees)'); ylabel('proportion of time');
subplot(1,2,2);plot(nBins,allThN(:,2),'b'); hold on; plot(nBins,allThA(:,2),'g');title('Left Eye Theta');axis square
ylim([0 .25]); xlabel('theta (degrees)'); ylabel('proportion of time');
legend('non approach','approach')
suptitle('Figure 2B: eyes are more centered in yaw during app')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


allThNP(:,1)=nansum(tNcountsP(:,:,1),1)./(length(useData))
allThNP(:,2)=nansum(tNcountsP(:,:,2),1)./(length(useData))

allThAP(:,1)=nansum(tAcountsP(:,:,1),1)./(length(useData))
allThAP(:,2)=nansum(tAcountsP(:,:,2),1)./(length(useData))

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
load('ACC_deInter_Analyzed_AllAnimals_011520_a.mat')

savePDF=1;
    psfilename = 'C:\analysisPS_B.ps';

figure
appTimes=find(appEpoch{40}(1:5400)==1)
plot(Rtheta{useData(40)}(1:5400)); hold on; %3 minute segment
plot(Ltheta{useData(40)}(1:5400)); hold on;
plot(appTimes,20,'og')
legend('r theta','l theta');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


range=3000:4800
figure
plot(Rtheta{useData(40)}(range,:)); hold on; % 1 minute segment
plot(Ltheta{useData(40)}(range,:)); hold on;
legend('r theta','l theta');
% plot(accelChannels{useData(40)}(range,1)); %roll 
% plot(Rtheta{useData(40)}(range,:)-Ltheta{useData(40)}(range,:))
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

range=3500:4400
figure
plot(Rtheta{useData(40)}(range,:)); hold on; % 30 second segment
plot(Ltheta{useData(40)}(range,:)); hold on;
legend('r theta','l theta');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% 
% close all
% range=4000:4600
% for vid = 40%(useData)
%     
%    if length(Rtheta{useData(vid)})>=100
% tr=medfilt1(Rtheta{useData(vid)}(range,:),8);
% pr=medfilt1(Rphi{useData(vid)}(range,:),8);
% 
% tl=medfilt1(Ltheta{useData(vid)}(range,:),8);
% pl=medfilt1(Lphi{useData(vid)}(range,:),8);
% 
% 
% figure(1)
% % subplot(1,2,1)
% plot(tr-nanmean(tr),pr-nanmean(pr),'k'); axis square; hold on; xlabel('eye theta'); ylabel('eye phi'); axis equal
% title('right eye'); xlim([-30 30]); ylim([-30 30]); colormap jet; colorbar
% figure(2)
% % subplot(1,2,2)
% plot(tl-nanmean(tl),pl-nanmean(pl),'k'); axis square; hold on; xlabel('eye theta'); ylabel('eye phi'); axis equal
% title('left eye'); xlim([-30 30]); ylim([-30 30]); colorbar
% 
% for i =1:length(tl)
%     figure(1)
% % subplot(1,2,1)
% plot(tr(i)-nanmean(tr),pr(i)-nanmean(pr),'.','Markersize',15,'Color', cmapVar(i,1,length(tl),jet)); axis square; hold on
% % subplot(1,2,2)
% figure(2)
%  plot(tl(i)-nanmean(tl),pl(i)-nanmean(pl),'.','Markersize',15,'Color', cmapVar(i,1,length(tl),jet)); axis square; hold on
% 
%  end
% 
% % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %    
%     end
% end
%%

clear all
% load('ACCAnalyzed_AllAnimals_010820_noDLS.mat')
load('ACC_AllAnimals_021520_a.mat')
savePDF=1;
    psfilename = 'C:\analysisPS_B.ps';
    
 % Figure 2D: example of convergence during approach
   
%  appRange=60:140;
appRange=97-45:151+30;
 appTimes=97-45:106;
% appRange=474-100:502+100

tR=Rtheta{useData(67)}(appRange,:); tR=tR-nanmean(tR);
tL=Ltheta{useData(67)}(appRange,:);tL=tL-nanmean(tL);
tilt=accelChannels{useData(67)}(appRange,2);  tilt=tilt-nanmean(tilt); tilt=medfilt1(tilt,8);
roll=accelChannels{useData(67)}(appRange,1);  roll=roll-nanmean(roll); roll=medfilt1(roll,5);
pR=Rphi{useData(67)}(appRange,:); pR=pR-nanmean(pR);
pL=Lphi{useData(67)}(appRange,:);pL=pL-nanmean(pL);
pVg=pR-pL; pVg=pVg-nanmean(pVg);
hT=thetaHead{useData(67)}(1,appRange); hT=hT-nanmean(hT);
gaze=(.5*(tR+tL))+hT'; %gaze=gaze-nanmean(gaze);

% figure
% subplot(2,1,1)
%  appTimes=find(appEpoch{67}(appRange)==1)
% % appTimes=97:151;
% plot(tR); hold on;
% plot(tL); hold on;
% plot([appTimes(1) appTimes(1)],[-30 30],'g'); plot([appTimes(end) appTimes(end)],[-30 30],'r');
% ylim([-60 60]);
% ylabel('degrees');xlabel('time'); 
% legend('right eye','left eye');
% 
% subplot(2,1,2)
% plot(roll);hold on;
% plot(tilt);
% plot([appTimes(1) appTimes(1)],[-30 30],'g'); plot([appTimes(end) appTimes(end)],[-30 30],'r');
% hold on;
% % plot(tR(appRange)-tL(appRange));
% ylim([-60 60]); title('acc tilt');
% legend('roll','tilt');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% 
% figure
% subplot(2,1,1)
% % appTimes=find(appEpoch{67}(appRange)==1)
% % appTimes=97:151;
% plot(tR-tL); hold on;
% plot(tilt); hold on;
% plot([appTimes(1) appTimes(1)],[-30 30],'g'); plot([appTimes(end) appTimes(end)],[-30 30],'r');
% ylim([-60 60]);
% ylabel('degrees');xlabel('time'); 
% legend('vergence','tilt');
% 
% subplot(2,1,2)
% plot(pVg);hold on;
% plot(roll);
% plot([appTimes(1) appTimes(1)],[-30 30],'g'); plot([appTimes(end) appTimes(end)],[-30 30],'r');
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
% plot([appTimes(1) appTimes(1)],[-30 30],'g'); plot([appTimes(end) appTimes(end)],[-30 30],'r');
% ylabel('degrees');xlabel('time'); 
% legend('gaze','head');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

% Figure 2D: example of convergence during approach
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
% plot([appTimes(1) appTimes(1)],[-30 30],'g'); plot([appTimes(end) appTimes(end)],[-30 30],'r');
% ylim([-60 60]);
% ylabel('degrees');xlabel('time'); 
% legend('right eye','left eye');
% 
% subplot(2,1,2)
% plot(roll(appRange));hold on;
% plot(tilt(appRange));
% plot([appTimes(1) appTimes(1)],[-30 30],'g'); plot([appTimes(end) appTimes(end)],[-30 30],'r');
% hold on;
% % plot(tR(appRange)-tL(appRange));
% ylim([-60 60]); title('acc tilt');
% legend('roll','tilt');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

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
% plot([appTimes(1) appTimes(1)],[-30 30],'g'); plot([appTimes(end) appTimes(end)],[-30 30],'r');
% % ylim([-40 40]);
% 
% subplot(2,1,2)
% plot(tilt(appRange));hold on;
% plot([appTimes(1) appTimes(1)],[-30 30],'g'); plot([appTimes(end) appTimes(end)],[-30 30],'r');
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
            vgHist = (hist(vg(useN),nBins))/sum(useN&~isnan(vg)');
        else
            vgHist=NaN;
        end
        vgDiff(vid,:,c+1)=vgHist;
    end
end

vergCounts(:,1)=nansum(vgDiff(:,:,1),1)./(length(useData))
vergCounts(:,2)=nansum(vgDiff(:,:,2),1)./(length(useData))
vergErr=nanstd(vgDiff)./(sqrt(length(vgDiff))); vergErr=squeeze(vergErr);

figure; shadedErrorBar(nBins,vergCounts(:,1),vergErr(:,1),'b',1); hold on
shadedErrorBar(nBins,vergCounts(:,2),vergErr(:,2),'g',1); hold on

title('Figure 2F: vergence (R - L eye)'); axis square
xlabel('vergence (deg)'); ylabel('proportion of time');
ylim([0 .22])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
[h p]=kstest2(vergCounts(:,1),vergCounts(:,2))

% Figure 2G: correlation of eye thetas & eye phis 
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
corrRAAll=[]
for vid=1:length(useData)
    nframe = min(length(dRtheta{useData(vid)}),length(dLtheta{useData(vid)}));
    dtR=dRtheta{useData(vid)}(1:nframe); dtL=dLtheta{useData(vid)}(1:nframe);
    nonapp=appEpoch{vid}==0;
    use = (nonapp==1)';%& ~isnan(dtR(1:nframe));
    if sum(use)>3
        [corrR lagsR]= nanxcorr(dtR(use),dtL(use),30,'coeff');
        uselagsR=(lagsR>=-30& lagsR<=30);

        clear nframe
        nframe = min(length(dRphi{useData(vid)}),length(dLphi{useData(vid)}));
        dpR=dRphi{useData(vid)}(1:nframe); dpL=dLphi{useData(vid)}(1:nframe);
        [corrL lagsL]= nanxcorr(dpR(use),dpL(use),30,'coeff');
        uselagsL=(lagsL>=-30 & lagsL<=30);

    else
    end
    use=appEpoch{vid}==1;
    if sum(use)>3 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dtR(use),dtL(use),30,'coeff');
        uselagsRA=(lagsRA>=-30& lagsRA<=30);

        nframe = min(length(dRphi{useData(vid)}),length(dLphi{useData(vid)}));
        dpR=dRphi{useData(vid)}(1:nframe); dpL=dLphi{useData(vid)}(1:nframe);
        [corrLA lagsLA]= nanxcorr(dpR(use),dpL(use),30,'coeff');
        uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    else
    end

    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
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
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.7 .7]);
xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dTheta non-app','dPhi non-app','dTh Approach','dPhi app'}); title('Figure 2G: mean between eye corr, dtheta & dphi');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%% FIGURE 3: during approaches, eyes are more centered (smaller vergence), due to less head movements in pitch and roll

% Figure 3A: overlaid traces of pitch & eye vergence

% Figure 3B: scatter & corr of pitch (tilt) and vergence
close all

for c=0:1
    clear use
use = find(appAll==c);

figure(1)
plot(tiltAll(use(1:20:end)),vergDlc(use(1:20:end)),'.'); axis equal; hold on; lsline
xlabel('acc tilt'); ylabel('eye vergence');
ylim([-90 90]); xlim([-90 90]);
R=corrcoef(tiltAll(use), vergDlc(use),'Rows','pairwise')
text(-80,(80-c*5), ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
title('Figure 3B: acc tilt and dlc eye vergence');
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end


figure(2)
[corr lags]=(nanxcorr(tiltAll(use),vergDlc(use),30,'coeff'));
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
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end
end

clear tiltHist
for c=0:1
    for vid=1:length(useData)
     clear useN tilt tiltHist
        tilt=accelChannels{useData(vid)}(:,2);
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

%this doesn't add to 1 for approaches 
tiltCounts(:,1)=nansum(tiltFull(:,:,1),1)./(length(useData))
tiltCounts(:,2)=nansum(tiltFull(:,:,2),1)./((length(useData))-sum(isnan(tiltFull(:,1,2))))
tiltErr=nanstd(tiltFull)./(sqrt(length(tiltFull))); tiltErr=squeeze(tiltErr);

figure; shadedErrorBar(nBins,tiltCounts(:,1),tiltErr(:,1),'b',1); hold on
shadedErrorBar(nBins,tiltCounts(:,2),tiltErr(:,2),'g',1); hold on

title('Figure 3B: pitch'); axis square
xlabel('pitch (deg)'); ylabel('proportion of time');
ylim([0 .16])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
[h p]=kstest2(tiltCounts(:,1),tiltCounts(:,2))


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
R=corrcoef(rollAll(use), dlcPhi(use),'Rows','pairwise')
text(-80,(80-c*5), ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
title('Figure 3F: acc roll & diff of eye phi)');
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

end

figure(5)
[corr lags]=(nanxcorr(rollAll(use),dlcPhi(use),30,'coeff'));
plot(lags,corr);axis square; hold on
title('Figure 3G: corr of acc roll & diff between eye phi'); ylabel('corr coeff')
ylim([-1 1]);
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end

figure(6)
clear nBins h
nBins=-90:5:90
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
plot(nBins, h/length(gyro3All(use))); hold on; title('head yaw from gyro')
title('Figure 4A: head angle is not diff between non-app and app')
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end
end
%%% stats here

% Figure 4B: hist of mean eye yaw, app and non-app
d_mnEyeAll=[];mnPhiAll=[];allDLChead=[];mnEyeAll=[]; mouseSpAll=[];
for vid=1:length(useData)
    nframe=min(length(dRtheta{useData(vid)}),length(dLtheta{useData(vid)}));
    nframe=min(nframe, length(d_Theta{useData(vid)}));
mnEye =.5*(Rtheta{useData(vid)}(1:nframe)+Ltheta{useData(vid)}(1:nframe));
mnEye=mnEye-nanmean(mnEye);
mnEyeAll=[mnEyeAll mnEye'];

mnEyeD =.5*(dRtheta{useData(vid)}(1:nframe)+dLtheta{useData(vid)}(1:nframe));
mnEyeD=mnEyeD-nanmean(mnEyeD);
d_mnEyeAll=[d_mnEyeAll mnEyeD'];

mnPhi =.5*(dRphi{useData(vid)}(1:nframe)+dLphi{useData(vid)}(1:nframe));
mnPhi=mnPhi-nanmean(mnPhi);
mnPhiAll=[mnPhiAll mnPhi'];
dHead=d_Theta{useData(vid)}(1:nframe); dHead=dHead-nanmean(dHead);
allDLChead = [allDLChead dHead'];

speed =mouseSp{useData(vid)}(1:nframe);
speed =speed-nanmean(speed);
mouseSpAll=[mouseSpAll speed];

end

figure
for c=0:1
    hold on
    clear nBins h use hp
    use = find(appAll==c);
    nBins= -60:2:60
    h=hist(mnEyeAll(use),nBins);
    subplot(1,2,1)
    plot(nBins,h/length(mnEyeAll(use))); hold on; title('mn Eye Yaw'); axis square
    plot([-20,-20],[.15,0],'k--');
    plot([20,20],[.15,0],'k--');
    
    nBins=-30:2:30
    hp=hist(d_mnEyeAll(use),nBins);
    subplot(1,2,2)
    plot(nBins, hp/length(d_mnEyeAll(use))); axis square; hold on; title('mn d Eye Theta')
    
    if c==1
        suptitle('Figure 4B: mn Eye yaw falls within binocular zone');
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
end


% Figure 4C: ex trace when dTh is 0, eyes are also 0, then hist

figure
for c=0:1
    hold on
    clear nBins h use hp
    use = (appAll==c);
    %still=(allDLChead<.5 & allDLChead>-.5);
    still=(gyro3All<.5 & gyro3All>-.5);
    nBins= -20:1:20;
    h=hist(d_mnEyeAll(use&still),nBins);
    subplot(1,2,1)
    plot(nBins,h/sum(use&still)); hold on; title('mn Eye Yaw'); axis square
    plot([-3,-3],[.6,0],'k--');
    plot([3,3],[.6,0],'k--');
    prop3(c+1,1) =1- sum(h(nBins<-3 | nBins>3))/sum(h);
   hp=hist(mnPhiAll(use& still),nBins);
    subplot(1,2,2)
    plot(nBins, hp/sum(use&still)); axis square; hold on; title('mn Eye Phi')
    prop3(c+1,2)=1- sum(hp(nBins<-3 | nBins>3))/sum(hp);

    plot([-3,-3],[.6,0],'k--');
    plot([3,3],[.6,0],'k--');
    if c==1
        suptitle('Figure 4C: when head angle is not changing, eyes are also still');
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
end
% 
figure
clear nBins h use hp
stationary = mouseSpAll<5;
nBins= -40:1:40;
hs=hist(d_mnEyeAll(stationary),nBins);
plot(nBins,hs/sum(stationary)); hold on; title('mn Eye Yaw, stationary'); axis square

hma=hist(d_mnEyeAll(stationary==0),nBins);
plot(nBins,hma/sum(stationary==0)); hold on; title('mn Eye Yaw, moving approach'); axis square
% h=hist(d_mnEyeAll(appAll==1),nBins);
% plot(nBins,h/sum(appAll==1)); hold on; title('mn Eye Yaw, approach'); axis square
legend('stationary','moving');
plot([-3,-3],[.3,0],'k--');
plot([3,3],[.3,0],'k--');
ylim([0 .32]); xlim([-22 22]);
title('Figure 4C: when mouse is stationary, so are eyes');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end



% Figure 4D: scatter of change in head yaw and change in eye yaw, shows
% mostly congruent but not all
figure
for c=1:2
    clear use 
    use = find(appAll==c-1);
    samp(c,:)=randsample(use,7375);
    plot(gyro3All(samp(c,:)), d_mnEyeAll(samp(c,:)),'.')
    axis([-35 35 -35 35]); axis square
    hold on
end
    xlabel('gyro yaw'); ylabel('d eye yaw');

title('Figure 4D: head & eye thetas, not all compensatory')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

% Figure 4E: corr of change in head and eye yaw, congruence
% (non-compensation) at short time scales


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
    [corrYaw lagsYaw]= nanxcorr(g3(use),mnEye(use),30,'coeff');
     uselags=(lagsYaw>=-30& lags<=30);
 end
    if sum(uselags)==61 
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

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); 
 ylim([-.5 .3]);
xlim([21 41]); %15 and 46 ==500 ms
xlim([8.5 53.5])
xlabel('time'); ylabel('correlation coeff');
axis square
title('Figure 4E: eye movements are mostly compensatory');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%% FIGURE 5: Non-compensatory eye movements
gyro3All = gyro3All - nanmedian(gyro3All);


%%% cluster on dHead vs dEye
%%%old way
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
gscatter(gyro3All(appAll==1),d_mnEyeAll(appAll==1),idx); axis equal
title('cluster on dHead vs dEye 3 clust')
xlim([-25 25]); ylim([-25 25]); %hold on; plot([-25 25], [25 -25],'r')
clear mvmts
mvmts=[gyro3All(appAll==1); d_mnEyeAll(appAll==1)] %; mnEyeAll(appAll==1)];
%mvmts = [gyro3All(appAll==1) + d_mnEyeAll(appAll==1)];
gm = fitgmdist(mvmts',3,'Replicates',10);
idx = cluster(gm,mvmts');

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
figure
gscatter(gyro3All(appAll==1),d_mnEyeAll(appAll==1),clust1); axis equal
title('clutering based on gaussian fit to dGz')


if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


figure
gscatter(gyro3All(appAll==1)-nanmedian(gyro3All),d_mnEyeAll(appAll==1),idx); axis equal; hold on
title('all approach pts only')
xlim([-25 25]); ylim([-25 25]); hold on; plot([-25 25],[25 -25],'g')

%% use ~clust1
% fullData=[gyro3All;d_mnEyeAll]%;mnEyeAll];
dGzFull= gyro3All + d_mnEyeAll;
fullData=[gyro3All; d_mnEyeAll];
idxAll=cluster(gm, fullData');
[sacc, clust]= max([nanmean(abs(dGzFull(idxAll==1))),nanmean(abs(dGzFull(idxAll==2))), nanmean(abs(dGzFull(idxAll==3)))])

% figure;
% [corr lags]=nanxcorr(gyro3All(idxAll==clust),d_mnEyeAll(idxAll==clust),30,'coeff')
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
maxlag = 30;
thbins = -60:5:60;
skip = 1; %%% only shows figures at this interval
nthresh  = 60;
headAll=[];
eyeAll=[];
% figure
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
    
    %%% get accelerometers

    if exist('accelData','var')
        tilt = accelChannels{vid}(:,2);
        roll = accelChannels{vid}(:,1);
        acc_dth = accelChannels{vid}(:,6);
        acc_dth=acc_dth-nanmedian(acc_dth);
            acc_hth = nancumsum([hth(1) acc_dth(1:end-1)'],[],2);
    else
        display('no acc')
    end

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
    gaze = hth + mnEyeTh';
    
    %%% find the longest approach, to use as example image
    appPts = find(app);
    newApp = [1 find(diff(appPts)>1)+1];
    endApp = [find(diff(appPts)>1)-1 length(appPts)];
    
    dur = endApp - newApp;
    [mx longest] = max(dur);
    mainApp = appPts(newApp(longest)  : endApp(longest)); %%% approach time only
    %%% add 2 secs on either side
    try
        appStart = max(appPts(newApp(longest))-60,1); appOffset = appPts(newApp(longest))-appStart;
        appEnd = min(appPts(endApp(longest))+60,length(app)); endOffset = appPts(endApp(longest))-appStart;
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
   hthApp = acc_hth(appRange) - nanmedian(acc_hth(mainApp));
    hthApp = mod(hthApp + 180,360)-180;
     
     
    gzApp = hthApp +0.5*( rth(appRange) +lth(appRange))';
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
        g3=accelChannels{vid}(:,6); g3=g3-nanmedian(g3);
        
        %dGaze=diff(gaze);
        dGaze = g3(1:end-1) + 0.5*(drth+dlth);
        
%         figure
%         plot(g3(1:end-1), 0.5*(drth+dlth),'.'); axis([-20 20 -20 20])

       % nBins=-20:0.25:20;
       
       %%% assign clusters based on previous GMM
       newData=[g3(1:end-1)' ; 0.5*(drth+dlth)'];
        idx2=cluster(gm,newData');



        %%% cluster based on gaussian fit for compensatory points
        %%% probability distribution from fit is define by pfit (from above)
        %%% see where new data lands on this (using interp1) then use
        %%% random probability to assign each one
           prob = interp1(nBins,pfit,dGaze);
            clust1 = rand(size(prob))<prob;
            clust1(isnan(dGaze))=1;
         idx2=~clust1;
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


figure('Renderer', 'painters', 'Position', [100 -100 600 900])
        subplot(4,1,1)
        plot(0.5*( rth(appRange) +lth(appRange))','k','LineWidth',2); hold on; 
        plot(rth(appRange)','r'); hold on; plot(lth(appRange)','b');
        plot(find(idx2(appRange)==clust),0.5*(rth(appRange(idx2(appRange)==clust)) +lth(appRange(idx2(appRange)==clust)))','og')
        ylim([-30 30]);  xlim([0 max(length(appRange),1)]);
%         xlim([40 110]);
        plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        legend('mean','right','left');
         title(sprintf('vid %d',vid));


        subplot(4,1,2)
        plot(0.5*( drth(appRange) +dlth(appRange))','k','LineWidth',2); hold on; 
        plot(drth(appRange)','r'); hold on; plot(dlth(appRange)','b');
        plot(find(idx2(appRange)==clust),0.5*(drth(appRange(idx2(appRange)==clust)) +dlth(appRange(idx2(appRange)==clust)))','og')
        ylim([-20 20]);  xlim([0 max(length(appRange),1)]);
        plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        title('d eye theta')
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

subplot(4,1,4)
hold on
plot(dhthApp','Color',[0 0.75 0],'LineWidth',2); hold on
plot(dhthApp' +0.5*( drth(appRange) +dlth(appRange))','k','LineWidth',.75);
        plot(find(idx2(appRange)==clust),dhthApp(idx2(appRange)==clust),'bo')
ylim([-50 50])
title('d head, dgaze')

if appOffset>0
    plot([appOffset appOffset],[-60 60],'g');end
plot([endOffset endOffset],[-60 60],'r');
if ~isempty(min(hthApp))
    %ylim([min(hthApp)-20 max(hthApp)+20]);
    ylim([-50 50])
else  ylim([-50 50]);end
        

gzApp = hthApp +0.5*( rth(appRange) +lth(appRange))';
% subplot(6,1,6);
subplot(4,1,3)
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
% %               [corr lags]=nanxcorr(head, eyeTrace(useT),30,'coeff');
% %               corrVid(i,:)=corr; lagsVid(i,:)=lags;
% %           end
%     headAll=[headAll head'];
%     eyeAll=[eyeAll eyeTrace(useT)'];
    end


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
    filen=sprintf('%s','PaperFigs_021620_dGazeFit ','.pdf')
%     filen=sprintf('%s','PaperFigs_011519_c','.pdf')
    pdfilename=fullfile(pSname,filen);
    dos(['ps2pdf ' psfilename ' ' pdfilename]);
    delete(psfilename);
else
    pSname='T:\PreyCaptureAnalysis\Data\';
end









