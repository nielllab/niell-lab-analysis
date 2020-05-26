clear all; close all;
dbstop if error
% load('J462a_AlACCSessions_120619_a.mat');
% load('J475c_DEINTERLACED_051920_a.mat');
% load('J475c_DEINTERLACED_052120_shiftTest.mat')
load('J475c_DEINTERLACED_052320_halfShift.mat')
set(groot,'defaultFigureVisible','on') %disable figure plotting
deInter=1;

if deInter
frRate=60;
else
frRate=30; 
end

savePDF=1;
if savePDF
    psfilename = 'C:\analysisPS3.ps';
    if exist(psfilename,'file')==2;delete(psfilename);end
end

mouse_xy=[];cricket_xy=[];EllipseParamsR=[];EllipseParamsL=[]; radR=[]; az=[];

for i=1:length(Data)
    
    animal(i,:)=Data(i).ani;
    sessionN(i) = Data(i).sessionnum;
    expdate(i)=Data(i).date;
    clip(i)=Data(i).clipnum;
    
    mouse_xyRaw{i,1,1}=Data(i).mouse_xyRaw;
    mouseVRaw{i,:}= Data(i).mouseVRaw; %in pix/frame
    mouseVRaw{i,:}=((Data(i).mouseVRaw)/27)*frRate; %cm/sec
    cricket_xyRaw{i,1}=Data(i).cricketxyRaw;
    cricketVRaw{i,:}= Data(i).cricketVRaw % pix/frame
    cricketVRaw{i,:} = ((Data(i).cricketVRaw)/27)*frRate; %now cm/sec
    thetaRaw{i,:}= Data(i).thetaRaw;
    dThetaRaw{i,:}= diff(Data(i).thetaRaw);
    rangeRaw{i,:}= (Data(i).rangeRaw)/27; %convert pixels to cm
    azRaw = mod(Data(i).azRaw,2*pi); azRaw(azRaw>pi) = azRaw(azRaw>pi)-2*pi;
    azTRaw{i,:,:}= azRaw;
    azDegRaw{i,:,:}=rad2deg(azTRaw{i,:,:});%az in deg
    %     crickHRaw{i,1,1} = Data(i).cricketHead;
    %     crickH{i,1,1} = Data(i).crickH;
    
    mouse_xy{i,1,1}=Data(i).mouse_xy;
    % mouseV{i,:}= Data(i).mouseV; %in pix/frame
    mouseV{i,:}=((Data(i).mouseV)/27)*frRate; %cm/sec
    cricket_xy{i,1}=Data(i).cricketxy;
    cricketAz{i,:}=Data(i).cricketTheta;
    %    cricketV{i,:}= Data(i).cricketV % pix/frame
    cricketV{i,:} = ((Data(i).cricketV)/27)*frRate; %now cm/sec
    theta{i,:}= rad2deg(Data(i).theta);
    dTheta{i,:}=rad2deg(Data(i).dth);
%     dThetaFilt{i,:}=medfilt1(dTheta{i,:},3);
    range{i,:}= (Data(i).range)/27; %convert pixels to cm
    az = mod(Data(i).az,2*pi); az(az>pi) = az(az>pi)-2*pi;
    azT{i,:,:}= az;
    azDeg{i,:,:}=rad2deg(azT{i,:,:});%az in deg
    slip{i,1} = Data(i).difTR;
    slip{i,2} = Data(i).difTL;
    slip{i,3} = Data(i).difRL;
    goodTheta(i,:) = Data(i).ThetaFract; %from un-interpolated theta
    TSused{i,:}=Data(i).usedTS;
    thMissing(i)= sum(~isnan(Data(i).theta))/(length(Data(i).theta)); %proportion of non-nans
    
    accelData{i,:} = Data(i).accShift;
  %  accelCorr(i) =Data(i).accXcorrMax;
%     accelDrift(i) =Data(i).accXcorrLag;
    accelDataRaw{i,:}= Data(i).rawAccShift;
    
    EllipseParamsR{i,1,1} = Data(i).EllipseParamsR;
    EllipseParamsL{i,1,1} = Data(i).EllipseParamsL;
    thetaR{i,:} = Data(i).Rtheta;
    thetaL{i,:} =Data(i).Ltheta;
    phiR{i,:} =Data(i).Rphi;
    phiL{i,:} =Data(i).Lphi;
    rMissing(i)=sum(~isnan(Data(i).Rtheta))/(length(Data(i).Rtheta)); %proportion of non-nans
    lMissing(i)=sum(~isnan(Data(i).Ltheta))/(length(Data(i).Ltheta));
    
    dthetaR{i,:} =(Data(i).dxRTheta);
    dthetaL{i,:} =(Data(i).dxLTheta);
    dphiR{i,:} =(Data(i).dxRPhi);
    dphiL{i,:} =(Data(i).dxLPhi);
    
    dXRcent{i,:} =Data(i).dxR;
    dXLcent{i,:} =Data(i).dxL;
    %     dYRcent{i,:} =Data(i).dyR;
    %     dYLcent{i,:} =Data(i).dyL;
    XRcent{i,:} =Data(i).XRcent;
    YRcent{i,:} =Data(i).YRcent;
    XLcent{i,1} = Data(i).XLcent;
    YLcent{i,:} =Data(i).YLcent;
    RRad{i,:}=Data(i).RRad;
    LRad{i,:}=Data(i).LRad;
%     
%     if deInter
%     goodR{:,i}=Data(i).goodReye;
%     else
%     goodR{i,:}=Data(i).goodReye; %1=all 8pts above likelihood .95
%     end
    Rngood{i,:}=Data(i).ngoodR; % num good DLC pts
    Rcc(i)=Data(i).RcalR;
    Rslope(i)=Data(i).RcalM;
    Rscale(i)=Data(i).scaleR;
    
%     goodL{i,:}=Data(i).goodLeye;
    Lngood{i,:}=Data(i).ngoodL;
    Lcc(i)=Data(i).LcalR;
    Lslope(i)=Data(i).LcalM;
    Lscale(i)=Data(i).scaleL;
    
    
    %     longR{i,1}=EllipseParamsR{i,1}(:,3);longL{i,1}=EllipseParamsL{i,1}(:,3);
    %     shortR{i,1}=EllipseParamsR{i,1}(:,4);shortL{i,1}=EllipseParamsL{i,1}(:,4);
    %
    %     radR{i,1}= (longR{i,1}+shortR{i,1})./length(longR{i,1});
    %     radL{i,1}= (longL{i,1}+shortL{i,1})./length(longL{i,1});
    %     pupilRvel{i,1}=diff(radR{i,1}); pupilLvel{i,1}=diff(radL{i,1})
    
    % tsData{i,1}=~isempty(Data(i).TopTs);
    
end
%
% tsData= cell2mat(tsData)
% delayFull=cell2mat(slip);
goodR= Rcc>.3 & rMissing>.75;
goodL=Lcc>.3 & lMissing>.75;
% useL = (delayFull(:,2)<=3 & delayFull(:,2)>=-3);
% useR = (delayFull(:,1)<=3 & delayFull(:,1)>=-3);
% useE = (delayFull(:,3)<=3 & delayFull(:,3)>=-3);

useTime = goodTheta>=.9 & goodR' & goodL' & thMissing'>.90; %tsData==1; %|(useL & useR)
useFilt=find(useTime); %useFilt=useFilt(1:4,6:end);

rownum=10 ; colnum=8;
% rownum=round(sqrt(length(useFilt)+4))
% colnum=round(sqrt(length(useFilt)));
%%
% figure
% subplot(3,2,1)
% hist(Rcc); title('R corrcoef'); ylim([0 frRate]); xlim([0 1]);
% subplot(3,2,3)
% hist(Rslope);title(' R Cal slope');ylim([0 frRate]); xlim([0 1]);
% subplot(3,2,2)
% hist(Lcc); title('L corrcoef');ylim([0 frRate]); xlim([0 1]);
% subplot(3,2,4)
% hist(Lslope);title('L Cal slope');ylim([0 frRate]); xlim([0 1]);
% subplot(3,2,5)
% hist(Rscale);title('R Cal scale');ylim([0 frRate]); xlim([20 80]);
% subplot(3,2,6)
% hist(Lscale);title('L Cal scale');ylim([0 frRate]); xlim([20 80]);
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
figure;
subplot(2,2,1)
plot(Rcc,Rslope,'o'); axis square; xlim([0 1]); ylim([0 1]); title('R');
xlabel('R corrcoef');ylabel('R Slope');
hold on; plot([0,1],[0,1]);
subplot(2,2,2)
plot(Lcc,Lslope,'o');axis square; xlim([0 1]); hold on;ylim([0 1]);plot([0,1],[0,1]);
title('L');
xlabel('L corrcoef');ylabel('L Slope');
subplot(2,2,3)
plot(Rcc,Rscale,'o'); axis square; xlim([0 1]); ylim([0 100]); title('R');
xlabel('R corrcoef');ylabel('R Scale');
hold on; plot([0,1],[0,100]);
subplot(2,2,4)
plot(Lcc,Lscale,'o');axis square; xlim([0 1]); hold on;ylim([0 100]);plot([0,1],[0,100]);
title('L');
xlabel('L corrcoef');ylabel('L Scale');

if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%
allTilt=[];allRoll=[];allYaw=[]; dlcVerg=[];dlcDverg=[];allGyro1=[];allGyro2=[];allGyro3=[];
dlcDverg=[]; dlcDphi=[]; dlcPhi=[];dlcDhth=[];dlcHth=[];
% figure
for vid = 1:length(useFilt)
    figure('units','normalized','outerposition',[0 0 1 1])
    roll = (accelData{useFilt(vid)}(:,1)); roll=roll-nanmean(roll);
    rollFilt = medfilt1(roll,2);
    tilt = (accelData{useFilt(vid)}(:,2)); tilt=tilt-nanmean(tilt);
    tiltFilt = medfilt1(tilt,8);
    yaw = (accelData{useFilt(vid)}(:,3)); yaw=yaw-nanmean(yaw);
    yawFilt = medfilt1(yaw,10);
    verg=(thetaR{useFilt(vid)}-thetaL{useFilt(vid)});
    dverg=(dthetaR{useFilt(vid)}-dthetaL{useFilt(vid)});
    gyro1=(accelData{useFilt(vid)}(:,4));
    gyro2=accelData{useFilt(vid)}(:,5);
    gyro3=(accelData{useFilt(vid)}(:,6));
    phi=(phiR{useFilt(vid)}-phiL{useFilt(vid)}); phi=phi-nanmean(phi);
    dphi=(dphiR{useFilt(vid)}-dphiL{useFilt(vid)}); dphi=dphi-nanmean(dphi);
    badE=dphi>15| dphi<-15 ;
    dphi(badE)=nan;
    hth=theta{useFilt(vid)}; hth=hth-nanmean(hth);
    dhth=dTheta{useFilt(vid)}; dhth=dhth-nanmean(dhth);
    bad=dhth>20 | dhth<-20;
    dhth(bad)=nan;
%     
%     subplot(4,3,1)
%     plot(rollFilt); hold on; axis square;
%     plot(phi);
%     title('acc 1 - roll, phi')
%     subplot(4,3,2)
%     plot(tiltFilt); hold on; axis square;plot(verg);
%     title('acc 2 - tilt, vg')
%     subplot(4,3,3)
%     plot(yawFilt); hold on; axis square; plot(verg);
%     title('acc 3 - yaw, vg')
%     subplot(4,3,4)
%     plot(gyro1); hold on; axis square; plot(verg);
%     title('gyro 1, vg')
%     subplot(4,3,5)
%     plot(gyro2); hold on; axis square; plot(dphi);
%     title('gyro 2, dPhi')
%     subplot(4,3,6)
%     plot(gyro3); hold on; axis square
%     plot(dhth); title('gyro 3, dHead Th')
%     subplot(4,3,7);
%     plot(-frRate:frRate,nanxcorr(rollFilt,phi,frRate,'coeff'))
%     title('roll, phi diff');ylim([-1 1]);
%     subplot(4,3,8);
%     plot(-frRate:frRate,nanxcorr(tiltFilt,verg,frRate,'coeff'))
%     title('tilt, vergence');ylim([-1 1]);
%     subplot(4,3,9);
%     plot(-frRate:frRate,nanxcorr(yawFilt,verg,frRate,'coeff'))
%     title('yaw, vergence');ylim([-1 1]);
%     subplot(4,3,10)
%     plot(-frRate:frRate,nanxcorr(gyro1(1:end-1),dverg,frRate,'coeff'))
%     title('gyro 1, dVg');ylim([-1 1]);
%     subplot(4,3,11)
%     plot(-frRate:frRate,nanxcorr(gyro2(1:end-1),dphi,frRate,'coeff'))
%     title('gyro 2, dPhi eyes');ylim([-1 1]);
%     subplot(4,3,12)
%     plot(-frRate:frRate,nanxcorr(gyro3,dhth,frRate,'coeff'))
%     title('gyro 3, dHead th');ylim([-1 1]); 
    
    allTilt=[allTilt tiltFilt(1:end-1)'];
    allRoll=[allRoll rollFilt(1:end-1)'];
    allYaw=[allYaw yawFilt(1:end-1)'];
    dlcVerg=[dlcVerg verg(1:end-1)];
    dlcDverg=[dlcDverg dverg];
    allGyro1=[allGyro1 gyro1(1:end-1)'];
    allGyro2=[allGyro2 gyro2(1:end-1)'];
    allGyro3=[allGyro3 gyro3(1:end-1)'];
    dlcPhi=[dlcPhi phi(1:end-1)];
    dlcDphi=[dlcDphi dphi];
    dlcHth=[dlcHth hth(1:end-1)];
    dlcDhth=[dlcDhth dhth(1:end-1)'];
    if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end
%%
% figure
% plot(allTilt,dlcVerg,'.'); axis equal; hold on; lsline
% mdl = fitlm(allTilt,dlcVerg)
% axis([-90 90 -90 90]);
% xlabel('acc pitch'); ylabel('vergence')
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% figure
% plot(allRoll,dlcPhi,'.'); axis equal; hold on; lsline
% axis([-90 90 -90 90]);
% xlabel('acc roll'); ylabel('eye phi')
% mdl = fitlm(allRoll,dlcPhi)
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
figure
subplot(2,3,1)
plot(allRoll,dlcPhi,'.'); axis equal; hold on; lsline
xlabel('acc roll'); ylabel('eye phi');
ylim([-90 90]); xlim([-90 90]);
R=corrcoef(allRoll, dlcPhi,'Rows','pairwise')
text(-80,80, ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
title('acc ch 1');

subplot(2,3,2)
plot(allTilt,dlcVerg,'.'); axis equal; hold on; lsline
xlabel('acc tilt'); ylabel('eye vergence');
ylim([-90 90]); xlim([-90 90]);
R=corrcoef(allTilt, dlcVerg,'Rows','pairwise')
text(-80,80, ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
title('acc ch 2');

subplot(2,3,3)
plot(allYaw,dlcVerg,'.'); axis equal; hold on; lsline
xlabel('acc yaw'); ylabel('eye verg');
ylim([-90 90]); xlim([-90 90]);
R=corrcoef(allYaw, dlcVerg,'Rows','pairwise')
text(-80,80, ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
title('acc ch 3');

subplot(2,3,4)
plot(allGyro1,dlcDverg,'.');axis equal; hold on; lsline %dlcVerg
xlabel('gyro 1'); ylabel('eye verg');
ylim([-frRate frRate]); xlim([-frRate frRate]);
R=corrcoef(allGyro1, dlcDverg,'Rows','pairwise')
text(-25,25, ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
title('gyro ch 1')

subplot(2,3,5)
plot(allGyro2,dlcDphi,'.');axis equal; hold on; lsline

xlabel('gyro 2'); ylabel('d phi');
ylim([-50 50]); xlim([-50 50]);
R=corrcoef(allGyro2, dlcDphi,'Rows','pairwise')
text(-40,40, ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
title('gyro ch 2')

test=(dlcDhth>10|dlcDhth<-10)&allGyro3>-5& allGyro3<5;

subplot(2,3,6)
plot(allGyro3(test==0),dlcDhth(test==0),'.');axis equal; hold on; lsline
plot(allGyro3(test==1),dlcDhth(test==1),'r.');axis equal; hold on; lsline
xlabel('gyro 3'); ylabel('d head theta');
ylim([-90 90]); xlim([-90 90]);
R=corrcoef(allGyro3(test==0), dlcDhth(test==0),'Rows','pairwise')
text(-80,80, ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
title('gyro ch 3')


if savePDF
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end
close(gcf)
% % test=allGyro3>-5& allGyro3<5  & (dlcDhth>15|dlcDhth<-15);
% plot(allGyro3(test),dlcDhth(test),'.')
% RS=corrcoef(allGyro3(test==0), dlcDhth(test==0),'Rows','pairwise')
% text(-80,70, ['corrcoef = ' num2str(RS(1,2),'%.2f')],'FontSize',10)

% dlcDhth(test)=NaN;
%%
test=(dlcDhth>10|dlcDhth<-10)&allGyro3>-5& allGyro3<5;

figure
for vid = 1:length(useFilt)
    subplot(5,6,vid)
    plot(dTheta{useFilt(vid)}); %blue =dlc
    hold on;
    plot(accelData{useFilt(vid)}(:,6));
    exptest=accelData{useFilt(vid)}(:,6)>-5& accelData{useFilt(vid)}(:,6)<5  & ((dTheta{useFilt(vid)})>15 |(dTheta{useFilt(vid)})<-15);
    plot(find(exptest),dTheta{useFilt(vid)}(exptest==1),'og')
end

if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%
% figure
% subplot(2,3,1)
% [corr lags]=(nanxcorr(allRoll,dlcPhi,frRate,'coeff'));
% plot(lags,corr);axis square
% title('acc roll, phi');
% ylim([-1 1]);

% subplot(2,3,2)
% [corr lags]=(nanxcorr(allTilt,dlcVerg,frRate,'coeff'));
% plot(lags,corr)
% axis square;ylim([-1 1]);
% title('acc tilt, vergence');
% 
% subplot(2,3,3)
% [corr lags]=(nanxcorr(allYaw,dlcVerg,frRate,'coeff'));
% plot(lags, corr)
% axis square;ylim([-1 1]);
% title('acc yaw, vergence');
% 
% subplot(2,3,4)
% [corr lags]=(nanxcorr(allGyro1,dlcDverg,frRate,'coeff'));
% plot(lags, corr)
% axis square;ylim([-1 1]);
% title('gyro1, d vergence');
% 
% subplot(2,3,5)
% [corr lags]=(nanxcorr(allGyro2,dlcDphi,frRate,'coeff'));
% plot(lags, corr)
% axis square;ylim([-1 1]);
% title('gyro2, d phi');
% 
% subplot(2,3,6)
% [corr lags]=(nanxcorr(allGyro3(test==0),dlcDhth(test==0),frRate,'coeff'));
% plot(lags, corr); hold on;
% % [corr lags]=(nanxcorr(allGyro3(test==0),dlcDhth(test==0),frRate,'coeff'));
% % plot(lags, corr,'r')
% axis square;ylim([-1 1]);
% title('gyro ch 3,d head')
% 
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%


for vid =1:length(useFilt)
    
   if length(thetaR{useFilt(vid)})>=100
tr=medfilt1(thetaR{useFilt(vid)}(1:(2*frRate)),8);
pr=medfilt1(phiR{useFilt(vid)}(1:(2*frRate)),8);

tl=medfilt1(thetaL{useFilt(vid)}(1:(2*frRate)),8);
pl=medfilt1(phiL{useFilt(vid)}(1:(2*frRate)),8);

figure
subplot(1,2,1)
plot(tr-nanmean(tr),pr-nanmean(pr),'k'); axis square; hold on; xlabel('eye theta'); ylabel('eye phi'); axis equal
title('right eye'); xlim([-40 40]); ylim([-40 40]);

subplot(1,2,2)
plot(tl-nanmean(tl),pl-nanmean(pl),'k'); axis square; hold on; xlabel('eye theta'); ylabel('eye phi'); axis equal
title('left eye'); xlim([-40 40]); ylim([-40 40]);

for i =1:length(tl)
subplot(1,2,1)
plot(tr(i)-nanmean(tr),pr(i)-nanmean(pr),'.','Markersize',15,'Color', cmapVar(i,1,length(tl),jet)); axis square; hold on
subplot(1,2,2)
plot(tl(i)-nanmean(tl),pl(i)-nanmean(pl),'.','Markersize',15,'Color', cmapVar(i,1,length(tl),jet)); axis square; hold on
 end


if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
   
   end
end

%%
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid = 1:length(useFilt)
%     subplot(rownum,colnum,vid)  ;
%     bar([mean(isnan(mouse_xy{useFilt(vid),1}(1,:))) mean(isnan(cricket_xy{useFilt(vid),1}(1,:))) mean(isnan(crickH{useFilt(vid),1}(1,:)))])
%     ylabel('% error'); xlim([0.5 3.5]); ylim([0 1])
%     set(gca,'XTick',[1 2 3])
%     set(gca,'XTickLabel',{'mouse','c body', 'c head'})
% end
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%
% figure
% for vid=1:length(useFilt)
%     clear approach speed dec heading nframe d2cr dDist
%    % subplot(rownum,colnum,vid)
%   % subplot(7,6,vid)
%
%     d2cr=(range{useFilt(vid)});
%     dDist=diff(d2cr);
% %  dDist=interpNan(dDist,3,'linear')
% %  dDist=medfilt2(dDist)
%     distThresh= d2cr<10; %long approaches can happen at up to 60cm away from cricket, which is diagonal through entire arena
%     dec = (dDist<-.2) %threshold for change in range to cricket
%     speed= mouseV{useFilt(vid)}>=5
%     az=azDeg{useFilt(vid)}; %too many missing points - heading looks weird
%     %azFilt=interpNan(az,3,'linear');
%     heading=(az<90) &(az>-90) ;
%     nframe=min(length(speed),length(heading));
%     nframe=min(nframe,length(dec));
%     speed=speed(1:nframe)';heading=heading(1:nframe)';dec=dec(1:nframe)'; distThresh=distThresh(1:nframe)';
%     approach =dec==1&(speed==1) &(distThresh==1) &heading==1
% %     plot(d2cr,'b'); hold on;
% %     plot(find(approach),d2cr(approach),'og');
%     appEpoch{vid,:}=(approach); axis square
%    % filtApp(vid,:)=sum(appEpoch{vid})
% %   plot(d2cr(mouseV{useFilt(vid)}>5),'g');
% end


%% identify approach!!!

for vid=1:length(useFilt)
    deltaR = diff(range{useFilt(vid)})*frRate;
    vsmooth = conv(mouseV{useFilt(vid)},ones(5,1)/5,'same');
    dRThresh=-10; %%%cm/sec
    vThresh=10;
    azThresh = pi/4;  %%% pi/4 = 45 deg
    approach = deltaR<dRThresh & vsmooth(1:end-1)>vThresh & abs(azT{useFilt(vid)}(1:end-1))<azThresh;
    approach(1)=0; approach(end)=0; %%% boundary conditions
    
    starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% find start;stop
    
    for j= 1:length(ends)-1;  %%% stitch together approachs with gap of less than 5 frames
        if (starts(j+1)-ends(j))<5 & (range{useFilt(vid)}(starts(j+1))- range{useFilt(vid)}(ends(j)))<3
            approach(ends(j) : starts(j+1))=1;
        end
    end
    starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% update start;stop
    
    for j = 1:length(starts);   %%% remove short approaches (less than 10 frames)
        if ends(j)-starts(j)<10
            approach(starts(j):ends(j))=0;
        end
    end
    starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% update start;stop
    appEpoch{vid,:}=approach;
end


% %%
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     clear dur
%     subplot(rownum,colnum,vid)
%     plot(mouse_xy{useFilt(vid),1}(1,:),mouse_xy{useFilt(vid),1}(2,:)); hold on;
%     plot(cricket_xy{useFilt(vid),1}(1,:),cricket_xy{useFilt(vid),1}(2,:))
%     dur=num2str(length(mouse_xy{useFilt(vid),1}(1,:))/frRate,'%.2f');
%     title([dur,'sec']);
%     xlim([300 1600]);
%     use=appEpoch{vid}==1
%     plot(mouse_xy{useFilt(vid),1}(1,use),mouse_xy{useFilt(vid),1}(2,use),'g'); hold on;
%     end
%
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% %%
% % clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll corrRA corrLA lagsRA lagsLA
% % figure('units','normalized','outerposition',[0 0 1 1])
% % for vid=1:length(useFilt)
% %     subplot(rownum,colnum,vid);
% %     nframe = min(length(range{useFilt(vid)}),length(RRad{useFilt(vid)}));
% %     nframe = min(nframe, length(LRad{useFilt(vid)}));
% %     nframe = min(nframe, size(appEpoch{vid},2))
% %     nonapp=squeeze(appEpoch{vid}(1:nframe)==0)
% %     dT=range{useFilt(vid)}(1:nframe); rR=(RRad{useFilt(vid)}(1:nframe))'; rL=LRad{useFilt(vid)}(1:nframe);
% %     clear use
% %     use =(nonapp==1)% & ~isnan(dT(1:nframe))
% %     if sum(use)>3 & (sum(~isnan(dT(use)))>20)
% %     [corrR lagsR]= nanxcorr(dT(use),rR(use),frRate,'zero');
% %     plot(lagsR/frRate,corrR,'b');xlim([-.3 .3])
% %     hold on;
% %     uselagsR=(lagsR>=-frRate& lagsR<=frRate);
% %
% %     [corrL lagsL]= nanxcorr(dT(use),rL(use),frRate,'zero');
% %     plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
% %     uselagsL=(lagsL>=-frRate & lagsL<=frRate);
% %     else
% %     end
% %     clear use;
% %     use=appEpoch{vid}==1% & ~isnan(dT(1:nframe))
% %     if sum(use)>4% & (sum(isnan(rR(use)))>20)
% %     [corrRA lagsRA]= nanxcorr(dT(use),rR(use),frRate,'zero');
% %     plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
% %     uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
% %
% %     [corrLA lagsLA]= nanxcorr(dT(use),rL(use),frRate,'zero');
% %     plot(lagsLA/frRate,corrLA,'c');xlim([-.3 .3]);
% %     uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
% %     else
% %     end
% %     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
% %         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
% %         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
% %
% %     else
% %     end
% % end
% % suptitle('range to cricket & pupil');
% % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %
% % figure('units','normalized','outerposition',[0 0 1 1])
% % errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% % errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
% %
% % shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% % shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% % shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% % shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
% %
% % plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]);
% % xlim([21 41]); axis square
% % L(1) = plot(nan, nan, 'b-');
% % L(2) = plot(nan, nan, 'r-');
% % L(3) = plot(nan, nan, 'g-');
% % L(4) = plot(nan, nan, 'c-');
% %
% % legend(L,{'R pupil','L pupil','R pupil Approach','L pupil Approach'}); title('range and pupil, both eyes');
% % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% % %%
% % clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% % figure('units','normalized','outerposition',[0 0 1 1])
% % for vid=1:length(useFilt)
% %     subplot(rownum,colnum,vid);
% %     nframe = min(length(cricketV{useFilt(vid)}),length(RRad{useFilt(vid)}));
% %     nframe = min(nframe, length(LRad{useFilt(vid)}));
% %     nframe = min(nframe, size(appEpoch{vid},2))
% %     nonapp=squeeze(appEpoch{vid}(1:nframe)==0)
% %     dT=cricketV{useFilt(vid)}(1:nframe); rR=RRad{useFilt(vid)}(1:nframe); rL=LRad{useFilt(vid)}(1:nframe);
% %     clear use
% %     use =(nonapp==1)%&~isnan(dT(1:nframe));
% %    if sum(use)>3 & sum(~isnan(dT(use)))>20
% %     [corrR lagsR]= nanxcorr(dT(use),rR(use),frRate,'zero');
% %     plot(lagsR/frRate,corrR,'b');xlim([-.3 .3])
% %     hold on;
% %     uselagsR=(lagsR>=-frRate& lagsR<=frRate);
% %
% %     [corrL lagsL]= nanxcorr(dT(use),rL(use),frRate,'zero');
% %     plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
% %     uselagsL=(lagsL>=-frRate & lagsL<=frRate);
% %     else
% %     end
% %     clear use;
% %     use=appEpoch{vid}==1% & ~isnan(dT(1:nframe))
% %     if sum(use)>4 & sum(isnan(rR(use))>20);
% %     [corrRA lagsRA]= nanxcorr(dT(use),rR(use),frRate,'zero');
% %     plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
% %     uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
% %     [corrLA lagsLA]= nanxcorr(dT(use),rL(use),frRate,'zero');
% %     plot(lagsLA/frRate,corrLA,'c');xlim([-.3 .3]);
% %     uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
% %     else
% %     end
% %     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
% %         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
% %         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
% %
% %     else
% %     end
% % end
% % suptitle('cricketV and pupil')
% % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %
% % figure('units','normalized','outerposition',[0 0 1 1])
% % errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% % errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
% %
% % shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% % shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);hold on;
% % shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% % shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);hold on
% %
% % plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]);
% % xlim([21 41]); axis square
% % L(1) = plot(nan, nan, 'b-');
% % L(2) = plot(nan, nan, 'r-');
% % L(3) = plot(nan, nan, 'g-');
% % L(4) = plot(nan, nan, 'c-');
% %
% % legend(L,{'R pupil','L pupil','R pupil Approach','L pupil Approach'}); title('cricket V and pupil, both eyes');
% % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% % %%
% %
% % clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% % figure('units','normalized','outerposition',[0 0 1 1])
% % for vid=1:length(useFilt)
% %     subplot(rownum,colnum,vid);
% %     nframe = min(length(mouseV{useFilt(vid)}),length(RRad{useFilt(vid)}));
% %     nframe = min(nframe, length(LRad{useFilt(vid)}));
% %     nframe = min(nframe, size(appEpoch{vid},2))
% %     nonapp=squeeze(appEpoch{vid}(1:nframe)==0)
% %     dT=mouseV{useFilt(vid)}(1:nframe); rR=RRad{useFilt(vid)}(1:nframe); rL=LRad{useFilt(vid)}(1:nframe);
% %     clear use
% %     use =  (nonapp==1) %& ~isnan(dT(1:nframe));
% %   %  if sum(use)>3
% %     [corrR lagsR]= nanxcorr(dT(use),rR(use),frRate,'zero');
% %     plot(lagsR/frRate,corrR,'b');xlim([-.3 .3])
% %     hold on;
% %     uselagsR=(lagsR>=-frRate& lagsR<=frRate);
% %
% %     [corrL lagsL]= nanxcorr(dT(use),rL(use),frRate,'zero');
% %     plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
% %     uselagsL=(lagsL>=-frRate & lagsL<=frRate);
% % %     else
% % %     end
% %     clear use;
% %     use=appEpoch{vid}==1 %& ~isnan(dT(1:nframe))
% %     if sum(use)>4
% %     [corrRA lagsRA]= nanxcorr(dT(use),rR(use),frRate,'zero');
% %     plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
% %     uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
% %     [corrLA lagsLA]= nanxcorr(dT(use),rL(use),frRate,'zero');
% %     plot(lagsLA/frRate,corrLA,'c');xlim([-.3 .3]);
% %     uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
% %     else
% %     end
% %     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
% %         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
% %         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
% %
% %     else
% %     end
% % end
% % suptitle('mouseV and pupil')
% % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %
% % figure('units','normalized','outerposition',[0 0 1 1])
% % errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% % errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
% %
% % shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% % shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% % shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% % shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
% %
% % plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]);
% % xlim([21 41]); axis square
% % L(1) = plot(nan, nan, 'b-');
% % L(2) = plot(nan, nan, 'r-');
% % L(3) = plot(nan, nan, 'g-');
% % L(4) = plot(nan, nan, 'c-');
% %
% % legend(L,{'R pupil','L pupil','R pupil Approach','L pupil Approach'}); title('mouse V and pupil, both eyes');
% % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %
% %
% %
%
% %%
% % clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% % figure('units','normalized','outerposition',[0 0 1 1])
% % for vid=1:length(useFilt)
% %     clear uselagsRA
% %     subplot(rownum,colnum,vid);
% %     nframe = min(length(dTheta{useFilt(vid)}),length(dthetaR{useFilt(vid)}));
% %     nframe = min(nframe, length(dthetaL{useFilt(vid)}));
% %     nonapp=appEpoch{vid}(1:nframe)==0
% %     dT=dTheta{useFilt(vid)}; dtR=dthetaR{useFilt(vid)}; dtL=dthetaL{useFilt(vid)};
% %     clear use
% %     use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
% %   %  if sum(use)>3
% %     [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
% %     plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
% %     hold on;
% %     uselagsR=(lagsR>=-frRate& lagsR<=frRate);
% %
% %     [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
% %     plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
% %     uselagsL=(lagsL>=-frRate & lagsL<=frRate);
% % %     else
% % %     end
% %     clear use;
% %     use=appEpoch{vid}==1
% %     if sum(use)>4 & sum(~isnan(dtR(use)))>20
% %     [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
% %     plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
% %     uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
% %
% %     else
% %         corrRA=NaN;uselagsRA=NaN;
% %
% %
% %     end
% %      if sum(use)>4 & sum(~isnan(dtL(use)))>20
% %     [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
% %     plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
% %     uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
% %     else
% %     end
% %     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
% %         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
% %         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
% %
% %     else
% %
% %     end
% % end
% %
% %
% % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %
% % figure('units','normalized','outerposition',[0 0 1 1])
% % errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% % errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
% %
% % shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% % shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% % shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% % shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
% %
% % plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
% % L(1) = plot(nan, nan, 'b-');
% % L(2) = plot(nan, nan, 'r-');
% % L(3) = plot(nan, nan, 'g-');
% % L(4) = plot(nan, nan, 'c-');
% %
% % legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('head theta, both eyes');
% % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %%
% % clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% % figure('units','normalized','outerposition',[0 0 1 1])
% % for vid=1:length(useFilt)
% %     subplot(rownum,colnum,vid);
% % %     nframe = min(length(dTheta{useFilt(vid)}),length(dthetaR{useFilt(vid)}));
% % %     nframe = min(nframe, length(dthetaL{useFilt(vid)}));
% %     nonapp=appEpoch{vid}==0
% %     dT=dTheta{useFilt(vid)}; dpR=dphiR{useFilt(vid)}; dpL=dphiL{useFilt(vid)};
% %     clear use
% %     use = (nonapp==1)'% & ~isnan(dT(1:length(dpR)));
% %     if sum(use)>3 &sum(~isnan(dpR(use)))>20
% %     [corrR lagsR]= nanxcorr(dT(use),dpR(use),frRate,'coeff');
% %     plot(lagsR/frRate,corrR,'b');xlim([-.3 .3]);hold on;
% %     uselagsR=(lagsR>=-frRate& lagsR<=frRate);
% %
% %     [corrL lagsL]= nanxcorr(dT(use),dpL(use),frRate,'coeff');
% %     plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
% %     uselagsL=(lagsL>=-frRate & lagsL<=frRate);
% %     else
% %     end
% %     clear use;
% %     use=appEpoch{vid}==1
% %     if sum(use)>3 & sum(~isnan(dpR(use)))>20
% %     [corrRA lagsRA]= nanxcorr(dT(use),dpR(use),frRate,'coeff');
% %     plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);hold on;
% %     uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
% %     [corrLA lagsLA]= nanxcorr(dT(use),dpL(use),frRate,'coeff');
% %     plot(lagsLA/frRate,corrLA,'c');xlim([-.3 .3]);
% %     uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
% %     else
% %         corrRA=NaN; uselagsRA=NaN;
% %     end
% %     if sum(uselagsR)==2*(frRate)+1) & sum(uselagsL)==2*(frRate)+1) &sum(uselagsRA)==2*(frRate)+1) & sum(uselagsLA)==2*(frRate)+1)
% %         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
% %         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
% %
% %     else
% %      corrRAAll(vid,:)=NaN; corrLAAll(vid,:)=NaN;
% %
% %     end
% % end
% % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %
% % figure('units','normalized','outerposition',[0 0 1 1])
% % errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% % errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
% %
% % shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% % shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% % shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% % shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
% %
% % plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
% % L(1) = plot(nan, nan, 'b-');
% % L(2) = plot(nan, nan, 'r-');
% % L(3) = plot(nan, nan, 'g-');
% % L(4) = plot(nan, nan, 'c-');
% %
% % legend(L,{'dPhi R non-app','dPhi L non-app','dPhi R Approach','dPhi L Approach'}); title('head Theta and Eye Phi, both eyes');
% % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll corrRA lagsRA corrLA lagsLA
% corrRAAll=[]
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
% clear uselagsRA corrRA
%     subplot(rownum,colnum,vid);
%     nframe = min(length(dthetaR{useFilt(vid)}),length(dthetaL{useFilt(vid)}));
%     dtR=dthetaR{useFilt(vid)}(1:nframe); dtL=dthetaL{useFilt(vid)}(1:nframe);
%     nonapp=appEpoch{vid}==0;
%     use = (nonapp==1)'%& ~isnan(dtR(1:nframe));
%     if sum(use)>3
%     [corrR lagsR]= nanxcorr(dtR(use),dtL(use),frRate,'coeff');
%     plot(lagsR/frRate,corrR,'b');xlim([-.3 .3]);hold on;
%     uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%     clear nframe
%     nframe = min(length(dphiR{useFilt(vid)}),length(dphiL{useFilt(vid)}));
%     dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
%     [corrL lagsL]= nanxcorr(dpR(use),dpL(use),frRate,'coeff');
%     plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%     uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%       % use = ~isnan(dT(1:nframe));
%        use=appEpoch{vid}==1;
%      if sum(use)>3 & sum(~isnan(dtR(use)))>20
%     [corrRA lagsRA]= nanxcorr(dtR(use),dtL(use),frRate,'coeff');
%     plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);hold on;
%     uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%
% %     clear nframe
%     nframe = min(length(dphiR{useFilt(vid)}),length(dphiL{useFilt(vid)}));
%     dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
%     [corrLA lagsLA]= nanxcorr(dpR(use),dpL(use),frRate,'coeff');
%     plot(lagsLA/frRate,corrLA,'c');xlim([-.3 .3]);
%     uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%      else
%        %  corrRA=NaN;
%        uselagsRA=NaN;
%      end
%
%     if sum(uselagsR)==2*(frRate)+1) & sum(uselagsL)==2*(frRate)+1) &sum(uselagsRA)==2*(frRate)+1) & sum(uselagsLA)==2*(frRate)+1)
%        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
% suptitle('between eye corr, dtheta & dphi - non app & approach');
%       if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
%       figure('units','normalized','outerposition',[0 0 1 1])
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]);
% xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'dTheta non-app','dPhi non-app','dTh Approach','dPhi app'}); title('mean between eye corr, dtheta & dphi');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% %%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% corrRAAll=[]
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     clear corrRA lagsRA
%
%     subplot(rownum,colnum,vid);
%     nframe = min(length(dthetaR{useFilt(vid)}),length(dphiR{useFilt(vid)}));
%     dtR=dthetaR{useFilt(vid)}(1:nframe); dpR=dphiR{useFilt(vid)}(1:nframe);
%     nonapp=appEpoch{vid}==0;
%     use = (nonapp==1)'%& ~isnan(dtR(1:nframe));
%     if sum(use)>3 &sum(~isnan(dpR(use)))>20
%     [corrR lagsR]= nanxcorr(dtR(use),dpR(use),frRate,'coeff');
%     plot(lagsR/frRate,corrR,'b');xlim([-.3 .3]);hold on;
%     uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%     clear nframe
%     nframe = min(length(dthetaL{useFilt(vid)}),length(dphiL{useFilt(vid)}));
%     dtL=dthetaL{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
%     [corrL lagsL]= nanxcorr(dtL(use),dpL(use),frRate,'coeff');
%     plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%     uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%       % use = ~isnan(dT(1:nframe));
%       clear use
%       use=appEpoch{vid}==1;
%      if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dpR(use)))>20
%     [corrRA lagsRA]= nanxcorr(dtR(use),dpR(use),frRate,'coeff');
%     plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);hold on;
%     uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%
% %     clear nframe
%     nframe = min(length(dphiR{useFilt(vid)}),length(dphiL{useFilt(vid)}));
%     dtL=dthetaL{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
%     [corrLA lagsLA]= nanxcorr(dtL(use),dpL(use),frRate,'coeff');
%     plot(lagsLA/frRate,corrLA,'c');xlim([-.3 .3]);
%     uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%      else
%          corrRA =NaN; uselagsRA=NaN;
%      end
%
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
% suptitle('dTheta vs dPhi corr, R and L eyes');
%       if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
%       figure('units','normalized','outerposition',[0 0 1 1])
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]);
% xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'R dTh/dPhi','L dTh/dPhi','R dTh/dPhi Approach','L dTh/dPhi Approach'}); title('mean dtheta & dphi, each eye');
%
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %%
%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     subplot(rownum,colnum,vid);
%     nonapp=appEpoch{vid}==0
%     nframe = min(length(dTheta{useFilt(vid)}),length(thetaR{useFilt(vid)}));
%     nframe = min(nframe, length(thetaL{useFilt(vid)}));
%     nframe=min(nframe,length(nonapp));
%
%     dT=dTheta{useFilt(vid)}(1:nframe); tR=thetaR{useFilt(vid)}(1:nframe); tL=thetaL{useFilt(vid)}(1:nframe);
%     clear use
%     use = ~isnan(dT(1:nframe))& (nonapp==1)';
%     if sum(use)>3
%     [corrR lagsR]= nanxcorr(dT(use),tR(use),frRate,'coeff');
%     plot(lagsR/frRate,corrR,'b');xlim([-.3 .3]);hold on;
%     uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%     [corrL lagsL]= nanxcorr(dT(use),tL(use),frRate,'coeff');
%     plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%     uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=(appEpoch{vid}==1)' %& ~isnan(dT(1:nframe));
%     if sum(use)>3 & sum(~isnan(dT(use)))>20
%       [corrRA lagsRA]= nanxcorr(dT(use),tR(use),frRate,'coeff');
%     plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);hold on;
%     uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%
%     [corrLA lagsLA]= nanxcorr(dT(use),tL(use),frRate,'coeff');
%     plot(lagsLA/frRate,corrLA,'c');xlim([-.3 .3]);
%     uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%
%     else
%     end
%
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
% title('dHead & theta Position');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% figure('units','normalized','outerposition',[0 0 1 1])
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]);
% xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
% title('dHead and theta position');
% legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'});
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %%
%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     subplot(rownum,colnum,vid);
%     nonapp=appEpoch{vid}==0
%     nframe = min(length(dTheta{useFilt(vid)}),length(phiR{useFilt(vid)}));
%     nframe = min(nframe, length(phiL{useFilt(vid)}));
%     nframe=min(nframe,length(nonapp));
%
%     dT=dTheta{useFilt(vid)}(1:nframe); pR=phiR{useFilt(vid)}(1:nframe); pL=phiL{useFilt(vid)}(1:nframe);
%     clear use
%     use = (nonapp==1)';%~isnan(dT(1:nframe))& (nonapp==1)';
%     if sum(use)>3
%     [corrR lagsR]= nanxcorr(dT(use),pR(use),frRate,'coeff');
%     plot(lagsR/frRate,corrR,'b');xlim([-.3 .3]);hold on;
%     uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%     [corrL lagsL]= nanxcorr(dT(use),pL(use),frRate,'coeff');
%     plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%     uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=(appEpoch{vid}==1)' %& ~isnan(dT(1:nframe));
%     if sum(use)>3 & sum(~isnan(dT(use)))>20
%       [corrRA lagsRA]= nanxcorr(dT(use),pR(use),frRate,'coeff');
%     plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);hold on;
%     uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%
%     [corrLA lagsLA]= nanxcorr(dT(use),pL(use),frRate,'coeff');
%     plot(lagsLA/frRate,corrLA,'c');xlim([-.3 .3]);
%     uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%
%     else
%     end
%
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
% title('dHead & phi Position');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% figure('units','normalized','outerposition',[0 0 1 1])
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]);
% xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
% title('dHead and phi position');
% legend(L,{'phi R non-app','phi L non-app','phi R Approach','phiL Approach'});
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% corrRAAll=[]
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     subplot(rownum,colnum,vid);
%     nframe = min(length(thetaR{useFilt(vid)}),length(thetaL{useFilt(vid)}))-1;
%     tR=thetaR{useFilt(vid)}(1:nframe); tL=thetaL{useFilt(vid)}(1:nframe);
%     nonapp=appEpoch{vid}(1:nframe)==0;
%
%     use = (nonapp==1)' %& ~isnan(tR(1:length(nonapp)));
%     if sum(use)>3
%     [corrR lagsR]= nanxcorr(tR(use),tL(use),frRate,'zero');
%     plot(lagsR/frRate,corrR,'b');xlim([-.3 .3]);hold on;
%     uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%     clear nframe
%     nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
%     pR=phiR{useFilt(vid)}(1:nframe); pL=phiL{useFilt(vid)}(1:nframe);
%     [corrL lagsL]= nanxcorr(pR(use),pL(use),frRate,'zero');
%     plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%     uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%       % use = ~isnan(dT(1:nframe));
%        use=appEpoch{vid}==1;
%      if sum(use)>3 & sum(~isnan(tR(use)))>20
%     [corrRA lagsRA]= nanxcorr(tR(use),tL(use),frRate,'zero');
%     plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);hold on;
%     uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%
% %     clear nframe
%     nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
%     pR=phiR{useFilt(vid)}(1:nframe); pL=phiL{useFilt(vid)}(1:nframe);
%     [corrLA lagsLA]= nanxcorr(pR(use),pL(use),frRate,'zero');
%     plot(lagsLA/frRate,corrLA,'c');xlim([-.3 .3]);
%     uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%
%      else
%      end
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
% title('between eye corr, theta & phi position - non app & approach');
%       if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
%       figure('units','normalized','outerposition',[0 0 1 1])
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-1 1]);
% xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'theta non-app','Phi non-app','Th Approach','Phi app'}); title('mean between eye corr, theta & phi position');
%
%
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% %%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% corrRAAll=[]
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     subplot(rownum,colnum,vid);
%     nframe = min(length(thetaR{useFilt(vid)}),length(phiR{useFilt(vid)}))-1;
%     tR=thetaR{useFilt(vid)}(1:nframe); pR=phiR{useFilt(vid)}(1:nframe);
%     nonapp=appEpoch{vid}(1:nframe)==0;
%
%     use = (nonapp==1)' %& ~isnan(tR(1:length(nonapp)));
%     if sum(use)>3
%     [corrR lagsR]= nanxcorr(tR(use),pR(use),frRate,'zero');
%     plot(lagsR/frRate,corrR,'b');xlim([-.3 .3]);hold on;
%     uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%     clear nframe
%     nframe = min(length(thetaL{useFilt(vid)}),length(phiL{useFilt(vid)}));
%     tL=thetaL{useFilt(vid)}(1:nframe); pL=phiL{useFilt(vid)}(1:nframe);
%     [corrL lagsL]= nanxcorr(tL(use),pL(use),frRate,'zero');
%     plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%     uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     % use = ~isnan(dT(1:nframe));
%     clear use
%     use=appEpoch{vid}==1;
%          if sum(use)>3 &sum(~isnan(tR(use)))>20
%     [corrRA lagsRA]= nanxcorr(tR(use),pR(use),frRate,'zero');
%     plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);hold on;
%     uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%
% %     clear nframe
%     nframe = min(length(thetaL{useFilt(vid)}),length(phiL{useFilt(vid)}));
%     tL=thetaL{useFilt(vid)}(1:nframe); pL=phiL{useFilt(vid)}(1:nframe);
%     [corrLA lagsLA]= nanxcorr(tL(use),pL(use),frRate,'zero');
%     plot(lagsLA/frRate,corrLA,'c');xlim([-.3 .3]);
%     uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%          else
%          end
%
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
% suptitle('theta vs phi corr each eye - non app & approach');
%       if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
%       figure('units','normalized','outerposition',[0 0 1 1])
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-1 1]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'theta & phi R eye','theta & phi L eye','th/phi R app','th/phi L app'}); title('mean theta & phi corr, each eye');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %%
%
% clear dTAll dTRALL dTLAll dPRAll dPLAll
% close all
% dTRAll=[];
%
% % figure%('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     clear dT dpR dpL dtR dtL
%     %     subplot(rownum,colnum,vid);
%     nframe = min(length(dTheta{useFilt(vid)}),length(dphiR{useFilt(vid)}));
%     nframe = min(nframe, length(dphiL{useFilt(vid)}));
%     nframe = min(nframe, length(dthetaR{useFilt(vid)}));
%     nframe = min(nframe,length(dthetaL{useFilt(vid)}));
%     dT=dTheta{useFilt(vid)}(1:nframe); dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
%     dtR=dthetaR{useFilt(vid)}(1:nframe); dtL=dthetaL{useFilt(vid)}(1:nframe);
%     useN= appEpoch{vid}==0;
%     use = (appEpoch{vid});
%
%     figure(1);
%     plot(dT(useN(1:15:end)),dpR(useN(1:15:end)),'b.');axis square; hold on
%     title('dHead theta, dPhi R, nonapp & app');
%     x = linspace(-40,40);
%     y = linspace(-60,60);
%     plot(-x,y); xlim([-40 40]); ylim([-60 60]);
%     plot(dT(use(1:15:end)),dpR(use(1:15:end)),'.g');
%     if vid==(useFilt(end))
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(2);
%     plot(dT(useN(1:15:end)),dpL(useN(1:15:end)),'b.'); axis square; hold on
%     title('dHead theta, dPhi L, nonapp & app');
%     x = linspace(-40,40);
%     y = linspace(-80,80);
%     plot(x,y);  xlim([-40 40]); ylim([-80 80]);
%     plot(dT(use(1:15:end)),dpL(use(1:15:end)),'.g');
%     if vid==(useFilt(end))
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(3);
%     plot(dT(useN(1:15:end)),dtR(useN(1:15:end)),'b.');axis square; hold on
%     title('dHead theta, dtheta R, nonapp & app');
%     x = linspace(-40,40);
%     y = linspace(-40,40);  xlim([-40 40]); ylim([-40 40]);
%     plot(-x,y);
%     plot(dT(use(1:15:end)),dtR(use(1:15:end)),'.g');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(4);
%     plot(dT(useN(1:15:end)),dtL(useN(1:15:end)),'b.');axis square; hold on;
%     title('dHead theta, dtheta L, nonapp & app');
%     x = linspace(-40,40);  xlim([-40 40]); ylim([-80 80]);
%     y = linspace(-80,80);
%     plot(x,y);
%     plot(dT(use(1:15:end)),dtL(use(1:15:end)),'.g');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(5);%subplot(rownum,colnum,vid);
%     plot(dtR(useN(1:15:end)),dtL(useN(1:15:end)),'b.');axis square; hold on;
%     title('r theta, l theta, nonapp & app');
%     x = linspace(-80,80);
%     y = linspace(-50,50);
%     plot(-x,y); xlim([-80 80]); ylim([-50 50]);
%     plot(dtR(use(1:15:end)),dtL(use(1:15:end)),'.g');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(6);
%     plot(dpR(useN(1:15:end)),dpL(useN(1:15:end)),'b.');axis square; hold on;
%     x = linspace(-80,80);
%     y = linspace(-80,80);
%     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
%     title('right phi, left phi, nonapp & app');
%     plot(dpR(use(1:15:end)),dpL(use(1:15:end)),'.g');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%
%      figure(7);
%     plot(dtR(useN(1:15:end)),dpR(useN(1:15:end)),'b.');axis square; hold on;
%     x = linspace(-80,80);
%     y = linspace(-80,80);
%     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
%     title('R eye dTheta dPhi - nonapp & app');
%     plot(dtR(use(1:15:end)),dpR(use(1:15:end)),'.g');
%        if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%
%     figure(8);
%     plot(dtL(useN(1:15:end)),dpL(useN(1:15:end)),'b.');axis square; hold on;
%     x = linspace(-80,80);
%     y = linspace(-80,80);
%     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
%     title('L eye dTheta dPhi - nonapp & app');
%     plot(dtL(use(1:15:end)),dpL(use(1:15:end)),'.g');
%        if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%
% end
%
%
%
% %%
% close all
%
% % figure%('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     clear dT dpR dpL dtR dtL
%     %     subplot(rownum,colnum,vid);
%     nframe = min(length(range{useFilt(vid)}),length(RRad{useFilt(vid)}));
%     nframe = min(nframe, length(LRad{useFilt(vid)}));
%     nframe = min(nframe, length(mouseV{useFilt(vid)}));
%     nframe = min(nframe,length(cricketV{useFilt(vid)}));
%     nframe = min(nframe,length(azT{useFilt(vid)}));
%
%    r=range{useFilt(vid)}(1:nframe); rR=RRad{useFilt(vid)}(1:nframe); lR=LRad{useFilt(vid)}(1:nframe);
%     mouseSp=mouseV{useFilt(vid)}(1:nframe); crSp=cricketV{useFilt(vid)}(1:nframe);
%     az=azT{useFilt(vid)}(1:nframe);
%     useN= appEpoch{vid}==0;
%     use = (appEpoch{vid});
%     figure(1);
%     plot(rR(useN(1:15:end)),lR(useN(1:15:end)),'bo');axis square; hold on
%     xlabel('r eye'); ylabel('l eye');
%     title('two eyes, rad');
% %     x = linspace(-40,40);
% %     y = linspace(-60,60);
% %     plot(-x,y);
%     xlim([15 36]); ylim([15 36]);
%     plot(rR(use(1:15:end)),lR(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(2);
%     plot(r(useN(1:15:end)),rR(useN(1:15:end)),'bo'); axis square; hold on
%     title('range and R Rad');
%     xlabel('range to cricket (cm)'); ylabel('R eye');
% %     x = linspace(-40,40);
% %     y = linspace(-80,80);
% %     plot(x,y);  xlim([-40 40]); ylim([-80 80]);
%     plot(r(use(1:15:end)),rR(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(3);
%     plot(r(useN(1:15:end)),lR(useN(1:15:end)),'bo');axis square; hold on
%     title('range and L Rad');
%     xlabel('range to cricket (cm)'); ylabel('L eye');
% %     x = linspace(-40,40);
% %     y = linspace(-40,40);  xlim([-40 40]); ylim([-40 40]);
% %     plot(-x,y);
%     plot(r(use(1:15:end)),lR(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(4);
%     plot(mouseSp(useN(1:15:end)),rR(useN(1:15:end)),'bo');axis square; hold on;
%     title('mouse speed & R rad');
%     xlabel('mouse Speed (cm/sec)'); ylabel('R Rad');
% %     x = linspace(-40,40);  xlim([-40 40]); ylim([-80 80]);
% %     y = linspace(-80,80);
% %     plot(x,y);
%     plot(mouseSp(use(1:15:end)),rR(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(5);%subplot(rownum,colnum,vid);
%     plot(mouseSp(useN(1:15:end)),lR(useN(1:15:end)),'bo');axis square; hold on;
%     title('mouse speed & L rad');
%     xlabel('mouse Speed (cm/sec)'); ylabel('L Rad');
% %     x = linspace(-80,80);
% %     y = linspace(-50,50);
% %     plot(-x,y); xlim([-80 80]); ylim([-50 50]);
%     plot(mouseSp(use(1:15:end)),lR(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(6);
%     plot(crSp(useN(1:15:end)),rR(useN(1:15:end)),'bo');axis square; hold on;
% %     x = linspace(-80,80);
% %     y = linspace(-80,80);
% %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
%     title('cricket speed & R rad');
%     xlabel('cricket Speed'); ylabel('R Rad');
%     plot(crSp(use(1:15:end)),rR(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%
%      figure(7);
%     plot(crSp(useN(1:15:end)),lR(useN(1:15:end)),'bo');axis square; hold on;
% %     x = linspace(-80,80);
% %     y = linspace(-80,80);
% %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
%   title('cricket speed & L rad');
%     xlabel('cricket Speed'); ylabel('l Rad');
%     plot(crSp(use(1:15:end)),lR(use(1:15:end)),'og');
%        if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%
%     figure(8);
%     plot(az(useN(1:15:end)),r(useN(1:15:end)),'b.');axis square; hold on;
% %     x = linspace(-80,80);
% %     y = linspace(-80,80);
% %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
%     title('az and range');
%     plot(az(use(1:15:end)),r(use(1:15:end)),'.g');
%        if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%
% end
%
%
%
%
%
%
% %%
% clear dTAll dTRALL dTLAll dPRAll dPLAll
%
% dTRAll=[];
% close all
% % figure%('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     clear dT dpR dpL dtR dtL
%     %     subplot(rownum,colnum,vid);
%     nframe = min(length(dTheta{useFilt(vid)}),length(phiR{useFilt(vid)}));
%     nframe = min(nframe, length(phiL{useFilt(vid)}));
%     nframe = min(nframe, length(thetaR{useFilt(vid)}));
%     nframe = min(nframe,length(thetaL{useFilt(vid)}));
%     dT=dTheta{useFilt(vid)}(1:nframe); pR=phiR{useFilt(vid)}(1:nframe); pL=phiL{useFilt(vid)}(1:nframe);
%     tR=thetaR{useFilt(vid)}(1:nframe); tL=thetaL{useFilt(vid)}(1:nframe);
%     useN= appEpoch{vid}==0;
%     use = (appEpoch{vid});
%
%
%     figure(1);
%     plot(dT(useN(1:15:end)),pR(useN(1:15:end)),'bo');axis square; hold on
%     title('head dtheta, phi R, nonapp & app');
%     x = linspace(-40,40);
%     y = linspace(-60,60);
%     plot(-x,y); xlim([-40 40]); ylim([-60 60]);
%     plot(dT(use(1:15:end)),pR(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(2);
%     plot(dT(useN(1:15:end)),pL(useN(1:15:end)),'bo'); axis square; hold on
%     title('head dtheta, phi L, nonapp & app');
%     x = linspace(-40,40);
%     y = linspace(-80,80);
%     plot(x,y);  xlim([-40 40]); ylim([-80 80]);
%     plot(dT(use(1:15:end)),pL(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(3);
%     plot(dT(useN(1:15:end)),tR(useN(1:15:end)),'bo');axis square; hold on
%     title('head dtheta, theta R, nonapp & app');
%     x = linspace(-40,40);
%     y = linspace(-40,40);  xlim([-40 40]); ylim([-40 40]);
%     plot(-x,y);
%     plot(dT(use(1:15:end)),tR(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(4);
%     plot(dT(useN(1:15:end)),tL(useN(1:15:end)),'bo');axis square; hold on;
%     title('head dtheta, theta L, nonapp & app');
%     x = linspace(-40,40);  xlim([-40 40]); ylim([-80 80]);
%     y = linspace(-80,80);
%     plot(x,y);
%     plot(dT(use(1:15:end)),tL(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(5);%subplot(rownum,colnum,vid);
%     plot(tR(useN(1:15:end)),tL(useN(1:15:end)),'bo');axis square; hold on;
%     title('r theta, l theta, nonapp & app');
%     x = linspace(-80,80);
%     y = linspace(-50,50);
%     plot(-x,y); xlim([-80 80]); ylim([-50 50]);
%     plot(tR(use(1:15:end)),tL(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(6);
%     plot(pR(useN(1:15:end)),pL(useN(1:15:end)),'bo');axis square; hold on;
%     x = linspace(-80,80);
%     y = linspace(-80,80);
%     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
%     title('right phi, left phi, nonapp & app');
%     plot(pR(use(1:15:end)),pL(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%
%      figure(7);
%     plot(tR(useN(1:15:end)),pR(useN(1:15:end)),'bo');axis square; hold on;
%     x = linspace(-80,80);
%     y = linspace(-80,80);
%     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
%     title('right eye - theta/phi pos - nonapp & app');
%     plot(tR(use(1:15:end)),pR(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%
%          figure(8);
%     plot(tL(useN(1:15:end)),pL(useN(1:15:end)),'bo');axis square; hold on;
%     x = linspace(-80,80);
%     y = linspace(-80,80);
%     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
%     title('left eye - theta/phi pos - nonapp & app');
%     plot(tL(use(1:15:end)),pL(use(1:15:end)),'og');
%     if vid==(useFilt(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%
% end

%%
pSname='T:\PreyCaptureAnalysis\Data\';
if savePDF
    filen=sprintf('%s',ani,'Deinterlaced_052320_analyzed_halfShift','.pdf')
    pdfilename=fullfile(pSname,filen);
    dos(['ps2pdf ' psfilename ' ' pdfilename]);
    delete(psfilename);
    
end


afilename=sprintf('%s',ani,'Deinterlaced_052320_analyzed_halfShift','.mat');
save(fullfile(pSname, afilename))