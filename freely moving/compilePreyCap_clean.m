clear all; close all
load('J475cAllVids.mat'); 
savePDF=1;
if savePDF
    psfilename = 'C:\analysisPS.ps';
    if exist(psfilename,'file')==2;delete(psfilename);end
end


% 
% Correlations:
% 
% DELTA
% - head & phi
%     - w/app
%     - w/app & no app
% - head & theta
%     - w/app
%     - w/ app & no app
% - both eyes theta & phi
%     - w app
%     - w app & no app
%     
% Phi & theta (no change)
% - head & phi
%     - w/app
%     - w/app & no app
% - head & theta
%     - w/app
%     - w/ app & no app
% - both eyes theta & phi
%     - w app
%     - w app & no app
%     
% SCATTER PLOTS
%    DELTA
%    - w/app
%    - w/app & no app
%    
%    ABS Position
%    - w/app
%    - w/app & no app

mouse_xy=[];cricket_xy=[];EllipseParamsR=[];EllipseParamsL=[]; radR=[]; az=[];

for i=1:length(Data)
    
    mouse_xyRaw{i,1,1}=Data(i).mouse_xyRaw;
    mouseVRaw{i,:}= Data(i).mouseVRaw; %in pix/frame
    mouseVRaw{i,:}=((Data(i).mouseVRaw)/27)*30; %cm/sec
    cricket_xyRaw{i,1}=Data(i).cricketxyRaw;
    cricketVRaw{i,:}= Data(i).cricketVRaw % pix/frame
    cricketVRaw{i,:} = ((Data(i).cricketVRaw)/27)*30; %now cm/sec
    thetaRaw{i,:}= Data(i).thetaRaw;
    dThetaRaw{i,:}= diff(Data(i).thetaRaw);
    rangeRaw{i,:}= (Data(i).rangeRaw)/27; %convert pixels to cm
    azRaw = mod(Data(i).azRaw,2*pi); azRaw(azRaw>pi) = azRaw(azRaw>pi)-2*pi;
    azTRaw{i,:,:}= azRaw;
    azDegRaw{i,:,:}=rad2deg(azTRaw{i,:,:});%az in deg
    cricketBRaw{i,1,1} = Data(i).cricketBody;
    cricketB{i,1,1} = Data(i).cricketBody;
            
    mouse_xy{i,1,1}=Data(i).mouse_xy;
    mouseV{i,:}= Data(i).mouseV; %in pix/frame
    mouseV{i,:}=((Data(i).mouseV)/27)*30; %cm/sec
    cricket_xy{i,1}=Data(i).cricketxy;
    cricketV{i,:}= Data(i).cricketV % pix/frame
    cricketV{i,:} = ((Data(i).cricketV)/27)*30; %now cm/sec
    theta{i,:}= rad2deg(Data(i).theta);
    dTheta{i,:}=rad2deg(Data(i).dth);
    range{i,:}= (Data(i).range)/27; %convert pixels to cm
    az = mod(Data(i).az,2*pi); az(az>pi) = az(az>pi)-2*pi;
    azT{i,:,:}= az;
    azDeg{i,:,:}=rad2deg(azT{i,:,:});%az in deg
    slip{i,1} = Data(i).difTR;
    slip{i,2} = Data(i).difTL;
    slip{i,3} = Data(i).difRL;
    goodTheta(i,:) = Data(i).ThetaFract; %from un-interpolated theta

    EllipseParamsR{i,1,1} = Data(i).EllipseParamsR;
    EllipseParamsL{i,1,1} = Data(i).EllipseParamsL;
    thetaR{i,:} =Data(i).Rtheta;
    thetaL{i,:} =Data(i).Ltheta;
    phiR{i,:} =Data(i).Rphi;
    phiL{i,:} =Data(i).Lphi;
    
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
    
    
    longR{i,1}=EllipseParamsR{i,1}(:,3);longL{i,1}=EllipseParamsL{i,1}(:,3);
    shortR{i,1}=EllipseParamsR{i,1}(:,4);shortL{i,1}=EllipseParamsL{i,1}(:,4);
    
    radR{i,1}= (longR{i,1}+shortR{i,1})./length(longR{i,1});
    radL{i,1}= (longL{i,1}+shortL{i,1})./length(longL{i,1});
    pupilRvel{i,1}=diff(radR{i,1}); pupilLvel{i,1}=diff(radL{i,1})
    
   % tsData{i,1}=~isempty(Data(i).TopTs);
    
end

% tsData= cell2mat(tsData)
delayFull=cell2mat(slip);

% useL = (delayFull(:,2)<=3 & delayFull(:,2)>=-3);
% useR = (delayFull(:,1)<=3 & delayFull(:,1)>=-3);
% useE = (delayFull(:,3)<=3 & delayFull(:,3)>=-3);

useTime = goodTheta>=.7 %& tsData==1; %|(useL & useR)
useFilt=find(useTime)

rownum=6; colnum=7
% rownum=round(sqrt(length(useFilt)+4))
% colnum=round(sqrt(length(useFilt)));
%%
figure('units','normalized','outerposition',[0 0 1 1])
for vid = 1:length(useFilt)
    subplot(rownum,colnum,vid)  ;
    bar([mean(isnan(mouse_xy{useFilt(vid),1}(1,:))) mean(isnan(cricket_xy{useFilt(vid),1}(1,:))) mean(isnan(cricketB{useFilt(vid),1}(1,:)))])
    ylabel('% error'); xlim([0.5 3.5]); ylim([0 1])
    set(gca,'XTick',[1 2 3])
    set(gca,'XTickLabel',{'mouse','c head', 'c body'})
end
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%
figure
for vid=1:length(useFilt)
    clear approach speed dec heading nframe d2cr dDist
   % subplot(rownum,colnum,vid)
  % subplot(7,6,vid)

    d2cr=(range{useFilt(vid)});
    dDist=diff(d2cr);
%  dDist=interpNan(dDist,3,'linear')
%  dDist=medfilt2(dDist)
    distThresh= d2cr<40; %long approaches can happen at up to 60cm away from cricket, which is diagonal through entire arena
    dec = (dDist<-.2) %threshold for change in range to cricket
    speed= mouseV{useFilt(vid)}>=5
    az=azDeg{useFilt(vid)}; %too many missing points - heading looks weird
    %azFilt=interpNan(az,3,'linear');
    heading=(az<90) &(az>-90) ;
    nframe=min(length(speed),length(heading));
    nframe=min(nframe,length(dec));
    speed=speed(1:nframe)';heading=heading(1:nframe)';dec=dec(1:nframe)'; distThresh=distThresh(1:nframe)';
    approach =dec==1&(speed==1)%&(distThresh==1)%&heading==1)
%     plot(d2cr,'b'); hold on;
%     plot(find(approach),d2cr(approach),'og');
    appEpoch{vid,:}=(approach); axis square
   % filtApp(vid,:)=sum(appEpoch{vid})
%   plot(d2cr(mouseV{useFilt(vid)}>5),'g');
end



%%
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    clear dur
    subplot(rownum,colnum,vid)
    plot(mouse_xy{useFilt(vid),1}(1,:),mouse_xy{useFilt(vid),1}(2,:)); hold on;
    plot(cricket_xy{useFilt(vid),1}(1,:),cricket_xy{useFilt(vid),1}(2,:))
    dur=num2str(length(mouse_xy{useFilt(vid),1}(1,:))/30,'%.2f');
    title([dur,'sec']); 
    xlim([300 1600]);
    use=appEpoch{vid};
    plot(mouse_xy{useFilt(vid),1}(1,use),mouse_xy{useFilt(vid),1}(2,use),'g'); hold on;
    end

if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%


clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(dTheta{useFilt(vid)}),length(dphiR{useFilt(vid)}));
    nframe = min(nframe, length(dphiL{useFilt(vid)}));
    dP=dTheta{useFilt(vid)}(1:nframe); dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);

    use = ~isnan(dP(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(dP(use),dpR(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(dP(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    use=appEpoch{vid}
      [corrRA lagsRA]= xcorr(dP(use),dpR(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
    [corrLA lagsLA]= xcorr(dP(use),dpL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    
    
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &&sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
          corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

    else
    end
end
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dPhi R','dPhi L','dPhi R Approach','dPhi L Approach'}); title('head phi, both eyes');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%% 
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
%     nframe = min(length(dTheta{useFilt(vid)}),length(dthetaR{useFilt(vid)}));
%     nframe = min(nframe, length(dthetaL{useFilt(vid)}));
    nonapp=appEpoch{useFilt(vid)}==0
  
    dP=dTheta{useFilt(vid)}; dpR=dthetaR{useFilt(vid)}; dpL=dthetaL{useFilt(vid)};
     clear use
    use = nonapp==1 & ~isnan(dP(1:length(dpR)));
    if sum(use)>3
    [corrR lagsR]= xcorr(dP(use),dpR(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(dP(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1
      [corrRA lagsRA]= xcorr(dP(use),dpR(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
    [corrLA lagsLA]= xcorr(dP(use),dpL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    
    
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        
    else
    end
end
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('head theta, both eyes');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(dTheta{useFilt(vid)}),length(dthetaR{useFilt(vid)}));
    nframe = min(nframe, length(dthetaL{useFilt(vid)}));
    dP=dTheta{useFilt(vid)}(1:nframe); dpR=dthetaR{useFilt(vid)}(1:nframe); dpL=dthetaL{useFilt(vid)}(1:nframe);

    use = ~isnan(dP(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(dP(use),dpR(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(dP(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    use=appEpoch{vid}
      [corrRA lagsRA]= xcorr(dP(use),dpR(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
    [corrLA lagsLA]= xcorr(dP(use),dpL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    
    
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        
    else
    end
end
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dtheta R','dtheta L','dtheta R Approach','dtheta L Approach'}); title('head theta, both eyes');


if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
corrRAAll=[]
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(dthetaR{useFilt(vid)}),length(dthetaL{useFilt(vid)}));
    dtR=dthetaR{useFilt(vid)}(1:nframe); dtL=dthetaL{useFilt(vid)}(1:nframe);
    nonapp=appEpoch{useFilt(vid)}==0;
    use = nonapp==1& ~isnan(dtR(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(dtR(use),dtL(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    clear nframe
    nframe = min(length(dphiR{useFilt(vid)}),length(dphiL{useFilt(vid)}));
    dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
    [corrL lagsL]= xcorr(dpR(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
      % use = ~isnan(dT(1:nframe));
       use=appEpoch{vid}==1;
%     if sum(use)>3
    [corrRA lagsRA]= xcorr(dtR(use),dtL(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
%     clear nframe
     nframe = min(length(dphiR{useFilt(vid)}),length(dphiL{useFilt(vid)}));
    dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
    [corrLA lagsLA]= xcorr(dpR(use),dpL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
   
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
       corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

    else
    end
end
title('between eye corr, dtheta & dphi - non app & approach');
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
      
      figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dTheta non-app','dPhi non-app','dTh Approach','dPhi app'}); title('mean between eye corr, dtheta & dphi');


if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end



%%
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
corrRAAll=[]
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(dthetaR{useFilt(vid)}),length(dthetaL{useFilt(vid)}));
    dtR=dthetaR{useFilt(vid)}(1:nframe); dtL=dthetaL{useFilt(vid)}(1:nframe);

    use = ~isnan(dtR(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(dtR(use),dtL(use),'coeff');
    plot(lagsR/30,corrR,'m');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    clear nframe
    nframe = min(length(dphiR{useFilt(vid)}),length(dphiL{useFilt(vid)}));
    dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
    [corrL lagsL]= xcorr(dpR(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'c');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
       % use = ~isnan(dT(1:nframe));
       use=appEpoch{vid};
%     if sum(use)>3
    [corrRA lagsRA]= xcorr(dtR(use),dtL(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
%     clear nframe
     nframe = min(length(dphiR{useFilt(vid)}),length(dphiL{useFilt(vid)}));
    dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
    [corrLA lagsLA]= xcorr(dpR(use),dpL(use),'coeff');
    plot(lagsLA/30,corrLA,'b');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
   
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
       corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

    else
    end
end
title('between eye corr, dtheta & dphi');
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
      
%%
% corrRAll = cell2mat(corrRAll)'; corrLAll = cell2mat(corrLAll)';
% corrRAAll = cell2mat(corrRAAll)'; corrLAAll = cell2mat(corrLAAll)';

figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-m',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-c',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-b',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'm-');
L(2) = plot(nan, nan, 'c-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'b-');

legend(L,{'dTheta','dPhi','dTh Approach','dPhi app'}); 
title('mean between eye corr, dtheta & dphi');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%


clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
corrRAAll=[]
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(dthetaR{useFilt(vid)}),length(dthetaL{useFilt(vid)}));
    dtR=dthetaR{useFilt(vid)}(1:nframe); dtL=dthetaL{useFilt(vid)}(1:nframe);
    nonapp=appEpoch{useFilt(vid)}==0;
    use = nonapp==1& ~isnan(dtR(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(dtR(use),dtL(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    clear nframe
    nframe = min(length(dphiR{useFilt(vid)}),length(dphiL{useFilt(vid)}));
    dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
    [corrL lagsL]= xcorr(dpR(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
      % use = ~isnan(dT(1:nframe));
       use=appEpoch{vid}==1;
%     if sum(use)>3
    [corrRA lagsRA]= xcorr(dtR(use),dtL(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
%     clear nframe
     nframe = min(length(dphiR{useFilt(vid)}),length(dphiL{useFilt(vid)}));
    dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
    [corrLA lagsLA]= xcorr(dpR(use),dpL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
   
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
       corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

    else
    end
end
title('between eye corr, dtheta & dphi - non app & approach');
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
      
      figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dTheta non-app','dPhi non-app','dTh Approach','dPhi app'}); title('mean between eye corr, dtheta & dphi');

%%


%%ABSOLUTE POSITION

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(theta{useFilt(vid)}),length(phiR{useFilt(vid)}));
    nframe = min(nframe, length(phiL{useFilt(vid)}));
    dP=theta{useFilt(vid)}(1:nframe); dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);

    use = ~isnan(dP(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(dP(use),dpR(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(dP(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    use=appEpoch{vid}
      [corrRA lagsRA]= xcorr(dP(use),dpR(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
    [corrLA lagsLA]= xcorr(dP(use),dpL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    
    
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &&sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
          corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

    else
    end
end
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dPhi R','dPhi L','dPhi R Approach','dPhi L Approach'}); title('head phi, both eyes');


if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%% 
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
%     nframe = min(length(theta{useFilt(vid)}),length(thetaR{useFilt(vid)}));
%     nframe = min(nframe, length(thetaL{useFilt(vid)}));
    nonapp=appEpoch{useFilt(vid)}==0
  
    dP=theta{useFilt(vid)}; dpR=thetaR{useFilt(vid)}; dpL=thetaL{useFilt(vid)};
     clear use
    use = nonapp==1 & ~isnan(dP(1:end-1))';
    if sum(use)>3
    [corrR lagsR]= xcorr(dP(use),dpR(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(dP(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1 & ~isnan(dP(1:end-1))';
      [corrRA lagsRA]= xcorr(dP(use),dpR(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
    [corrLA lagsLA]= xcorr(dP(use),dpL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    
    
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        
    else
    end
end
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'}); title('head theta, both eyes');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(theta{useFilt(vid)}),length(thetaR{useFilt(vid)}));
    nframe = min(nframe, length(thetaL{useFilt(vid)}));
    dP=theta{useFilt(vid)}(1:nframe); dpR=thetaR{useFilt(vid)}(1:nframe); dpL=thetaL{useFilt(vid)}(1:nframe);

    use = ~isnan(dP(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(dP(use),dpR(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(dP(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    use=appEpoch{vid}
      [corrRA lagsRA]= xcorr(dP(use),dpR(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
    [corrLA lagsLA]= xcorr(dP(use),dpL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    
    
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        
    else
    end
end
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'theta R','theta L','theta R Approach','theta L Approach'}); title('head theta, both eyes');


if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
corrRAAll=[]
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(thetaR{useFilt(vid)}),length(thetaL{useFilt(vid)}))-1;
    dtR=thetaR{useFilt(vid)}(1:nframe); dtL=thetaL{useFilt(vid)}(1:nframe);
    nonapp=appEpoch{useFilt(vid)}(1:nframe)==0;
    
    use = nonapp==1 & ~isnan(dtR(1:length(nonapp)));
    if sum(use)>3
    [corrR lagsR]= xcorr(dtR(use),dtL(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    clear nframe
    nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
    dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
    [corrL lagsL]= xcorr(dpR(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
      % use = ~isnan(dT(1:nframe));
       use=appEpoch{vid}==1;
%     if sum(use)>3
    [corrRA lagsRA]= xcorr(dtR(use),dtL(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
%     clear nframe
     nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
    dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
    [corrLA lagsLA]= xcorr(dpR(use),dpL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
   
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
       corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

    else
    end
end
title('between eye corr, theta & dphi - non app & approach');
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
      
      figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'theta non-app','dPhi non-app','dTh Approach','dPhi app'}); title('mean between eye corr, theta & dphi');


if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end



%%
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
corrRAAll=[]
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(thetaR{useFilt(vid)}),length(thetaL{useFilt(vid)}));
    dtR=thetaR{useFilt(vid)}(1:nframe); dtL=thetaL{useFilt(vid)}(1:nframe);

    use = ~isnan(dtR(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(dtR(use),dtL(use),'coeff');
    plot(lagsR/30,corrR,'m');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    clear nframe
    nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
    dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
    [corrL lagsL]= xcorr(dpR(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'c');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
       % use = ~isnan(dT(1:nframe));
       use=appEpoch{vid};
%     if sum(use)>3
    [corrRA lagsRA]= xcorr(dtR(use),dtL(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
%     clear nframe
     nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
    dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
    [corrLA lagsLA]= xcorr(dpR(use),dpL(use),'coeff');
    plot(lagsLA/30,corrLA,'b');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
   
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
       corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

    else
    end
end
title('between eye corr, theta & dphi');
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
      
%%
% corrRAll = cell2mat(corrRAll)'; corrLAll = cell2mat(corrLAll)';
% corrRAAll = cell2mat(corrRAAll)'; corrLAAll = cell2mat(corrLAAll)';

figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-m',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-c',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-b',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'm-');
L(2) = plot(nan, nan, 'c-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'b-');

legend(L,{'theta','dPhi','dTh Approach','dPhi app'}); 
title('mean between eye corr, theta & dphi');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%


clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
corrRAAll=[]
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(thetaR{useFilt(vid)}),length(thetaL{useFilt(vid)}))-1;
    dtR=thetaR{useFilt(vid)}(1:nframe); dtL=thetaL{useFilt(vid)}(1:nframe);
    nonapp=appEpoch{useFilt(vid)}==0;
    use = nonapp==1& ~isnan(dtR(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(dtR(use),dtL(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    clear nframe
    nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
    dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
    [corrL lagsL]= xcorr(dpR(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
      % use = ~isnan(dT(1:nframe));
       use=appEpoch{vid}==1;
%     if sum(use)>3
    [corrRA lagsRA]= xcorr(dtR(use),dtL(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
%     clear nframe
     nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
    dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
    [corrLA lagsLA]= xcorr(dpR(use),dpL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
   
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
       corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

    else
    end
end
title('between eye corr, theta & dphi - non app & approach');
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
      
      figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'theta non-app','dPhi non-app','dTh Approach','dPhi app'}); title('mean between eye corr, theta & dphi');

      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%
%SCATTER PLOTS
clear dTAll dTRALL dTLAll dPRAll dPLAll
close all

dTRAll=[];

% figure%('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    clear dT dpR dpL dtR dtL
%     subplot(rownum,colnum,vid);
    nframe = min(length(dTheta{useFilt(vid)}),length(dphiR{useFilt(vid)}));
    nframe = min(nframe, length(dphiL{useFilt(vid)}));
    nframe = min(nframe, length(dthetaR{useFilt(vid)}));
    nframe = min(nframe,length(dthetaL{useFilt(vid)}));
    dT=dTheta{useFilt(vid)}(1:nframe); dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
    dtR=dthetaR{useFilt(vid)}(1:nframe); dtL=dthetaL{useFilt(vid)}(1:nframe);
    figure(1);
    plot(dT(1:15:end),dpR(1:15:end),'b.');axis square; hold on
    title('dHead theta, dPhi R');
        x = linspace(-40,40);
    y = linspace(-60,60);
    plot(-x,y); xlim([-40 40]); ylim([-60 60]);
        use = (appEpoch{vid});
plot(dT(use(1:15:end)),dpR(use(1:15:end)),'.g');
      if vid==(useFilt(end))        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(2);
    plot(dT(1:15:end),dpL(1:15:end),'b.'); axis square; hold on
    title('dHead theta, dPhi L');
    x = linspace(-40,40);
   y = linspace(-80,80);
   plot(x,y);  xlim([-40 40]); ylim([-80 80]);
    plot(dT(use(1:15:end)),dpL(use(1:15:end)),'.g');
  if vid==(useFilt(end))       
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(3);
    plot(dT(1:15:end),dtR(1:15:end),'b.');axis square; hold on
   title('dHead theta, dtheta R');
       x = linspace(-40,40);
    y = linspace(-40,40);  xlim([-40 40]); ylim([-40 40]);
    plot(-x,y);
    plot(dT(use(1:15:end)),dtR(use(1:15:end)),'.g');
  if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(4);
    plot(dT(1:15:end),dtL(1:15:end),'b.');axis square; hold on;
    title('dHead theta, dtheta L');
       x = linspace(-40,40);  xlim([-40 40]); ylim([-80 80]);
    y = linspace(-80,80);
    plot(x,y);
        plot(dT(use(1:15:end)),dtL(use(1:15:end)),'.g');
  if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
        figure(5);%subplot(rownum,colnum,vid);
    plot(dtR(1:15:end),dtL(1:15:end),'b.');axis square; hold on;
    title('r theta, l theta');
    x = linspace(-80,80);
    y = linspace(-50,50);
    plot(-x,y); xlim([-80 80]); ylim([-50 50]);
            plot(dtR(use(1:15:end)),dtL(use(1:15:end)),'.g');
  if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
        figure(6);
    plot(dpR(1:15:end),dpL(1:15:end),'b.');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('right phi, left phi');
                plot(dpR(use(1:15:end)),dpL(use(1:15:end)),'.g');
  if vid==(useFilt(end) )       
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end

end

%%
%app and non app
clear dTAll dTRALL dTLAll dPRAll dPLAll
close all
dTRAll=[];

% figure%('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    clear dT dpR dpL dtR dtL
    %     subplot(rownum,colnum,vid);
    nframe = min(length(dTheta{useFilt(vid)}),length(dphiR{useFilt(vid)}));
    nframe = min(nframe, length(dphiL{useFilt(vid)}));
    nframe = min(nframe, length(dthetaR{useFilt(vid)}));
    nframe = min(nframe,length(dthetaL{useFilt(vid)}));
    dT=dTheta{useFilt(vid)}(1:nframe); dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
    dtR=dthetaR{useFilt(vid)}(1:nframe); dtL=dthetaL{useFilt(vid)}(1:nframe);
    useN= appEpoch{vid}==0;
    use = (appEpoch{vid});

    
    figure(1);
    plot(dT(useN(1:15:end)),dpR(useN(1:15:end)),'b.');axis square; hold on
    title('dHead theta, dPhi R, nonapp & app');
    x = linspace(-40,40);
    y = linspace(-60,60);
    plot(-x,y); xlim([-40 40]); ylim([-60 60]);
    plot(dT(use(1:15:end)),dpR(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(2);
    plot(dT(useN(1:15:end)),dpL(useN(1:15:end)),'b.'); axis square; hold on
    title('dHead theta, dPhi L, nonapp & app');
    x = linspace(-40,40);
    y = linspace(-80,80);
    plot(x,y);  xlim([-40 40]); ylim([-80 80]);
    plot(dT(use(1:15:end)),dpL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(3);
    plot(dT(useN(1:15:end)),dtR(useN(1:15:end)),'b.');axis square; hold on
    title('dHead theta, dtheta R, nonapp & app');
    x = linspace(-40,40);
    y = linspace(-40,40);  xlim([-40 40]); ylim([-40 40]);
    plot(-x,y);
    plot(dT(use(1:15:end)),dtR(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(4);
    plot(dT(useN(1:15:end)),dtL(useN(1:15:end)),'b.');axis square; hold on;
    title('dHead theta, dtheta L, nonapp & app');
    x = linspace(-40,40);  xlim([-40 40]); ylim([-80 80]);
    y = linspace(-80,80);
    plot(x,y);
    plot(dT(use(1:15:end)),dtL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(5);%subplot(rownum,colnum,vid);
    plot(dtR(useN(1:15:end)),dtL(useN(1:15:end)),'b.');axis square; hold on;
    title('r theta, l theta, nonapp & app');
    x = linspace(-80,80);
    y = linspace(-50,50);
    plot(-x,y); xlim([-80 80]); ylim([-50 50]);
    plot(dtR(use(1:15:end)),dtL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(6);
    plot(dpR(useN(1:15:end)),dpL(useN(1:15:end)),'b.');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('right phi, left phi, nonapp & app');
    plot(dpR(use(1:15:end)),dpL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
end

%%
%ABSOLUTE POSITON

clear dTAll dTRALL dTLAll dPRAll dPLAll
close all
dTRAll=[];

% figure%('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    clear dT dpR dpL dtR dtL
%     subplot(rownum,colnum,vid);
    nframe = min(length(theta{useFilt(vid)}),length(phiR{useFilt(vid)}));
    nframe = min(nframe, length(phiL{useFilt(vid)}));
    nframe = min(nframe, length(thetaR{useFilt(vid)}));
    nframe = min(nframe,length(thetaL{useFilt(vid)}));
    dT=theta{useFilt(vid)}(1:nframe); dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
    dtR=thetaR{useFilt(vid)}(1:nframe); dtL=thetaL{useFilt(vid)}(1:nframe);
    figure(1);
    plot(dT(1:15:end),dpR(1:15:end),'b.');axis square; hold on
    title('head theta, phi R');
        x = linspace(-40,40);
    y = linspace(-60,60);
    plot(-x,y); %xlim([-40 40]); ylim([-60 60]);
        use = (appEpoch{vid});
plot(dT(use(1:15:end)),dpR(use(1:15:end)),'.g');
      if vid==(useFilt(end))        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(2);
    plot(dT(1:15:end),dpL(1:15:end),'b.'); axis square; hold on
    title('head theta, phi L');
    x = linspace(-40,40);
   y = linspace(-80,80);
   plot(x,y);  %xlim([-40 40]); ylim([-80 80]);
    plot(dT(use(1:15:end)),dpL(use(1:15:end)),'.g');
  if vid==(useFilt(end))       
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(3);
    plot(dT(1:15:end),dtR(1:15:end),'b.');axis square; hold on
   title('head theta, theta R');
       x = linspace(-40,40);
    y = linspace(-40,40); % xlim([-40 40]); ylim([-40 40]);
    plot(-x,y);
    plot(dT(use(1:15:end)),dtR(use(1:15:end)),'.g');
  if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(4);
    plot(dT(1:15:end),dtL(1:15:end),'b.');axis square; hold on;
    title('Head theta, theta L');
%        x = linspace(-40,40); % xlim([-40 40]); ylim([-80 80]);
%     y = linspace(-80,80);
    plot(x,y);
        plot(dT(use(1:15:end)),dtL(use(1:15:end)),'.g');
  if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
        figure(5);%subplot(rownum,colnum,vid);
    plot(dtR(1:15:end),dtL(1:15:end),'b.');axis square; hold on;
    title('r theta, l theta');
    x = linspace(-80,80);
    y = linspace(-50,50);
    plot(-x,y); xlim([-80 80]); ylim([-50 50]);
            plot(dtR(use(1:15:end)),dtL(use(1:15:end)),'.g');
  if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
        figure(6);
    plot(dpR(1:15:end),dpL(1:15:end),'b.');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('right phi, left phi');
                plot(dpR(use(1:15:end)),dpL(use(1:15:end)),'.g');
  if vid==(useFilt(end))       
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end

end

%%
%app and non app, abs position
clear dTAll dTRALL dTLAll dPRAll dPLAll

dTRAll=[];
close all
% figure%('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    clear dT dpR dpL dtR dtL
    %     subplot(rownum,colnum,vid);
    nframe = min(length(theta{useFilt(vid)}),length(phiR{useFilt(vid)}));
    nframe = min(nframe, length(phiL{useFilt(vid)}));
    nframe = min(nframe, length(thetaR{useFilt(vid)}));
    nframe = min(nframe,length(thetaL{useFilt(vid)}));
    dT=theta{useFilt(vid)}(1:nframe); dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
    dtR=thetaR{useFilt(vid)}(1:nframe); dtL=thetaL{useFilt(vid)}(1:nframe);
    useN= appEpoch{vid}==0;
    use = (appEpoch{vid});

    
    figure(1);
    plot(dT(useN(1:15:end)),dpR(useN(1:15:end)),'b.');axis square; hold on
    title('head theta, phi R, nonapp & app');
    x = linspace(-40,40);
    y = linspace(-60,60);
    plot(-x,y); xlim([-40 40]); ylim([-60 60]);
    plot(dT(use(1:15:end)),dpR(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(2);
    plot(dT(useN(1:15:end)),dpL(useN(1:15:end)),'b.'); axis square; hold on
    title('head theta, phi L, nonapp & app');
    x = linspace(-40,40);
    y = linspace(-80,80);
    plot(x,y);  xlim([-40 40]); ylim([-80 80]);
    plot(dT(use(1:15:end)),dpL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(3);
    plot(dT(useN(1:15:end)),dtR(useN(1:15:end)),'b.');axis square; hold on
    title('head theta, theta R, nonapp & app');
    x = linspace(-40,40);
    y = linspace(-40,40);  xlim([-40 40]); ylim([-40 40]);
    plot(-x,y);
    plot(dT(use(1:15:end)),dtR(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(4);
    plot(dT(useN(1:15:end)),dtL(useN(1:15:end)),'b.');axis square; hold on;
    title('head theta, theta L, nonapp & app');
    x = linspace(-40,40);  xlim([-40 40]); ylim([-80 80]);
    y = linspace(-80,80);
    plot(x,y);
    plot(dT(use(1:15:end)),dtL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(5);%subplot(rownum,colnum,vid);
    plot(dtR(useN(1:15:end)),dtL(useN(1:15:end)),'b.');axis square; hold on;
    title('r theta, l theta, nonapp & app');
    x = linspace(-80,80);
    y = linspace(-50,50);
    plot(-x,y); xlim([-80 80]); ylim([-50 50]);
    plot(dtR(use(1:15:end)),dtL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(6);
    plot(dpR(useN(1:15:end)),dpL(useN(1:15:end)),'b.');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('right phi, left phi, nonapp & app');
    plot(dpR(use(1:15:end)),dpL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
end

%%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     subplot(rownum,colnum,vid);
%     nframe = min(length(theta{useFilt(vid)}),length(YRcent{useFilt(vid)}));
%     nframe = min(nframe, length(YLcent{useFilt(vid)}));
%     dP=theta{useFilt(vid)}(1:nframe); dpR=YRcent{useFilt(vid)}(1:nframe); dpL=YLcent{useFilt(vid)}(1:nframe);
% 
%     use = ~isnan(dP(1:nframe));
%     if sum(use)>3
%     [corrR lagsR]= xcorr(dP(use),dpR(use),'coeff');
%     plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
%     uselagsR=(lagsR>=-30& lagsR<=30);
%     
%     [corrL lagsL]= xcorr(dP(use),dpL(use),'coeff');
%     plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
%     uselagsL=(lagsL>=-30 & lagsL<=30);
%     else
%     end
%     use=appEpoch{vid}
%       [corrRA lagsRA]= xcorr(dP(use),dpR(use),'coeff');
%     plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
%     uselagsRA=(lagsRA>=-30& lagsRA<=30);
%     
%     [corrLA lagsLA]= xcorr(dP(use),dpL(use),'coeff');
%     plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
%     uselagsLA=(lagsLA>=-30 & lagsLA<=30);
%     
%     
%     
%     if sum(uselagsR)==61 & sum(uselagsL)==61 &&sum(uselagsRA)==61 & sum(uselagsLA)==61
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%           corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
% 
%     else
%     end
% end
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
% plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
% 
% legend(L,{'dPhi R','dPhi L','dPhi R Approach','dPhi L Approach'}); title('head phi, both eyes');
% 
% 
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% 
% %%
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     subplot(rownum,colnum,vid);
% %     nframe = min(length(theta{useFilt(vid)}),length(XRcent{useFilt(vid)}));
% %     nframe = min(nframe, length(XLcent{useFilt(vid)}));
%     nonapp=appEpoch{useFilt(vid)}==0
%   
%     dP=theta{useFilt(vid)}; dpR=XRcent{useFilt(vid)}; dpL=XLcent{useFilt(vid)};
%      clear use
%     use = nonapp==1 & ~isnan(dP(1:end-1))';
%     if sum(use)>3
%     [corrR lagsR]= xcorr(dP(use),dpR(use),'coeff');
%     plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
%     uselagsR=(lagsR>=-30& lagsR<=30);
%     
%     [corrL lagsL]= xcorr(dP(use),dpL(use),'coeff');
%     plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
%     uselagsL=(lagsL>=-30 & lagsL<=30);
%     else
%     end
%     clear use;
%     use=appEpoch{vid}==1 & ~isnan(dP(1:end-1))';
%       [corrRA lagsRA]= xcorr(dP(use),dpR(use),'coeff');
%     plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
%     uselagsRA=(lagsRA>=-30& lagsRA<=30);
%     
%     [corrLA lagsLA]= xcorr(dP(use),dpL(use),'coeff');
%     plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
%     uselagsLA=(lagsLA>=-30 & lagsLA<=30);
%     
%     
%     
%     if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%         
%     else
%     end
% end
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
% plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
% 
% legend(L,{'Xcent R non-app','X cent L non-app','X cent R Approach','X cent L Approach'}); title('head theta, both eyes x position');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% 
% %%
% 
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     subplot(rownum,colnum,vid);
%     nframe = min(length(theta{useFilt(vid)}),length(XRcent{useFilt(vid)}));
%     nframe = min(nframe, length(XLcent{useFilt(vid)}));
%     dP=theta{useFilt(vid)}(1:nframe); dpR=XRcent{useFilt(vid)}(1:nframe); dpL=XLcent{useFilt(vid)}(1:nframe);
% 
%     use = ~isnan(dP(1:nframe));
%     if sum(use)>3
%     [corrR lagsR]= xcorr(dP(use),dpR(use),'coeff');
%     plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
%     uselagsR=(lagsR>=-30& lagsR<=30);
%     
%     [corrL lagsL]= xcorr(dP(use),dpL(use),'coeff');
%     plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
%     uselagsL=(lagsL>=-30 & lagsL<=30);
%     else
%     end
%     use=appEpoch{vid}
%       [corrRA lagsRA]= xcorr(dP(use),dpR(use),'coeff');
%     plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
%     uselagsRA=(lagsRA>=-30& lagsRA<=30);
%     
%     [corrLA lagsLA]= xcorr(dP(use),dpL(use),'coeff');
%     plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
%     uselagsLA=(lagsLA>=-30 & lagsLA<=30);
%     
%     
%     
%     if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%         
%     else
%     end
% end
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
% plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
% 
% legend(L,{'xcent R','xcent L','xcent R Approach','xcent L Approach'}); title('head theta & X pos both eyes');
% 
% 
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% 
% %%
% 
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% corrRAAll=[]
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     subplot(rownum,colnum,vid);
%     nframe = min(length(XRcent{useFilt(vid)}),length(XLcent{useFilt(vid)}))-1;
%     dtR=XRcent{useFilt(vid)}(1:nframe); dtL=XLcent{useFilt(vid)}(1:nframe);
%     nonapp=appEpoch{useFilt(vid)}(1:nframe)==0;
%     
%     use = nonapp==1 & ~isnan(dtR(1:length(nonapp)));
%     if sum(use)>3
%     [corrR lagsR]= xcorr(dtR(use),dtL(use),'coeff');
%     plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
%     uselagsR=(lagsR>=-30& lagsR<=30);
%     
%     clear nframe
%     nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
%     dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
%     [corrL lagsL]= xcorr(dpR(use),dpL(use),'coeff');
%     plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
%     uselagsL=(lagsL>=-30 & lagsL<=30);
%     else
%     end
%       % use = ~isnan(dT(1:nframe));
%        use=appEpoch{vid}==1;
% %     if sum(use)>3
%     [corrRA lagsRA]= xcorr(dtR(use),dtL(use),'coeff');
%     plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
%     uselagsRA=(lagsRA>=-30& lagsRA<=30);
%     
% %     clear nframe
%      nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
%     dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
%     [corrLA lagsLA]= xcorr(dpR(use),dpL(use),'coeff');
%     plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
%     uselagsLA=(lagsLA>=-30 & lagsLA<=30);
%    
%     
%     if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
% 
%     else
%     end
% end
% title('between eye corr, x & y position - non app & approach');
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
% plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
% 
% legend(L,{'X cent non-app','y Cent non-app','X-cent Approach','y-cent app'}); title('mean between eye corr, x cent & y cent');
% 
% 
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% 
% 
% 
% %%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% corrRAAll=[]
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     subplot(rownum,colnum,vid);
%     nframe = min(length(XRcent{useFilt(vid)}),length(XLcent{useFilt(vid)}));
%     dtR=XRcent{useFilt(vid)}(1:nframe); dtL=XLcent{useFilt(vid)}(1:nframe);
% 
%     use = ~isnan(dtR(1:nframe));
%     if sum(use)>3
%     [corrR lagsR]= xcorr(dtR(use),dtL(use),'coeff');
%     plot(lagsR/30,corrR,'m');xlim([-.3 .3]);hold on;
%     uselagsR=(lagsR>=-30& lagsR<=30);
%     
%     clear nframe
%     nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
%     dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
%     [corrL lagsL]= xcorr(dpR(use),dpL(use),'coeff');
%     plot(lagsL/30,corrL,'c');xlim([-.3 .3]);
%     uselagsL=(lagsL>=-30 & lagsL<=30);
%     else
%     end
%        % use = ~isnan(dT(1:nframe));
%        use=appEpoch{vid};
% %     if sum(use)>3
%     [corrRA lagsRA]= xcorr(dtR(use),dtL(use),'coeff');
%     plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
%     uselagsRA=(lagsRA>=-30& lagsRA<=30);
%     
% %     clear nframe
%      nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
%     dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
%     [corrLA lagsLA]= xcorr(dpR(use),dpL(use),'coeff');
%     plot(lagsLA/30,corrLA,'b');xlim([-.3 .3]);
%     uselagsLA=(lagsLA>=-30 & lagsLA<=30);
%    
%     if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
% 
%     else
%     end
% end
% title('between eye corr, x cent and y cent');
%       if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%       
% %%
% % corrRAll = cell2mat(corrRAll)'; corrLAll = cell2mat(corrLAll)';
% % corrRAAll = cell2mat(corrRAAll)'; corrLAAll = cell2mat(corrLAAll)';
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
% 
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-m',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-c',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-b',1);
% 
% plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'm-');
% L(2) = plot(nan, nan, 'c-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'b-');
% 
% legend(L,{'X cent','Y cent','X cent Approach','Y cent app'}); 
% title('mean between eye corr, theta & dphi');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %%
% 
% 
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% corrRAAll=[]
% figure('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useFilt)
%     subplot(rownum,colnum,vid);
%     nframe = min(length(XRcent{useFilt(vid)}),length(XLcent{useFilt(vid)}))-1;
%     dtR=XRcent{useFilt(vid)}(1:nframe); dtL=XLcent{useFilt(vid)}(1:nframe);
%     nonapp=appEpoch{useFilt(vid)}==0;
%     use = nonapp==1& ~isnan(dtR(1:nframe));
%     if sum(use)>3
%     [corrR lagsR]= xcorr(dtR(use),dtL(use),'coeff');
%     plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
%     uselagsR=(lagsR>=-30& lagsR<=30);
%     
%     clear nframe
%     nframe = min(length(YRcent{useFilt(vid)}),length(YLcent{useFilt(vid)}));
%     dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
%     [corrL lagsL]= xcorr(dpR(use),dpL(use),'coeff');
%     plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
%     uselagsL=(lagsL>=-30 & lagsL<=30);
%     else
%     end
%       % use = ~isnan(dT(1:nframe));
%        use=appEpoch{vid}==1;
% %     if sum(use)>3
%     [corrRA lagsRA]= xcorr(dtR(use),dtL(use),'coeff');
%     plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
%     uselagsRA=(lagsRA>=-30& lagsRA<=30);
%     
% %     clear nframe
%      nframe = min(length(YRcent{useFilt(vid)}),length(YLcent{useFilt(vid)}));
%     dpR=phiR{useFilt(vid)}(1:nframe); dpL=phiL{useFilt(vid)}(1:nframe);
%     [corrLA lagsLA]= xcorr(dpR(use),dpL(use),'coeff');
%     plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
%     uselagsLA=(lagsLA>=-30 & lagsLA<=30);
%    
%     
%     if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
% 
%     else
%     end
% end
% title('between eye corr, X & Y cent - non app & approach');
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
% plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
% 
% legend(L,{'X cent non-app','Y cent non-app','X-cent Approach','Y-cent app'}); title('mean between eye corr, x & y cent');
% 
%       if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% 
% 
% 

if savePDF
    pSname='T:\PreyCaptureAnalysis\Data\';
    filen=sprintf('%s',ani,'AnalyzedTS_AllVideos_withApproach_test082219a','.pdf')
    pdfilename=fullfile(pSname,filen);
    dos(['ps2pdf ' psfilename ' ' pdfilename]);
    delete(psfilename);
end

afilename=sprintf('%s',ani,'Analyzed_Test082219a','.mat')
save(fullfile(pSname, afilename))