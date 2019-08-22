clear all; close all
load('J463cAllVids.mat'); 
savePDF=1;
if savePDF
    psfilename = 'C:\analysisPS.ps';
    if exist(psfilename,'file')==2;delete(psfilename);end
end


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
    theta{i,:}= Data(i).theta;
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
    
    XRcent{i,:} =Data(i).dxR;
  %  YRcent{i,:} =Data(i).YRcent;
    XLcent{i,:} =Data(i).dxL;
  %  YLcent{i,:} =Data(i).YLcent;
 
    
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
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    clear dur
    subplot(rownum,colnum,vid)
    plot(mouse_xy{useFilt(vid),1}(1,:),mouse_xy{useFilt(vid),1}(2,:)); hold on;
    plot(cricket_xy{useFilt(vid),1}(1,:),cricket_xy{useFilt(vid),1}(2,:))
    dur=num2str(length(mouse_xy{useFilt(vid),1}(1,:))/30,'%.2f');
    title([dur,'sec']); 
    xlim([300 1600]);
  %  use=appEpoch{vid};
  %  plot(mouse_xy{useFilt(vid),1}(1,use),mouse_xy{useFilt(vid),1}(2,use),'og'); hold on;
    end

if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
figure
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid)
    plot(rangeRaw{useFilt(vid)});hold on
    plot(range{useFilt(vid)});hold on
    axis square; xlabel('frame number'); ylabel('range (cm)');
   ylim([0 60]); %%% max possible distance
end
legend('raw range','interp range');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%
figure
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid)
%     plot(normalize(theta{useFilt(vid)}));hold on
    plot(normalize(phiR{useFilt(vid)}));hold on
    plot(normalize(phiL{useFilt(vid)}));hold on
     plot(normalize(theta{useFilt(vid)}));hold on

    axis square; xlabel('frame number'); ylabel('norm');

end
legend('phiR','phiL','head theta');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%
figure
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid)
     plot(normalize(thetaRaw{useFilt(vid)}));hold on
   
     plot(normalize(theta{useFilt(vid)}));hold on
    axis square; xlabel('frame number'); ylabel('head theta');
    ylim([-pi pi]);

end
legend('raw theta', 'interp theta');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%
figure
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid)
    plot(normalize(thetaR{useFilt(vid)}));hold on
    plot(normalize(thetaL{useFilt(vid)}));hold on
    plot(normalize(theta{useFilt(vid)}));hold on
    axis square; xlabel('frame number'); ylabel('norm');
    
end
legend('thetaR','thetaL','head theta');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL

figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(dTheta{useFilt(vid)}),length(dthetaR{useFilt(vid)}));
    nframe = min(nframe, length(dthetaL{useFilt(vid)}));
    dT=dTheta{useFilt(vid)}(1:nframe); dtR=dthetaR{useFilt(vid)}(1:nframe); dtL=dthetaL{useFilt(vid)}(1:nframe);

    use = ~isnan(dT(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(dT(use),dtR(use),'coeff');
    plot(lagsR/30,corrR);xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(dT(use),dtL(use),'coeff');
    plot(lagsL/30,corrL);xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    if sum(uselagsR)==61 & sum(uselagsL)==61
        corrRAll{vid}=corrR(uselagsR); corrLAll{vid}=corrL(uselagsL);
    else
    end
end
%%
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

corrRAll = cell2mat(corrRAll)'; corrLAll = cell2mat(corrLAll)';

figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]);
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
legend(L,{'right eye','left eye'}); title('mean dTheta Head and dTheta eye corr');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL

figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(dTheta{useFilt(vid)}),length(dphiR{useFilt(vid)}));
    nframe = min(nframe, length(dphiL{useFilt(vid)}));
    dP=dTheta{useFilt(vid)}(1:nframe); dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);

    use = ~isnan(dP(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(dP(use),dpR(use),'coeff');
    plot(lagsR/30,corrR);xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(dP(use),dpL(use),'coeff');
    plot(lagsL/30,corrL);xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    if sum(uselagsR)==61 & sum(uselagsL)==61
        corrRAll{vid}=corrR(uselagsR); corrLAll{vid}=corrL(uselagsL);
    else
    end
end
%%
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

corrRAll = cell2mat(corrRAll)'; corrLAll = cell2mat(corrLAll)';

figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]);
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
legend(L,{'right eye','left eye'}); title('mean dTheta Head and dphi eye corr');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL

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
    if sum(uselagsR)==61 & sum(uselagsL)==61
        corrRAll{vid}=corrR(uselagsR); corrLAll{vid}=corrL(uselagsL);
    else
    end
end
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

corrRAll = cell2mat(corrRAll)'; corrLAll = cell2mat(corrLAll)';

figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-m',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-c',1);
plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'm-');
L(2) = plot(nan, nan, 'c-');
legend(L,{'dTheta','dPhi'}); title('mean between eye corr, dtheta & dphi');

if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%

clear dTAll dTRALL dTLAll dPRAll dPLAll

dTRAll=[];

figure%('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    clear dT dpR dpL dtR dtL
%     subplot(rownum,colnum,vid);
    nframe = min(length(dTheta{useFilt(vid)}),length(dphiR{useFilt(vid)}));
    nframe = min(nframe, length(dphiL{useFilt(vid)}));
    nframe = min(nframe, length(dthetaR{useFilt(vid)}));
    nframe = min(nframe,length(dthetaL{useFilt(vid)}));
    dT=dTheta{useFilt(vid)}(1:nframe); dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
    dtR=dthetaR{useFilt(vid)}(1:nframe); dtL=dthetaL{useFilt(vid)}(1:nframe);
%     figure(1);
%     subplot(rownum,colnum,vid);
%     plot(dT,dpR,'.');axis square
%     title('dHead theta, dPhi R');
%     figure(2);
%        subplot(rownum,colnum,vid);
%     plot(dT,dpL,'.'); axis square
%     title('dHead theta, dPhi L');
%     figure(3);
subplot(rownum,colnum,vid);
    plot(dT,dtR,'.');axis square
  %  title('dHead theta, dtheta R');
%     figure(4);subplot(rownum,colnum,vid);
%     plot(dT,dtL,'.');axis square
%     title('dHead theta, dtheta L');
% if length(dT)>500
% dTAll(vid,:)=datasample(dT,500);
% dTRAll(vid,:)=datasample(dtR,500);
% dTLAll(vid,:)=datasample(dtL,500);
% dPRAll(vid,:)=datasample(dpR,500);
% dPLAll(vid,:)=datasample(dpL,500);
% else
% end
end
    title('dHead theta, dtheta R');
   % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%

figure;plot(nanmean(dTAll),nanmean(dTRAll),'.'); axis square; hold on
plot(mean(dTAll), mean(dTLAll),'.');
x = linspace(-.0314,.0314);
y = linspace(-pi,pi);
plot(x,y);

xlabel('eye theta'); ylabel('eye theta');
title('mean headTheta & eye theta');
legend('right theta','left theta');
    if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


figure;plot(nanmean(dTAll),nanmean(dPRAll),'.'); axis square; hold on
plot(nanmean(dTAll), nanmean(dPLAll),'.');
x = linspace(-.0314,.0314);
y = linspace(-pi,pi);
plot(x,y);

xlabel('eye theta'); ylabel('eye phi');
title('mean headTheta & deye phi');
legend('right dphi','left dphi');
    if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%
figure;plot(nanmean(dTRAll),nanmean(dTLAll),'.'); axis square; hold on
y = linspace(-pi,pi);
plot(y,y);

xlabel('right eye theta'); ylabel('left eye theta');
title('mean eye theta');
    if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

figure;plot(nanmean(dPRAll),nanmean(dPLAll),'.'); axis square; hold on
y = linspace(-pi,pi);
plot(y,y);

xlabel('right eye phi'); ylabel('left eye phi');
title('mean eye phi');
    if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
clear dTAll dTRALL dTLAll dPRAll dPLAll

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
%     subplot(rownum,colnum,vid);
    plot(dT(1:15:end),dpR(1:15:end),'.');axis square; hold on
    title('dHead theta, dPhi R');
        x = linspace(-.6,.60);
    y = linspace(-80,80);
    plot(x,y); xlim([-.6 .6]); ylim([-80 80]);
    if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(2);
%        subplot(rownum,colnum,vid);
    plot(dT(1:15:end),dpL(1:15:end),'.'); axis square; hold on
    title('dHead theta, dPhi L');
        x = linspace(-.6,.60);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-.6 .6]); ylim([-80 80]);
      if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(3);
%subplot(rownum,colnum,vid);
    plot(dT(1:15:end),dtR(1:15:end),'.');axis square; hold on
   title('dHead theta, dtheta R');
       x = linspace(-.6,.60);
    y = linspace(-80,80);  xlim([-.6 .6]); ylim([-80 80]);
    plot(-x,y);
      if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(4);%subplot(rownum,colnum,vid);
    plot(dT(1:15:end),dtL(1:15:end),'.');axis square; hold on;
    title('dHead theta, dtheta L');
       x = linspace(-.6,.60);  xlim([-.6 .6]); ylim([-50 50]);
    y = linspace(-50,50);
    plot(x,y);
      if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
        figure(5);%subplot(rownum,colnum,vid);
    plot(dtR(1:15:end),dtL(1:15:end),'.');axis square; hold on;
    title('r theta, l theta');
    x = linspace(-80,80);
    y = linspace(-50,50);
    plot(-x,y); xlim([-80 80]); ylim([-50 50]);
      if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
        figure(6);%subplot(rownum,colnum,vid);
    plot(dpR(1:15:end),dpL(1:15:end),'.');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('right phi, left phi');
      if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
% if length(dT)>500
% dTAll(vid,:)=datasample(dT,500);
% dTRAll(vid,:)=datasample(dtR,500);
% dTLAll(vid,:)=datasample(dtL,500);
% dPRAll(vid,:)=datasample(dpR,500);
% dPLAll(vid,:)=datasample(dpL,500);
% else
% end
end
%     title('dHead theta, dtheta R');

%%



figure
for vid=1:length(useFilt)
    clear approach speed dec heading nframe d2cr dDist
    subplot(rownum,colnum,vid)

    d2cr=(range{useFilt(vid)});
    dDist=diff(d2cr);
%  dDist=interpNan(dDist,3,'linear')
%  dDist=medfilt2(dDist)
    distThresh= d2cr<40; %long approaches can happen at up to 60cm away from cricket, which is diagonal through entire arena
    dec = (dDist<.10) %threshold for change in range to cricket
    speed= mouseV{useFilt(vid)}>=5;
    az=azDeg{useFilt(vid)}; %too many missing points - heading looks weird
    %azFilt=interpNan(az,3,'linear');
    heading=(az<90) &(az>-90) ;
    nframe=min(length(speed),length(heading));
    nframe=min(nframe,length(dec));
    speed=speed(1:nframe)';heading=heading(1:nframe)';dec=dec(1:nframe)'; distThresh=distThresh(1:nframe)';
    approach =dec==1&(speed==1)&heading==1
    plot(d2cr,'b'); hold on;
    plot(find(approach),d2cr(approach),'og');
    appEpoch{vid,:}=(approach); axis square
    filtApp(vid,:)=sum(appEpoch{vid})
%   plot(d2cr(mouseV{useFilt(vid)}>5),'g');
end


%%
clear dTAll dTRALL dTLAll dPRAll dPLAll

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
%     subplot(rownum,colnum,vid);
    plot(dT(1:15:end),dpR(1:15:end),'.');axis square; hold on
    title('dHead theta, dPhi R');
        x = linspace(-.6,.60);
    y = linspace(-80,80);
    plot(-x,y); xlim([-.6 .6]); ylim([-80 80]);
        use = (appEpoch{vid});
plot(dT(use(1:15:end)),dpR(use(1:15:end)),'.g');
      if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(2);
%        subplot(rownum,colnum,vid);
    plot(dT(1:15:end),dpL(1:15:end),'.'); axis square; hold on
    title('dHead theta, dPhi L');
        x = linspace(-.6,.60);
    y = linspace(-80,80);
    plot(x,y);  xlim([-.6 .6]); ylim([-80 80]);
    plot(dT(use(1:15:end)),dpL(use(1:15:end)),'.g');
  if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(3);
%subplot(rownum,colnum,vid);
    plot(dT(1:15:end),dtR(1:15:end),'.');axis square; hold on
   title('dHead theta, dtheta R');
       x = linspace(-.6,.60);
    y = linspace(-80,80);  xlim([-.6 .6]); ylim([-80 80]);
    plot(-x,y);
    plot(dT(use(1:15:end)),dtR(use(1:15:end)),'.g');
  if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(4);%subplot(rownum,colnum,vid);
    plot(dT(1:15:end),dtL(1:15:end),'.');axis square; hold on;
    title('dHead theta, dtheta L');
       x = linspace(-.6,.60);  xlim([-.6 .6]); ylim([-50 50]);
    y = linspace(-50,50);
    plot(x,y);
        plot(dT(use(1:15:end)),dtL(use(1:15:end)),'.g');
  if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
        figure(5);%subplot(rownum,colnum,vid);
    plot(dtR(1:15:end),dtL(1:15:end),'.');axis square; hold on;
    title('r theta, l theta');
    x = linspace(-80,80);
    y = linspace(-50,50);
    plot(-x,y); xlim([-80 80]); ylim([-50 50]);
            plot(dtR(use(1:15:end)),dtL(use(1:15:end)),'.g');
  if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
        figure(6);%subplot(rownum,colnum,vid);
    plot(dpR(1:15:end),dpL(1:15:end),'.');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('right phi, left phi');
                plot(dpR(use(1:15:end)),dpL(use(1:15:end)),'.g');
  if vid==(useFilt(end))
        
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
% if length(dT)>500
% dTAll(vid,:)=datasample(dT,500);
% dTRAll(vid,:)=datasample(dtR,500);
% dTLAll(vid,:)=datasample(dtL,500);
% dPRAll(vid,:)=datasample(dpR,500);
% dPLAll(vid,:)=datasample(dpL,500);
% else
% end
end


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
       use=appEpoch{vid};
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
title('between eye correlation, dTheta & dPhi');
      if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
% corrRAll = cell2mat(corrRAll)'; corrLAll = cell2mat(corrLAll)';
% corrRAAll = cell2mat(corrRAAll)'; corrLAAll = cell2mat(corrLAAll)';

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

legend(L,{'dTheta','dPhi','dTh Approach','dPhi app'}); title('mean between eye corr, dtheta & dphi');


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

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll use nonapp


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
%%
clear dTAll dTRALL dTLAll dPRAll dPLAll

dTRAll=[];

% figure%('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    clear dT dpR dpL dtR dtL
    %     subplot(rownum,colnum,vid);
    nframe = min(length(dTheta{useFilt(vid)}),length(dphiR{useFilt(vid)}));
    nframe = min(nframe, length(dphiL{useFilt(vid)}));
    nframe = min(nframe, length(dthetaR{useFilt(vid)}));
    nframe = min(nframe,length(dthetaL{useFilt(vid)}));
    
    nonapp=appEpoch{useFilt(vid)}==0;
    
    dT=dTheta{useFilt(vid)}(nonapp==1); dpR=dphiR{useFilt(vid)}(nonapp==1); dpL=dphiL{useFilt(vid)}(nonapp==1);
    dtR=dthetaR{useFilt(vid)}(nonapp==1); dtL=dthetaL{useFilt(vid)}(nonapp==1);
    figure(1);
    %     subplot(rownum,colnum,vid);
    plot(dT(1:15:end),dpR(1:15:end),'b.');axis square; hold on
    title('dHead theta, dPhi R');
    x = linspace(-.6,.60);
    y = linspace(-80,80);
    plot(-x,y); %xlim([-.6 .6]); ylim([-80 80]);
    use = (appEpoch{vid});
    plot(dT(use(1:15:end)),dpR(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(2);
    %        subplot(rownum,colnum,vid);
    plot(dT(1:15:end),dpL(1:15:end),'b.'); axis square; hold on
    title('dHead theta, dPhi L');
    x = linspace(-.6,.60);
    y = linspace(-80,80);
    plot(x,y);  %xlim([-.6 .6]); ylim([-80 80]);
    plot(dT(use(1:15:end)),dpL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(3);
    %subplot(rownum,colnum,vid);
    plot(dT(1:15:end),dtR(1:15:end),'b.');axis square; hold on
    title('dHead theta, dtheta R');
    x = linspace(-.6,.60);
    y = linspace(-80,80); % xlim([-.6 .6]); ylim([-80 80]);
    plot(-x,y);
    plot(dT(use(1:15:end)),dtR(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(4);%subplot(rownum,colnum,vid);
    plot(dT(1:15:end),dtL(1:15:end),'b.');axis square; hold on;
    title('dHead theta, dtheta L');
    x = linspace(-.6,.60); % xlim([-.6 .6]); ylim([-50 50]);
    y = linspace(-50,50);
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
    plot(-x,y);% xlim([-80 80]); ylim([-50 50]);
    plot(dtR(use(1:15:end)),dtL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(6);%subplot(rownum,colnum,vid);
    plot(dpR(1:15:end),dpL(1:15:end),'b.');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y); % xlim([-80 80]); ylim([-80 80]);
    title('right phi, left phi');
    plot(dpR(use(1:15:end)),dpL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
end
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

%%clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(dTheta{useFilt(vid)}),length(dphiR{useFilt(vid)}));
    nframe = min(nframe, length(dphiL{useFilt(vid)}));
    dP=dTheta{useFilt(vid)}(1:nframe); dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
    nonapp=appEpoch{useFilt(vid)}==0;
    use = nonapp& ~isnan(dP(1:nframe));
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
title('nonapp and app, head theta and dphi');
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

legend(L,{'dPhi R non app','dPhi L non app','dPhi R Approach','dPhi L Approach'});
title('head phi, both eyes - non app and app');


if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(dTheta{useFilt(vid)}),length(dphiR{useFilt(vid)}));
    nframe = min(nframe, length(dphiL{useFilt(vid)}));
    dP=dTheta{useFilt(vid)}(1:nframe); dpR=dphiR{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
    nonapp=appEpoch{useFilt(vid)}==0;
    use = nonapp& ~isnan(dP(1:nframe));
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
title('nonapp and app, head theta and dphi');
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

legend(L,{'dPhi R non app','dPhi L non app','dPhi R Approach','dPhi L Approach'});
title('head phi, both eyes - non app and app');


if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll


figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(dTheta{useFilt(vid)}),length(dthetaR{useFilt(vid)}));
    nframe = min(nframe, length(dthetaL{useFilt(vid)}));
    dP=dTheta{useFilt(vid)}(1:nframe); dpR=dthetaR{useFilt(vid)}(1:nframe); dpL=dthetaL{useFilt(vid)}(1:nframe);
    nonapp=appEpoch{useFilt(vid)}==0;
    use = nonapp==1 & ~isnan(dP(1:nframe));
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

legend(L,{'dtheta R non app','dtheta L non app','dtheta R Approach','dtheta L Approach'}); title('head theta, both eyes, nonapp and app');


if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%
if savePDF
    pSname='T:\PreyCaptureAnalysis\Data\';
    filen=sprintf('%s',ani,'AnalyzedTS_AllVideos_withApproach_test2','.pdf')
    pdfilename=fullfile(pSname,filen);
    dos(['ps2pdf ' psfilename ' ' pdfilename]);
    delete(psfilename);
end

afilename=sprintf('%s',ani,'Analyzed','.mat')
save(fullfile(pSname, afilename))
