clear all; close all
load('J462aAllVids_082919_a.mat'); 
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
    RRad{i,:}=Data(i).RRad;
    LRad{i,:}=Data(i).LRad;
    
    
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

rownum=10; colnum=6
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
% 

%%% identify approach!!!

for vid=1:length(useFilt)
deltaR = diff(range{vid})*30;
vsmooth = conv(mouseV{vid},ones(5,1)/5,'same');
dRThresh=-10; %%%cm/sec
vThresh=10;
azThresh = pi/4;  %%% pi/4 = 45 deg
approach = deltaR<dRThresh & vsmooth(1:end-1)>vThresh & abs(azT{vid}(1:end-1))<azThresh;
approach(1)=0; approach(end)=0; %%% boundary conditions

starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% find start;stop

for j= 1:length(ends)-1;  %%% stitch together approachs with gap of less than 5 frames
    if (starts(j+1)-ends(j))<5 & (range{vid}(starts(j+1))- range{vid}(ends(j)))<3
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
    nframe = min(length(dTheta{useFilt(vid)}),length(dthetaR{useFilt(vid)}));
    nframe = min(nframe, length(dthetaL{useFilt(vid)}));
    nonapp=appEpoch{useFilt(vid)}(1:nframe)==0
    dT=dTheta{useFilt(vid)}; dtR=dthetaR{useFilt(vid)}; dtL=dthetaL{useFilt(vid)};
    clear use
    use =  ~isnan(dT(1:nframe))&(nonapp==1)';
  %  if sum(use)>3
    [corrR lagsR]= xcorr(dT(use),dtR(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(dT(use),dtL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
%     else
%     end
    clear use;
    use=appEpoch{vid}==1
    [corrRA lagsRA]= xcorr(dT(use),dtR(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    [corrLA lagsLA]= xcorr(dT(use),dtL(use),'coeff');
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
%     nframe = min(length(dTheta{useFilt(vid)}),length(dthetaR{useFilt(vid)}));
%     nframe = min(nframe, length(dthetaL{useFilt(vid)}));
    nonapp=appEpoch{useFilt(vid)}==0
    dT=dTheta{useFilt(vid)}; dpR=dphiR{useFilt(vid)}; dpL=dphiL{useFilt(vid)};
    clear use
    use = (nonapp==1)' & ~isnan(dT(1:length(dpR)));
    if sum(use)>3
    [corrR lagsR]= xcorr(dT(use),dpR(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(dT(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1
    [corrRA lagsRA]= xcorr(dT(use),dpR(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    [corrLA lagsLA]= xcorr(dT(use),dpL(use),'coeff');
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

legend(L,{'dPhi R non-app','dPhi L non-app','dPhi R Approach','dPhi L Approach'}); title('head Theta and Eye Phi, both eyes');
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
    use = (nonapp==1)'& ~isnan(dtR(1:nframe));
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
    nframe = min(length(dthetaR{useFilt(vid)}),length(dphiR{useFilt(vid)}));
    dtR=dthetaR{useFilt(vid)}(1:nframe); dpR=dphiR{useFilt(vid)}(1:nframe);
    nonapp=appEpoch{useFilt(vid)}==0;
    use = (nonapp==1)'& ~isnan(dtR(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(dtR(use),dpR(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    clear nframe
    nframe = min(length(dthetaL{useFilt(vid)}),length(dphiL{useFilt(vid)}));
    dtL=dthetaL{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
    [corrL lagsL]= xcorr(dtL(use),dpL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
      % use = ~isnan(dT(1:nframe));
      clear use
      use=appEpoch{vid}==1;
%     if sum(use)>3
    [corrRA lagsRA]= xcorr(dtR(use),dpR(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
%     clear nframe
    nframe = min(length(dphiR{useFilt(vid)}),length(dphiL{useFilt(vid)}));
    dtL=dthetaL{useFilt(vid)}(1:nframe); dpL=dphiL{useFilt(vid)}(1:nframe);
    [corrLA lagsLA]= xcorr(dtL(use),dpL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
   
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
       corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

    else
    end
end
title('dTheta vs dPhi corr, R and L eyes');
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

legend(L,{'R dTh/dPhi','L dTh/dPhi','R dTh/dPhi Approach','L dTh/dPhi Approach'}); title('mean dtheta & dphi, each eye');

if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nonapp=appEpoch{useFilt(vid)}==0
    nframe = min(length(dTheta{useFilt(vid)}),length(thetaR{useFilt(vid)}));
    nframe = min(nframe, length(thetaL{useFilt(vid)}));
    nframe=min(nframe,length(nonapp));
     
    dT=dTheta{useFilt(vid)}(1:nframe); tR=thetaR{useFilt(vid)}(1:nframe); tL=thetaL{useFilt(vid)}(1:nframe);
    clear use
    use = ~isnan(dT(1:nframe))& (nonapp==1)';
    if sum(use)>3
    [corrR lagsR]= xcorr(dT(use),tR(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(dT(use),tL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=(appEpoch{vid}==1)' & ~isnan(dT(1:nframe));
      [corrRA lagsRA]= xcorr(dT(use),tR(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
    [corrLA lagsLA]= xcorr(dT(use),tL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    
    
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        
    else
    end
end
title('dHead & theta Position');
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
title('dHead and theta position');
legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'});
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nonapp=appEpoch{useFilt(vid)}==0
    nframe = min(length(dTheta{useFilt(vid)}),length(phiR{useFilt(vid)}));
    nframe = min(nframe, length(phiL{useFilt(vid)}));
    nframe=min(nframe,length(nonapp));
     
    dT=dTheta{useFilt(vid)}(1:nframe); pR=phiR{useFilt(vid)}(1:nframe); pL=phiL{useFilt(vid)}(1:nframe);
    clear use
    use = ~isnan(dT(1:nframe))& (nonapp==1)';
    if sum(use)>3
    [corrR lagsR]= xcorr(dT(use),pR(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(dT(use),pL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=(appEpoch{vid}==1)' & ~isnan(dT(1:nframe));
      [corrRA lagsRA]= xcorr(dT(use),pR(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
    [corrLA lagsLA]= xcorr(dT(use),pL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    
    
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        
    else
    end
end
title('dHead & phi Position');
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
title('dHead and phi position');
legend(L,{'phi R non-app','phi L non-app','phi R Approach','phiL Approach'});
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
corrRAAll=[]
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(thetaR{useFilt(vid)}),length(thetaL{useFilt(vid)}))-1;
    tR=thetaR{useFilt(vid)}(1:nframe); tL=thetaL{useFilt(vid)}(1:nframe);
    nonapp=appEpoch{useFilt(vid)}(1:nframe)==0;
    
    use = (nonapp==1)' & ~isnan(tR(1:length(nonapp)));
    if sum(use)>3
    [corrR lagsR]= xcorr(tR(use),tL(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    clear nframe
    nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
    pR=phiR{useFilt(vid)}(1:nframe); pL=phiL{useFilt(vid)}(1:nframe);
    [corrL lagsL]= xcorr(pR(use),pL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
      % use = ~isnan(dT(1:nframe));
       use=appEpoch{vid}==1;
%     if sum(use)>3
    [corrRA lagsRA]= xcorr(tR(use),tL(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
%     clear nframe
    nframe = min(length(phiR{useFilt(vid)}),length(phiL{useFilt(vid)}));
    pR=phiR{useFilt(vid)}(1:nframe); pL=phiL{useFilt(vid)}(1:nframe);
    [corrLA lagsLA]= xcorr(pR(use),pL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
   
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
       corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

    else
    end
end
title('between eye corr, theta & phi position - non app & approach');
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

legend(L,{'theta non-app','Phi non-app','Th Approach','Phi app'}); title('mean between eye corr, theta & phi position');


if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
corrRAAll=[]
figure('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    subplot(rownum,colnum,vid);
    nframe = min(length(thetaR{useFilt(vid)}),length(phiR{useFilt(vid)}))-1;
    tR=thetaR{useFilt(vid)}(1:nframe); pR=phiR{useFilt(vid)}(1:nframe);
    nonapp=appEpoch{useFilt(vid)}(1:nframe)==0;
    
    use = (nonapp==1)' & ~isnan(tR(1:length(nonapp)));
    if sum(use)>3
    [corrR lagsR]= xcorr(tR(use),pR(use),'coeff');
    plot(lagsR/30,corrR,'b');xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    clear nframe
    nframe = min(length(thetaL{useFilt(vid)}),length(phiL{useFilt(vid)}));
    tL=thetaL{useFilt(vid)}(1:nframe); pL=phiL{useFilt(vid)}(1:nframe);
    [corrL lagsL]= xcorr(tL(use),pL(use),'coeff');
    plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    % use = ~isnan(dT(1:nframe));
    clear use
    use=appEpoch{vid}==1;
    %     if sum(use)>3
    [corrRA lagsRA]= xcorr(tR(use),pR(use),'coeff');
    plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);hold on;
    uselagsRA=(lagsRA>=-30& lagsRA<=30);
    
%     clear nframe
    nframe = min(length(thetaL{useFilt(vid)}),length(phiL{useFilt(vid)}));
    tL=thetaL{useFilt(vid)}(1:nframe); pL=phiL{useFilt(vid)}(1:nframe);
    [corrLA lagsLA]= xcorr(tL(use),pL(use),'coeff');
    plot(lagsLA/30,corrLA,'c');xlim([-.3 .3]);
    uselagsLA=(lagsLA>=-30 & lagsLA<=30);
   
    
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
       corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

    else
    end
end
title('theta vs phi corr each eye - non app & approach');
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

legend(L,{'theta & phi R eye','theta & phi L eye','th/phi R app','th/phi L app'}); title('mean theta & phi corr, each eye');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%

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
    
     figure(7);
    plot(dtR(useN(1:15:end)),dpR(useN(1:15:end)),'b.');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('R eye dTheta dPhi - nonapp & app');
    plot(dtR(use(1:15:end)),dpR(use(1:15:end)),'.g');
       if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
    figure(8);
    plot(dtL(useN(1:15:end)),dpL(useN(1:15:end)),'b.');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('L eye dTheta dPhi - nonapp & app');
    plot(dtL(use(1:15:end)),dpL(use(1:15:end)),'.g');
       if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
end

%%
clear dTAll dTRALL dTLAll dPRAll dPLAll

dTRAll=[];
close all
% figure%('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useFilt)
    clear dT dpR dpL dtR dtL
    %     subplot(rownum,colnum,vid);
    nframe = min(length(dTheta{useFilt(vid)}),length(phiR{useFilt(vid)}));
    nframe = min(nframe, length(phiL{useFilt(vid)}));
    nframe = min(nframe, length(thetaR{useFilt(vid)}));
    nframe = min(nframe,length(thetaL{useFilt(vid)}));
    dT=dTheta{useFilt(vid)}(1:nframe); pR=phiR{useFilt(vid)}(1:nframe); pL=phiL{useFilt(vid)}(1:nframe);
    tR=thetaR{useFilt(vid)}(1:nframe); tL=thetaL{useFilt(vid)}(1:nframe);
    useN= appEpoch{vid}==0;
    use = (appEpoch{vid});

    
    figure(1);
    plot(dT(useN(1:15:end)),pR(useN(1:15:end)),'b.');axis square; hold on
    title('head dtheta, phi R, nonapp & app');
    x = linspace(-40,40);
    y = linspace(-60,60);
    plot(-x,y); xlim([-40 40]); ylim([-60 60]);
    plot(dT(use(1:15:end)),pR(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(2);
    plot(dT(useN(1:15:end)),pL(useN(1:15:end)),'b.'); axis square; hold on
    title('head dtheta, phi L, nonapp & app');
    x = linspace(-40,40);
    y = linspace(-80,80);
    plot(x,y);  xlim([-40 40]); ylim([-80 80]);
    plot(dT(use(1:15:end)),pL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(3);
    plot(dT(useN(1:15:end)),tR(useN(1:15:end)),'b.');axis square; hold on
    title('head dtheta, theta R, nonapp & app');
    x = linspace(-40,40);
    y = linspace(-40,40);  xlim([-40 40]); ylim([-40 40]);
    plot(-x,y);
    plot(dT(use(1:15:end)),tR(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(4);
    plot(dT(useN(1:15:end)),tL(useN(1:15:end)),'b.');axis square; hold on;
    title('head dtheta, theta L, nonapp & app');
    x = linspace(-40,40);  xlim([-40 40]); ylim([-80 80]);
    y = linspace(-80,80);
    plot(x,y);
    plot(dT(use(1:15:end)),tL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(5);%subplot(rownum,colnum,vid);
    plot(tR(useN(1:15:end)),tL(useN(1:15:end)),'b.');axis square; hold on;
    title('r theta, l theta, nonapp & app');
    x = linspace(-80,80);
    y = linspace(-50,50);
    plot(-x,y); xlim([-80 80]); ylim([-50 50]);
    plot(tR(use(1:15:end)),tL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(6);
    plot(pR(useN(1:15:end)),pL(useN(1:15:end)),'b.');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('right phi, left phi, nonapp & app');
    plot(pR(use(1:15:end)),pL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
     figure(7);
    plot(tR(useN(1:15:end)),pR(useN(1:15:end)),'b.');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('right eye - theta/phi pos - nonapp & app');
    plot(tR(use(1:15:end)),pR(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
         figure(8);
    plot(tL(useN(1:15:end)),pL(useN(1:15:end)),'b.');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('left eye - theta/phi pos - nonapp & app');
    plot(tL(use(1:15:end)),pL(use(1:15:end)),'.g');
    if vid==(useFilt(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
end


if savePDF
    pSname='T:\PreyCaptureAnalysis\Data\';
    filen=sprintf('%s',ani,'AnalyzedTS_cleaned_082819_oldApproach','.pdf')
    pdfilename=fullfile(pSname,filen);
    dos(['ps2pdf ' psfilename ' ' pdfilename]);
    delete(psfilename);
end


afilename=sprintf('%s',ani,'Analyzed_oldApproach','.mat')
save(fullfile(pSname, afilename))