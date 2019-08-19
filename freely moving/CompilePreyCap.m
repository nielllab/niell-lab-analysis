clear all; close all
load('J463b.mat'); 
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
    
            
    mouse_xy{i,1,1}=Data(i).mouse_xy;
    mouseV{i,:}= Data(i).mouseV; %in pix/frame
    mouseV{i,:}=((Data(i).mouseV)/27)*30; %cm/sec
    cricket_xy{i,1}=Data(i).cricketxy;
    cricketV{i,:}= Data(i).cricketV % pix/frame
    cricketV{i,:} = ((Data(i).cricketV)/27)*30; %now cm/sec
    theta{i,:}= Data(i).theta;
    dTheta{i,:}= diff(Data(i).theta);
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
    
    dthetaR{i,:} =diff(Data(i).Rtheta);
    dthetaL{i,:} =diff(Data(i).Ltheta);
    dphiR{i,:} =diff(Data(i).Rphi);
    dphiL{i,:} =diff(Data(i).Lphi);
    
    XRcent{i,:} =Data(i).XRcent;
    YRcent{i,:} =Data(i).YRcent;
    XLcent{i,:} =Data(i).XLcent;
    YLcent{i,:} =Data(i).YLcent;
 
    
    longR{i,1}=EllipseParamsR{i,1}(:,3);longL{i,1}=EllipseParamsL{i,1}(:,3);
    shortR{i,1}=EllipseParamsR{i,1}(:,4);shortL{i,1}=EllipseParamsL{i,1}(:,4);
    
    radR{i,1}= (longR{i,1}+shortR{i,1})./length(longR{i,1});
    radL{i,1}= (longL{i,1}+shortL{i,1})./length(longL{i,1});
    pupilRvel{i,1}=diff(radR{i,1}); pupilLvel{i,1}=diff(radL{i,1})
    
    tsData{i,1}=~isempty(Data(i).TopTs);
    
end

tsData= cell2mat(tsData)
delayFull=cell2mat(slip);

% useL = (delayFull(:,2)<=3 & delayFull(:,2)>=-3);
% useR = (delayFull(:,1)<=3 & delayFull(:,1)>=-3);
% useE = (delayFull(:,3)<=3 & delayFull(:,3)>=-3);

useTime = (tsData==1) & goodTheta>=.7;%|(useL & useR)
useFilt=find(useTime)

%rownum=6; colnum=7
rownum=round(sqrt(length(useFilt)+4))
colnum=round(sqrt(length(useFilt)));
%%
figure('units','normalized','outerposition',[0 0 1 1])
for vid = 1:length(useFilt)
    subplot(rownum,colnum,vid)  ;
    bar([mean(isnan(mouse_xy{useFilt(vid),1}(1,:))) mean(isnan(cricket_xy{useFilt(vid),1}(1,:)))])
    ylabel('% error'); xlim([0.5 2.5]); ylim([0 1])
    set(gca,'XTick',[1 2])
    set(gca,'XTickLabel',{'mouse','crick'})
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
    HA=dTheta{useFilt(vid)}(1:nframe); velPhiR=dthetaR{useFilt(vid)}(1:nframe); velPhiL=dthetaL{useFilt(vid)}(1:nframe);
    use = ~isnan(HA(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(HA(use),velPhiR(use),'coeff');
    plot(lagsR/30,corrR);xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(HA(use),velPhiL(use),'coeff');
    plot(lagsL/30,corrL);xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else end
    if sum(uselagsR)==61 & sum(uselagsL)==61
        corrRAll{vid,1}=corrR(uselagsR); corrLAll{vid,1}=corrL(uselagsL);
    else
    end
end
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

corrRAll = cell2mat(corrRAll); corrLAll = cell2mat(corrLAll);
figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll),errL,'-r',1);
plot([30,30],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.25 .25]); xlim([20 40]);
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
    HA=dTheta{useFilt(vid)}(1:nframe); velPhiR=dphiR{useFilt(vid)}(1:nframe); velPhiL=dphiL{useFilt(vid)}(1:nframe);
    use = ~isnan(HA(1:nframe));
    if sum(use)>3
    [corrR lagsR]= xcorr(HA(use),velPhiR(use),'coeff');
    plot(lagsR/30,corrR);xlim([-.3 .3]);hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= xcorr(HA(use),velPhiL(use),'coeff');
    plot(lagsL/30,corrL);xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    if sum(uselagsR)==61 & sum(uselagsL)==61
        corrRAll{vid,1}=corrR(uselagsR); corrLAll{vid,1}=corrL(uselagsL);
    else
    end
end
corrRAll = cell2mat(corrRAll); corrLAll = cell2mat(corrLAll);
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
figure('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll),errL,'-r',1); 
plot([30,30],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.1 .1]); xlim([20 40]);
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
legend(L,{'right eye','left eye'}); title('mean dTheta Head and dPhi eye corr');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%

if savePDF
    pSname='T:\PreyCaptureAnalysis\Data\';
    filen=sprintf('%s',ani,'AnalyzedTS_newNetwork_AllVids','.pdf')
    pdfilename=fullfile(pSname,filen)
    dos(['ps2pdf ' psfilename ' ' pdfilename]);
    delete(psfilename);
end

