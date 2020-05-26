clear all; close all;
set(groot,'defaultFigureVisible','on') %disable figure plotting
deInter=1;
if deInter
    frRate=60;
else frRate=30; 
end

savePDF=1;
if savePDF
    psfilename = 'C:\analysisPS_all.ps';
    if exist(psfilename,'file')==2;delete(psfilename);end
end


pname = 'T:\PreyCaptureAnalysis\Data\';

analyzeAcc = input('use acc data? (1==yes) ');


% for acc sessions
if analyzeAcc==1
    files={'J462aDeinterlaced_052320_analyzed_halfShift.mat','J462bDeinterlaced_052320_analyzed_halfShift.mat','J462cDeinterlaced_052320_analyzed_halfShift.mat',...
        'J463bDeinterlaced_052320_analyzed_halfShift.mat','J463cDeinterlaced_052320_analyzed_halfShift.mat','J470cDeinterlaced_052320_analyzed_halfShift.mat',...
        'J475cDeinterlaced_052320_analyzed_halfShift.mat'};

%         files={'J462aAnalyzed_010819_ACCSessions_a','J462bAnalyzed_010819_ACCSessions_a.mat','J462cAnalyzed_010819_ACCSessions_a',...
%         'J463bAnalyzed_010819_ACCSessions_a','J463cAnalyzed_010819_ACCSessions_a','J470cAnalyzed_010819_ACCSessions_a',...
%         'J475cAnalyzed_010819_ACCSessions_a'}; %analysis redone with daylight savings data on 010819
%     
else
    files={'J462aAnalyzed_120219_allSessions.mat','J462bAnalyzed_120219_allSessions.mat','J462cAnalyzed_120219_allSessions.mat',...
        'J463bAnalyzed_120219_allSessions.mat','J463cAnalyzed_120219_allSessions.mat','J465dAnalyzed_120219_allSessions.mat'...
        'J470aAnalyzed_120219_allSessions.mat','J470bAnalyzed_120219_allSessions.mat'}
    
end

fname=[];aniname=[];cricketPos=[];mousePos=[];cricketSp=[];mouseSp=[];thetaHead=[];
d_Theta=[];az=[];dist=[];Rtheta=[];Ltheta=[];Rphi=[];Lphi=[];Rcent=[];Lcent=[];n_vid=[]; mousePos=[];theta_good=[];approach=[];useData=[];
fitParamsR=[];fitParamsL=[]; dRtheta=[];dLtheta=[]; dRphi=[];dLphi=[]; dxR=[];dyR=[];Rcent=[];Lcent=[];radR=[];radL=[];Rgood=[];ccR=[];slopeR=[];
scaleR=[]; Lgood=[];ccL=[];slopeL=[];scaleL=[]; sess=[]; date=[];clipnum=[]; cricketPosHead=[]; cricketSpHead=[];cricketTh=[]; missingR=[];missingL=[];
tiltPts=[]; rollPts=[];yawPts=[]; dlcHead=[];dlcDhead=[];dlcVg=[];dlcDvg=[];gyroCh1=[];gyroCh2=[];gyroCh3=[];dlcEyePhi=[];dlcDeyePhi=[];accelChannels=[];
accelChannelsRaw=[];accelCorrelation=[];accDrift=[];


n=0;
for i = 1:length(files)
    fname= fullfile(pname,files{i});
    load(fname,'cricket_xy','mouse_xy','cricketV','cricketH','theta', 'dTheta','azT','mouseV','goodTheta','range','theta','dTheta','appEpoch','useData','EllipseParamsR',...
        'EllipseParamsL','thetaR','thetaL','phiR','phiL','dthetaR','dthetaL','dphiR','dphiL','dXRcent','XRcent','YRcent','YLcent','XLcent','RRad','LRad',...
        'goodR','Rcc','Rslope','Rscale','Rngood','goodL','Lcc','Lslope','Lscale','Lngood','animal','sessionN','clip','expdate','cricketAz','rMissing','lMissing');
   if analyzeAcc==1  
       load(fname,'accelCorr','accelDrift','allTilt','allRoll','allYaw','dlcHth','dlcDhth','dlcVerg','dlcDverg','allGyro1','allGyro2','allGyro3','dlcPhi','dlcDphi','accelData','accelDataRaw');
   end
    nc = length(cricket_xy(:,1))'; vidrange= n+1:n+nc;
    nvid{i} = length(cricket_xy(:,1))';
    %  nvidTot = nvid(i+1)
    %   delayFull{end+1,1} = cell2mat(slip);
    
    for j = 1:length(vidrange)
        n_vid(end+1) = j;
        aniname{end+1} = animal{j};
        sess{end+1} = sessionN{j};
        date{end+1}=expdate{j};
        clipnum{end+1} =clip{j};
        
        mousePos{end+1}=mouse_xy{j,:};
        mouseSp{end+1} = mouseV{j,1};
        %         cricketPosHead{end+1} = cricket_xy{j,:};
        %         cricketSpHead{end+1} =cricketV{j,:};
        
        % cricketPos{end+1} = cricket_xy{j,:};
        cricketSp{end+1} =cricketV{j,:};
        cricketPos{end+1}=cricket_xy{j,:};
        cricketTh{end+1}=cricketAz{j,:};
        
        thetaHead{end+1,1} = theta{j,:};
        d_Theta{end+1,1,1} = dTheta{j,:};
        az{end+1,1}=azT{j,:};
        dist{end+1,1}=range{j,:};
        theta_good(end+1)=goodTheta(j);
        %         useData(end+1)=useData(j);
        %         approach{end+1}=appEpoch{j,:};
        %         Rcent{end+1,1}=x_centR{j,:}; Rcent{end,2}=y_centR{j,:};
        %         Lcent{end+1,1}=x_centL{j,:}; Lcent{end,2}=y_centL{j,:};
        if analyzeAcc==1
        accelChannels{end+1,1}  = accelData{j,:};
    %    accelCorrelation  = accelCorr(j)
   %     accDrift  = accelDrift(j)
        accelChannelsRaw{end+1,1} = accelDataRaw{j,:};
        end
        
        
        fitParamsR{end+1,1,1} = EllipseParamsR{j,:};
        fitParamsL{end+1,1,1} = EllipseParamsL{j,:};
        Rtheta{end+1,1}=thetaR{j,:};
        Ltheta{end+1,1}=thetaL{j,:};
        Rphi{end+1,1}=phiR{j,:};
        Lphi{end+1,1}=phiL{j,:};
        
        dRtheta{end+1,1}=dthetaR{j,:};
        dLtheta{end+1,1}=dthetaL{j,:};
        dRphi{end+1,1}=dphiR{j,:};
        dLphi{end+1,1}=dphiL{j,:};
        %         dxR{end+1,1} =dXRcent{j,:}
        %       %  dxR{end+1,2} =dYRcent{j,:}; %undefined from single animal compile
        
        Rcent{end+1,1} =XRcent{j,:};
        Rcent{end,2} =YRcent{j,:};
        Lcent{end+1,1} =XLcent{j,:};
        Lcent{end,2} =YLcent{j,:};
        radR{end+1,1}=RRad{j,:};
        radL{end+1,1}=LRad{j,:};
        
        Rgood{end+1,1}= Rngood{j,:};% goodR(j);%1=all 8pts above likelihood .95 %doesn't makes sense here...
        ccR(end+1)= Rcc(j);
        slopeR(end+1)= Rslope(j);
        scaleR(end+1)= Rscale(j);
        missingR(end+1)=rMissing(j);
        
        Lgood{end+1,1}= Lngood{j,:};
        ccL(end+1)= Lcc(j);
        slopeL(end+1)= Lslope(j);
        scaleL(end+1)= Lscale(j);
        missingL(end+1)=lMissing(j);
        
        
    end
    
    if analyzeAcc==1
    tiltPts{end+1,1} = allTilt;
    rollPts{end+1,1}= allRoll;
    yawPts{end+1,1}=allYaw;
    dlcHead{end+1,1}=dlcHth;
    dlcDhead{end+1,1}=dlcDhth;
    dlcVg{end+1,1}=dlcVerg;
    dlcDvg{end+1,1}=dlcDverg;
    gyroCh1{end+1,1}=allGyro1;
    gyroCh2{end+1,1}=allGyro2;
    gyroCh3{end+1,1}=allGyro3;
    dlcEyePhi{end+1,1}=dlcPhi;
    dlcDeyePhi{end+1,1}=dlcDphi;
    end
    
    
end
clear goodR goodL
goodR=ccR>.3 & missingR>.75; %prop of non-nans is more than .75
goodL=ccL>.3& missingL>.75;
useTime = theta_good>=.7 & goodR & goodL;
useData=find(useTime);
%%
figure;
subplot(2,2,1)
plot(ccR,slopeR,'o'); axis square; xlim([0 1]); ylim([0 1]); title('R');
xlabel('R corrcoef');ylabel('R Slope');
hold on; plot([0,1],[0,1]);
subplot(2,2,2)
plot(ccL,slopeL,'o');axis square; xlim([0 1]); hold on;ylim([0 1]);plot([0,1],[0,1]);
title('L');
xlabel('L corrcoef');ylabel('L Slope');
subplot(2,2,3)
plot(ccR,scaleR,'o'); axis square; xlim([0 1]); ylim([0 100]); title('R');
xlabel('R corrcoef');ylabel('R Scale');
hold on; plot([0,1],[0,100]);
subplot(2,2,4)
plot(ccL,scaleL,'o');axis square; xlim([0 1]); hold on;ylim([0 100]);plot([0,1],[0,100]);
title('L');
xlabel('L corrcoef');ylabel('L Scale');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
clear appEpoch

for vid=1:length(useData)
    deltaR = diff(dist{useData(vid)})*frRate;
    badDist=(isnan(deltaR));
    vsmooth = conv(mouseSp{useData(vid)},ones(5,1)/5,'same');
    dRThresh=-10; %%%cm/sec
    vThresh=5;
    azThresh = pi/4;  %%% pi/4 = 45 deg
    %  azThresh = pi/2;  %%% pi/2 = 90 deg
    
    approach = deltaR<dRThresh & vsmooth(1:end-1)>vThresh & abs(az{useData(vid)}(1:end-1))<azThresh;
    approach(1)=0; approach(end)=0; %%% boundary conditions
    
    if sum(approach(badDist))>1
        approach(badDist)=NaN;
    end
    
    starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% find start;stop
    
    for j= 1:length(ends)-1;  %%% stitch together approachs with gap of less than 5 frames
        if (starts(j+1)-ends(j))<5 & (dist{useData(vid)}(starts(j+1))- dist{useData(vid)}(ends(j)))<3
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
%%

if analyzeAcc==1

tiltAll=[];rollAll=[];yawAll=[]; vergDlc=[];dvergDlc=[];gyro1All=[];gyro2All=[];gyro3All=[];
dvergDlc=[]; dlcDphi=[]; dlcPhi=[];dlcDhth=[];dlcHth=[];appAll=[];
% figure
for i = 1:length(useData)
%     figure('units','normalized','outerposition',[0 0 1 1])
    appT=appEpoch{(i)};
    appAll=[appAll appT];
    roll = (accelChannels{useData(i)}(:,1)); roll=roll-nanmean(roll);
    rollFilt = medfilt1(roll,2);
    tilt = (accelChannels{useData(i)}(:,2)); tilt=tilt-nanmean(tilt);
    tiltFilt = medfilt1(tilt,8);
    yaw = (accelChannels{useData(i)}(:,3)); yaw=yaw-nanmean(yaw);
    yawFilt = medfilt1(yaw,10);
    verg=(Rtheta{useData(i)}-Ltheta{useData(i)});
    dverg=(dRtheta{useData(i)}-dLtheta{useData(i)});
  
    gyro3=(accelChannels{useData(i)}(:,6));
    phi=(Rphi{useData(i)}-Lphi{useData(i)}); phi=phi-nanmean(phi);
    dphi=(dRphi{useData(i)}-dLphi{useData(i)}); dphi=dphi-nanmean(dphi);
    badE=dphi>15| dphi<-15 ;
    dphi(badE)=nan;
    hth=thetaHead{useData(i)}; hth=hth-nanmean(hth);
    dhth=d_Theta{useData(i)}; dhth=dhth-nanmean(dhth);
    bad=dhth>20 | dhth<-20;
    dhth(bad)=nan;
%     
%     subplot(2,3,1)
%     plot(rollFilt); hold on; axis square;
%     plot(phi);
%     title('acc 1 - roll, phi')
%     subplot(2,3,2)
%     plot(tiltFilt); hold on; axis square;plot(verg);
%     title('acc 2 - tilt, vg')
%     subplot(2,3,3)
%     plot(gyro3); hold on; axis square
%     plot(dhth); title('gyro 3, dHead Th')
%     subplot(2,3,4);
%     plot(-frRate:frRate,nanxcorr(rollFilt,phi,frRate,'coeff'))
%     title('roll, phi diff');ylim([-1 1]);
%     subplot(2,3,5);
%     plot(-frRate:frRate,nanxcorr(tiltFilt,verg,frRate,'coeff'))
%     title('tilt, vergence');ylim([-1 1]);
% 
%     subplot(2,3,6)
%     plot(-frRate:frRate,nanxcorr(gyro3,dhth,frRate,'coeff'))
%     title('gyro 3, dHead th');ylim([-1 1]); 
    
    tiltAll=[tiltAll tiltFilt(1:end-1)'];
    rollAll=[rollAll rollFilt(1:end-1)'];
    yawAll=[yawAll yawFilt(1:end-1)'];
    if deInter
        vergDlc=[vergDlc verg(1:end-1)];
        dvergDlc=[dvergDlc dverg];
        dlcPhi=[dlcPhi phi(1:end-1)];
        dlcDphi=[dlcDphi dphi];
    else
        vergDlc=[vergDlc verg(1:end-1)'];
        dvergDlc=[dvergDlc dverg'];
        dlcPhi=[dlcPhi phi(1:end-1)'];
        dlcDphi=[dlcDphi dphi'];
    end
    gyro3All=[gyro3All gyro3(1:end-1)'];

    dlcHth=[dlcHth hth(1:end-1)];
    dlcDhth=[dlcDhth dhth(1:end-1)'];
%     if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end
% 
% figure
% plot(tiltAll,vergDlc,'.'); axis equal; hold on; lsline
% mdl = fitlm(tiltAll,vergDlc)
% axis([-90 90 -90 90]);
% xlabel('acc pitch'); ylabel('vergence')
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% figure
% plot(rollAll,dlcPhi,'.'); axis equal; hold on; lsline
% axis([-90 90 -90 90]);
% xlabel('acc roll'); ylabel('eye phi')
% mdl = fitlm(rollAll,dlcPhi)
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
% use=find(randsample(rollAll,15000)); 
%%
close all
for c=0:1
    clear use
use = find(appAll==c);
figure(1)
plot(rollAll(use),dlcPhi(use),'.'); axis equal; hold on; lsline
xlabel('acc roll'); ylabel('eye phi');
ylim([-90 90]); xlim([-90 90]);
R=corrcoef(rollAll(use), dlcPhi(use),'Rows','pairwise')
text(-80,(80-c*5), ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
title('acc ch 1');
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

end
figure(2)
plot(tiltAll(use),vergDlc(use),'.'); axis equal; hold on; lsline
xlabel('acc tilt'); ylabel('eye vergence');
ylim([-90 90]); xlim([-90 90]);
R=corrcoef(tiltAll(use), vergDlc(use),'Rows','pairwise')
text(-80,(80-c*5), ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
title('acc ch 2');
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end
% test=(dlcDhth>10|dlcDhth<-10)&gyro3All>-5& gyro3All<5;

figure(3)
plot(gyro3All(use),dlcDhth(use),'.');axis equal; hold on; lsline
% % plot(gyro3All(test==1),dlcDhth(test==1),'r.');axis equal; hold on; lsline
xlabel('gyro 3'); ylabel('d head theta');
ylim([-40 40]); xlim([-40 40]);
R=corrcoef(gyro3All(use), dlcDhth(use),'Rows','pairwise')
text(-frRate,(frRate-c*5), ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
title('gyro ch 3')
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

end
% % test=gyro3All>-5& gyro3All<5  & (dlcDhth>15|dlcDhth<-15);
% plot(gyro3All(test),dlcDhth(test),'.')
% RS=corrcoef(gyro3All(test==0), dlcDhth(test==0),'Rows','pairwise')
% text(-80,70, ['corrcoef = ' num2str(RS(1,2),'%.2f')],'FontSize',10)

figure(4)
[corr lags]=(nanxcorr(rollAll(use),dlcPhi(use),frRate,'coeff'));
plot(lags,corr);axis square; hold on
title('acc roll, phi');
ylim([-1 1]);
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end

figure(5)
[corr lags]=(nanxcorr(tiltAll(use),vergDlc(use),frRate,'coeff'));
plot(lags,corr); hold on
axis square;ylim([-1 1]);
title('acc tilt, vergence');
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end


figure(6)
[corr lags]=(nanxcorr(gyro3All(use),dlcDhth(use),frRate,'coeff'));
plot(lags, corr); hold on;
% [corr lags]=(nanxcorr(gyro3All(test==0),dlcDhth(test==0),frRate,'coeff'));
% plot(lags, corr,'r')
axis square;ylim([-1 1]);
title('gyro ch 3,d head theta')
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

end
end
%%
figure
for c=0:1
    use = find(appAll==c);

    rollPhi=[rollAll(use);dlcPhi(use)];
    z1 = hist3(rollPhi',{-90:1.5:90, -90:1.5:90});
    subplot(3,2,c+1)
    imagesc(imresize(z1,10,'bilinear'));
    axis square
    xlabel('acc roll'); ylabel('diff in phi');
    if c==0 title('non-approach')
    else
        title('approach')
    end
    
    tiltVerg=[tiltAll(use);vergDlc(use)];
    z2 = hist3(tiltVerg',{-90:1.5:90, -90:1.5:90});
    subplot(3,2,c+3)
    imagesc(imresize(z2,10,'bilinear'));
    axis square
    xlabel('acc tilt'); ylabel('eye vergence');

    
    gyroYaw=[gyro3All(use);dlcDhth(use)];
    z3 = hist3(gyroYaw',{-25:1.5:25, -25:1.5:25});
    subplot(3,2,c+5)
    imagesc(imresize(z3,10,'bilinear'));
    axis square
    xlabel('gyro 3'); ylabel('dlc head th');
end
end


%%
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    %subplot(rownum,colnum,vid);
    nframe = min(length(d_Theta{useData(vid)}),length(Rtheta{useData(vid)}));
    nframe = min(nframe, length(Ltheta{useData(vid)}));
    nframe=min(nframe, length(appEpoch{vid}));
    nonapp=appEpoch{vid}(1:nframe)==0;
    dT=d_Theta{useData(vid)}; dtR=Rtheta{useData(vid)}; dtL=Ltheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    %  if sum(use)>3
    [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
    %  plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
    hold on;
    uselagsR=(lagsR>=-frRate& lagsR<=frRate);

    [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
    % plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-frRate & lagsL<=frRate);
    %     else
    %     end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
        %  plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
        %  plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
        uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
    else
    end
    if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

    else
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
% subplot(1,2,1)
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]);
ylim([-.3 .3]); xlim([21 41]);
axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('d_head theta, both eyes');
% %%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% for vid=1:length(useData)
%     %     nframe = min(length(dTheta{useData(vid)}),length(dthetaR{useData(vid)}));
%     %     nframe = min(nframe, length(dthetaL{useData(vid)}));
%     nonapp=appEpoch{vid}==0;
%     dT=d_Theta{useData(vid)}; dpR=dRphi{useData(vid)}; dpL=dLphi{useData(vid)};
%     clear use
%     use = (nonapp==1)'% & ~isnan(dT(1:length(dpR)));
%     if sum(use)>3
%         [corrR lagsR]= nanxcorr(dT(use),dpR(use),frRate,'coeff');
%         uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%         [corrL lagsL]= nanxcorr(dT(use),dpL(use),frRate,'coeff');
%         uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=appEpoch{vid}==1;
%     if sum(use)>3 &sum(~isnan(dpR(use)))>20
%         [corrRA lagsRA]= nanxcorr(dT(use),dpR(use),frRate,'coeff');
%         uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%         [corrLA lagsLA]= nanxcorr(dT(use),dpL(use),frRate,'coeff');
%         uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%     else
%     end
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
% % if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %subplot(1,2,2);
% figure
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'dPhi R non-app','dPhi L non-app','dPhi R Approach','dPhi L Approach'}); title('d_head Theta and Eye Phi, both eyes');
%  if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %%
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

    if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

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

plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.7 .7]);
xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dTheta non-app','dPhi non-app','dTh Approach','dPhi app'}); title('mean between eye corr, dtheta & dphi');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
% 
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    nframe = min(length(az{useData(vid)}),length(dRtheta{useData(vid)}));
    nframe = min(nframe, length(dLtheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe)==0;
    dT=az{useData(vid)}; dtR=dRtheta{useData(vid)}; dtL=dLtheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
        [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
         plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
        hold on;
        uselagsR=(lagsR>=-frRate& lagsR<=frRate);

        [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
        plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
        uselagsL=(lagsL>=-frRate & lagsL<=frRate);
    else
    end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
         plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
         plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
        uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
    else
    end
    if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);

    else
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.25 .2]); xlim([1 (2*(frRate)+1)]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('az, eye theta both eyes');

title('az and dEye Theta')

% %%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
%
% for vid=1:length(useData)
%     nframe = min(length(diff(az{useData(vid)})),length(dRtheta{useData(vid)}));
%     nframe = min(nframe, length(dLtheta{useData(vid)}));
%     nonapp=appEpoch{vid}(1:nframe)==0
%     dT=diff(az{useData(vid)}); dtR=dRtheta{useData(vid)}; dtL=dLtheta{useData(vid)};
%     clear use
%     use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
%     if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
%         [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
%         hold on;
%         uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%         [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         % plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%         uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=appEpoch{vid}==1;
%     if sum(use)>4 & sum(~isnan(dtR(use)))>20
%         [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
%         uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%         [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         %  plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
%         uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%     else
%     end
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
% subplot(1,2,2)
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.3 .3]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('d_az, eye theta both eyes');
%  if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %%
%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
%
% for vid=1:length(useData)
%     nframe = min(length(az{useData(vid)}),length(Rtheta{useData(vid)}));
%     nframe = min(nframe, length(Ltheta{useData(vid)}));
%     nonapp=appEpoch{vid}(1:nframe-1)==0;
%     dT=az{useData(vid)}; dtR=Rtheta{useData(vid)}; dtL=Ltheta{useData(vid)};
%     clear use
%     use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
%     if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
%         [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
%         hold on;
%         uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%         [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         % plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%         uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=appEpoch{vid}==1;
%     if sum(use)>4 & sum(~isnan(dtR(use)))>20
%         [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
%         uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%         [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         %  plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
%         uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%     else
%     end
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
%
% figure('units','normalized','outerposition',[0 0 1 1])
% % subplot(1,2,1)
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.2 .2]); xlim([1 60]);
% axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'}); title('az, eye theta pos both eyes');
%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% %%
% for vid=1:length(useData)
%     nframe = min(length(diff(az{useData(vid)})),length(Rtheta{useData(vid)}));
%     nframe = min(nframe, length(Ltheta{useData(vid)}));
%     nonapp=appEpoch{vid}(1:nframe)==0
%     dT=diff(az{useData(vid)}); dtR=Rtheta{useData(vid)}; dtL=Ltheta{useData(vid)};
%     clear use
%     use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
%     if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
%         [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
%         hold on;
%         uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%         [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         % plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%         uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=appEpoch{vid}==1;
%     if sum(use)>4 & sum(~isnan(dtR(use)))>20
%         [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
%         uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%         [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         %  plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
%         uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%     else
%     end
%     if sum(uselagsR)==(2*(frRate)+1) && sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
% subplot(1,2,2)
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.3 .3]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'}); title('d_az, eye theta position both eyes');
%  if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%
% % %%
% % close all
% % for vid=1:length(useData)
% %
% % nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
% % nframe = min(nframe, length(az{useData(vid)}));
% %
% % tR=Rtheta{useData(vid)}(21:end); tL=Ltheta{useData(vid)}(21:end); azC=(az{useData(vid)}(1:end-20));
% % tR=tR-nanmean(tR); tL=tL-nanmean(tL); azC=azC-nanmean(azC);
% % % tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe); azC=(az{useData(vid)}(1:nframe));
% % useN= appEpoch{vid}==0;
% %
% % %useN= appEpoch{vid}(1:end-4)==0;
% % figure(1)
% % plot(azC(useN(1:15:end)),tR(useN(1:15:end)),'bo'); axis square; hold on;xlim([-pi pi]);
% % title('az, R eye Theta');
% % figure(2);
% % plot(azC(useN(1:15:end)),tL(useN(1:15:end)),'bo');axis square; hold on;xlim([-pi pi]);
% % title('az, L eye Theta');
% % end
% %
% % for vid=1:length(useData)
% %
% % nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
% % nframe = min(nframe, length(az{useData(vid)}));
% % % tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe); azC=(az{useData(vid)}(1:nframe));
% % tR=Rtheta{useData(vid)}(21:end); tL=Ltheta{useData(vid)}(21:end); azC=(az{useData(vid)}(1:end-20));
% %
% % tR=tR-nanmean(tR); tL=tL-nanmean(tL); azC=azC-nanmean(azC);
% %
% % % tR=Rtheta{useData(vid)}(5:end); tL=Ltheta{useData(vid)}(5:end); azC=(az{useData(vid)}(1:end-4));
% % use= appEpoch{vid}==1;
% % % use= appEpoch{vid}(1:end-4)==1;
% %
% % figure(1);
% % plot(azC(use(1:15:end)),tR(use(1:15:end)),'go');axis square; hold on; %xlim([-pi pi]);
% % title('az, R eye Theta');
% % figure(2);
% % plot(azC(use(1:15:end)),tL(use(1:15:end)),'go');axis square; hold on;xlim([-pi pi]);
% % title('az, L eye Theta');
% % end
% %
%
%
%
%
% %%
%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
%
% for vid=1:length(useData)
%     nframe = min(length(cricketSp{useData(vid)}),length(dRtheta{useData(vid)}));
%     nframe = min(nframe, length(dLtheta{useData(vid)}));
%     nonapp=appEpoch{vid}(1:nframe)==0;
%     dT=cricketSp{useData(vid)}; dtR=dRtheta{useData(vid)}; dtL=dLtheta{useData(vid)};
%     clear use
%     use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
%     if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
%         [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
%         hold on;
%         uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%         [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         % plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%         uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=appEpoch{vid}==1;
%     if sum(use)>4 & sum(~isnan(dtR(use)))>20
%         [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
%         uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%         [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         %  plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
%         uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%     else
%     end
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
%
% figure('units','normalized','outerposition',[0 0 1 1])
% subplot(1,2,1)
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.1 .1]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('cricketSp, eye theta both eyes');
%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
%
% for vid=1:length(useData)
%     nframe = min(length(diff(cricketSp{useData(vid)})),length(dRtheta{useData(vid)}));
%     nframe = min(nframe, length(dLtheta{useData(vid)}));
%     nonapp=appEpoch{vid}(1:nframe)==0;
%     dT=diff(cricketSp{useData(vid)}); dtR=dRtheta{useData(vid)}; dtL=dLtheta{useData(vid)};
%     clear use
%     use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
%     if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
%         [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
%         hold on;
%         uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%         [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         % plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%         uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=appEpoch{vid}==1;
%     if sum(use)>4 & sum(~isnan(dtR(use)))>20
%         [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
%         uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%         [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         %  plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
%         uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%     else
%     end
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
% subplot(1,2,2)
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.1 .1]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('d_cricketSp, eye theta both eyes');
%  if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %%
%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
%
% for vid=1:length(useData)
%     nframe = min(length(cricketSp{useData(vid)}),length(Rtheta{useData(vid)}));
%     nframe = min(nframe, length(Ltheta{useData(vid)}));
%     nonapp=appEpoch{vid}(1:nframe-1)==0
%     dT=cricketSp{useData(vid)}; dtR=Rtheta{useData(vid)}; dtL=Ltheta{useData(vid)};
%     clear use
%     use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
%     if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
%         [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
%         hold on;
%         uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%         [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         % plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%         uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=appEpoch{vid}==1;
%     if sum(use)>4 & sum(~isnan(dtR(use)))>20
%         [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
%         uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%         [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         %  plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
%         uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%     else
%     end
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
%
% figure('units','normalized','outerposition',[0 0 1 1])
% subplot(1,2,1)
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.3 .3]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'}); title('cricketSp, eye theta pos both eyes');
%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
%
% for vid=1:length(useData)
%     nframe = min(length(diff(cricketSp{useData(vid)})),length(Rtheta{useData(vid)}));
%     nframe = min(nframe, length(Ltheta{useData(vid)}));
%     nonapp=appEpoch{vid}(1:nframe)==0
%     dT=diff(cricketSp{useData(vid)}); dtR=Rtheta{useData(vid)}; dtL=Ltheta{useData(vid)};
%     clear use
%     use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
%     if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
%         [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
%         hold on;
%         uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%         [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         % plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%         uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=appEpoch{vid}==1;
%     if sum(use)>4 & sum(~isnan(dtR(use)))>20
%         [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
%         uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%         [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         %  plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
%         uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%     else
%     end
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
% subplot(1,2,2)
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.1 .1]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'}); title('d_cricketSp, eye theta position both eyes');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %%
%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
%
% for vid=1:length(useData)
%     nframe = min(length(mouseSp{useData(vid)}),length(dRtheta{useData(vid)}));
%     nframe = min(nframe, length(dLtheta{useData(vid)}));
%     nonapp=appEpoch{vid}(1:nframe)==0;
%     dT=mouseSp{useData(vid)}; dtR=dRtheta{useData(vid)}; dtL=dLtheta{useData(vid)};
%     clear use
%     use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
%     if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
%         [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
%         hold on;
%         uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%         [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         % plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%         uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=appEpoch{vid}==1;
%     if sum(use)>4 & sum(~isnan(dtR(use)))>20
%         [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
%         uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%         [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         %  plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
%         uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%     else
%     end
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
%
% figure('units','normalized','outerposition',[0 0 1 1])
% subplot(1,2,1)
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.1 .1]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
% %
% legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('mouseSp, eye theta both eyes');
%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
%
% for vid=1:length(useData)
%     nframe = min(length(diff(mouseSp{useData(vid)})),length(dRtheta{useData(vid)}));
%     nframe = min(nframe, length(dLtheta{useData(vid)}));
%     nonapp=appEpoch{vid}(1:nframe)==0;
%     dT=diff(mouseSp{useData(vid)}); dtR=dRtheta{useData(vid)}; dtL=dLtheta{useData(vid)};
%     clear use
%     use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
%     if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
%         [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
%         hold on;
%         uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%         [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         % plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%         uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=appEpoch{vid}==1;
%     if sum(use)>4 & sum(~isnan(dtR(use)))>20
%         [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
%         uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%         [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         %  plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
%         uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%     else
%     end
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
% subplot(1,2,2)
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.2 .2]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('d_mouseSp, eye theta both eyes');
%  if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %%
%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
%
% for vid=1:length(useData)
%     nframe = min(length(mouseSp{useData(vid)}),length(Rtheta{useData(vid)}));
%     nframe = min(nframe, length(Ltheta{useData(vid)}));
%     nonapp=appEpoch{vid}(1:nframe-1)==0;
%     dT=mouseSp{useData(vid)}; dtR=Rtheta{useData(vid)}; dtL=Ltheta{useData(vid)};
%     clear use
%     use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
%     if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
%         [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
%         hold on;
%         uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%         [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         % plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%         uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=appEpoch{vid}==1;
%     if sum(use)>4 & sum(~isnan(dtR(use)))>20
%         [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
%         uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%         [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         %  plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
%         uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%     else
%     end
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
%
% figure('units','normalized','outerposition',[0 0 1 1])
% subplot(1,2,1)
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.5 .5]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'}); title('mouseSp, eye theta pos both eyes');
%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
%
% for vid=1:length(useData)
%     nframe = min(length(diff(mouseSp{useData(vid)})),length(Rtheta{useData(vid)}));
%     nframe = min(nframe, length(Ltheta{useData(vid)}));
%     nonapp=appEpoch{vid}(1:nframe)==0;
%     dT=diff(mouseSp{useData(vid)}); dtR=Rtheta{useData(vid)}; dtL=Ltheta{useData(vid)};
%     clear use
%     use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
%     if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
%         [corrR lagsR]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsR/frRate,corrR,'b');%xlim([-.3 .3])
%         hold on;
%         uselagsR=(lagsR>=-frRate& lagsR<=frRate);
%
%         [corrL lagsL]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         % plot(lagsL/frRate,corrL,'r');xlim([-.3 .3]);
%         uselagsL=(lagsL>=-frRate & lagsL<=frRate);
%     else
%     end
%     clear use;
%     use=appEpoch{vid}==1;
%     if sum(use)>4 & sum(~isnan(dtR(use)))>20
%         [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),frRate,'coeff');
%         %  plot(lagsRA/frRate,corrRA,'g');xlim([-.3 .3]);
%         uselagsRA=(lagsRA>=-frRate& lagsRA<=frRate);
%         [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),frRate,'coeff');
%         %  plot(lagsLA/frRate,corrLA,'c');%xlim([-.3 .3]);
%         uselagsLA=(lagsLA>=-frRate & lagsLA<=frRate);
%     else
%     end
%     if sum(uselagsR)==(2*(frRate)+1) & sum(uselagsL)==(2*(frRate)+1) &sum(uselagsRA)==(2*(frRate)+1) & sum(uselagsLA)==(2*(frRate)+1)
%         corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
%         corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
%
%     else
%     end
% end
% subplot(1,2,2)
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
% errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.2 .2]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'}); title('d_mouseSp, eye theta position both eyes');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
% for vid = 1:length(useData)
%
%     rRad_non(vid)= nanmean(normalize(radR{useData(vid)}(appEpoch{vid}==0,1)));
%     lRad_non(vid)=nanmean(normalize(radL{useData(vid)}(appEpoch{vid}==0,1)));
%
%     rRad_app(vid)= nanmean(normalize(radR{useData(vid)}(appEpoch{vid}==1,1)));
%     lRad_app(vid)= nanmean(normalize(radL{useData(vid)}(appEpoch{vid}==1,1)));
% end




% figure;
% plot(rRad_non,lRad_non,'o');hold on;
% plot(rRad_app,lRad_app,'go'); axis square
%
% figure;plot(radR{useData(200)}); hold on; plot(find(appEpoch{200}==1),radR{useData(200)}(appEpoch{200}==1,1),'g')
% figure('units','normalized','outerposition',[0 0 1 1])
% subplot(1,2,1)
% bar([nanmean(rRad_non) nanmean(rRad_app)])
% subplot(1,2,2)
% bar([nanmean(lRad_non) nanmean(lRad_app)])

%%
clear dTAll dTRALL dTLAll dPRAll dPLAll

% dTRAll=[];
% close all
% % figure%('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useData)
%     clear dT dpR dpL dtR dtL
%     %     subplot(rownum,colnum,vid);
%     nframe = min(length(d_Theta{useData(vid)}),length(dRphi{useData(vid)}));
%     nframe = min(nframe, length(dLphi{useData(vid)}));
%     nframe = min(nframe, length(dRtheta{useData(vid)}));
%     nframe = min(nframe,length(dLtheta{useData(vid)}));
%
%     dT=d_Theta{useData(vid)}(1:nframe); pR=dRphi{useData(vid)}(1:nframe); pL=dLphi{useData(vid)}(1:nframe);
%     tR=dRtheta{useData(vid)}(1:nframe); tL=dLtheta{useData(vid)}(1:nframe);
%     useN= appEpoch{vid}==0; azC=az{useData(vid)}(1:nframe);
%     use = (appEpoch{vid});
%
%
%     figure(1);
%     plot(dT(useN(1:15:end)),pR(useN(1:15:end)),'bo');axis square; hold on
%     title('head d_Theta, phi R, nonapp & app');
%     x = linspace(-40,40);
%     y = linspace(-60,60);
%     plot(-x,y); xlim([-40 40]); ylim([-60 60]);
%     plot(dT(use(1:15:end)),pR(use(1:15:end)),'og');
%     if vid==(useData(end))
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(2);
%     plot(dT(useN(1:15:end)),pL(useN(1:15:end)),'bo'); axis square; hold on
%     title('head d_Theta, phi L, nonapp & app');
%     x = linspace(-40,40);
%     y = linspace(-80,80);
%     plot(x,y);  xlim([-40 40]); ylim([-80 80]);
%     plot(dT(use(1:15:end)),pL(use(1:15:end)),'og');
%     if vid==(useData(end))
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(3);
%     plot(dT(useN(1:15:end)),tR(useN(1:15:end)),'bo');axis square; hold on
%     title('head d_Theta, theta R, nonapp & app');
%     x = linspace(-40,40);
%     y = linspace(-40,40);  xlim([-40 40]); ylim([-40 40]);
%     plot(-x,y);
%     plot(dT(use(1:15:end)),tR(use(1:15:end)),'og');
%     if vid==(useData(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(4);
%     plot(dT(useN(1:15:end)),tL(useN(1:15:end)),'bo');axis square; hold on;
%     title('head d_Theta, theta L, nonapp & app');
%     x = linspace(-40,40);  xlim([-40 40]); ylim([-80 80]);
%     y = linspace(-80,80);
%     plot(x,y);
%     plot(dT(use(1:15:end)),tL(use(1:15:end)),'og');
%     if vid==(useData(end))
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
%     if vid==(useData(end))
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
%     if vid==(useData(end))
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
%     if vid==(useData(end))
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
%     if vid==(useData(end))
%
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%
% end

%%

% close all
%
% % figure%('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useData)
%     clear dT dpR dpL dtR dtL
%     %     subplot(rownum,colnum,vid);
%     nframe = min(length(dist{useData(vid)}),length(radR{useData(vid)}));
%     nframe = min(nframe, length(radL{useData(vid)}));
%     nframe = min(nframe, length(mouseSp{useData(vid)}));
%     nframe = min(nframe,length(cricketSp{useData(vid)}));
%     nframe = min(nframe,length(az{useData(vid)}));
%
%     r=dist{useData(vid)}(1:nframe);
%     rR=radR{useData(vid)}(1:nframe);%rR = rR - min(rR(:)); rR = rR ./ max(rR(:))
%     lR=radL{useData(vid)}(1:nframe);%lR = lR - min(lR(:)); lR = lR ./ max(lR(:))
%
%     mouseV=mouseSp{useData(vid)}(1:nframe); crSp=cricketSp{useData(vid)}(1:nframe);
%     azC=az{useData(vid)}(1:nframe);
%     useN= appEpoch{vid}==0;
%     use = (appEpoch{vid});
%     figure(1);
%     plot(rR(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on
%     xlabel('r eye'); ylabel('l eye');
%     title('two eyes, rad');
% %     x = linspace(-40,40);
% %     y = linspace(-60,60);
% %     plot(-x,y);
%   %  xlim([0 1]); ylim([0 1]);
%    % plot(rR(use(1:22:end)),lR(use(1:22:end)),'og');
%    % plot(rR(use),lR(use),'og');
%
%
% %     figure(2);
% %     plot(r(useN(1:22:end)),rR(useN(1:22:end)),'bo'); axis square; hold on
% %     title('range and R Rad');
% %     xlabel('range to cricket (cm)'); ylabel('R eye');
% % %     x = linspace(-40,40);
% % %     y = linspace(-80,80);
% % %     plot(x,y);  xlim([-40 40]); ylim([-80 80]);
% % %     figure(3);
% %     plot(r(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on
% %     title('range and L Rad');
% %     xlabel('range to cricket (cm)'); ylabel('L eye');
% % %     x = linspace(-40,40);
% % %     y = linspace(-40,40);  xlim([-40 40]); ylim([-40 40]);
% % %     plot(-x,y);
% %     figure(4);
% %     plot(mouseV(useN(1:22:end)),rR(useN(1:22:end)),'bo');axis square; hold on;
% %     title('mouse speed & R rad');
% %     xlabel('mouse Speed (cm/sec)'); ylabel('R Rad');
% % %     x = linspace(-40,40);  xlim([-40 40]); ylim([-80 80]);
% % %     y = linspace(-80,80);
% % %     plot(x,y);
% %     plot(mouseV(use(1:22:end)),rR(use(1:22:end)),'og');
% %    figure(5);%subplot(rownum,colnum,vid);
% %     plot(mouseV(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on;
% %     title('mouse speed & L rad');
% %     xlabel('mouse Speed (cm/sec)'); ylabel('L Rad');
% % %     x = linspace(-80,80);
% % %     y = linspace(-50,50);
% % %     plot(-x,y); xlim([-80 80]); ylim([-50 50]);
% %     plot(mouseV(use(1:22:end)),lR(use(1:22:end)),'og');
% %       figure(6);
% %     plot(crSp(useN(1:22:end)),rR(useN(1:22:end)),'bo');axis square; hold on;
% % %     x = linspace(-80,80);
% % %     y = linspace(-80,80);
% % %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
% %     title('cricket speed & R rad');
% %      xlabel('cricket Speed'); ylabel('R Rad');
% %     plot(crSp(use(1:22:end)),rR(use(1:22:end)),'og');
% %     if vid==(useData(end))
% %
% %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %     end
% %
% %      figure(7);
% %     plot(crSp(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on;
% % % %     x = linspace(-80,80);
% % % %     y = linspace(-80,80);
% % % %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
% % %   title('cricket speed & L rad');
% %      xlabel('cricket Speed'); ylabel('l Rad');
% % %     plot(crSp(use(1:22:end)),lR(use(1:22:end)),'og');
% % %
% %     figure(8);
% %     plot(azC(useN(1:22:end)),r(useN(1:22:end)),'b.');axis square; hold on;
% % % %     x = linspace(-80,80);
% % % %     y = linspace(-80,80);
% % % %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
% %     title('az and range');
% % %     plot(azC(use(1:22:end)),r(use(1:22:end)),'.g');
% %        if vid==(useData(end))
% %
% %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %     end
%
% end



% figure%('units','normalized','outerposition',[0 0 1 1])
% for vid=1:length(useData)
%     clear dT dpR dpL dtR dtL
%     %     subplot(rownum,colnum,vid);
%     nframe = min(length(dist{useData(vid)}),length(radR{useData(vid)}));
%     nframe = min(nframe, length(radL{useData(vid)}));
%     nframe = min(nframe, length(mouseSp{useData(vid)}));
%     nframe = min(nframe,length(cricketSp{useData(vid)}));
%     nframe = min(nframe,length(az{useData(vid)}));
%
%     r=dist{useData(vid)}(1:nframe);
%     rR=radR{useData(vid)}(1:nframe);%rR = rR - min(rR(:)); rR = rR ./ max(rR(:))
%     lR=radL{useData(vid)}(1:nframe);%lR = lR - min(lR(:)); lR = lR ./ max(lR(:))
%
%     mouseV=mouseSp{useData(vid)}(1:nframe); crSp=cricketSp{useData(vid)}(1:nframe);
%     azC=az{useData(vid)}(1:nframe);
%   %  useN= appEpoch{vid}==0;
%     use = (appEpoch{vid});
%     figure(1);
%    % plot(rR(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on
%     xlabel('r eye'); ylabel('l eye');
%     title('two eyes, rad');
% %     x = linspace(-40,40);
% %     y = linspace(-60,60);
% %     plot(-x,y);
%     xlim([15 46]); ylim([15 46]);
%     plot(rR(use(1:22:end)),lR(use(1:22:end)),'og');
%     if vid==(useData(end))
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%    % figure(2);
% %     plot(r(useN(1:22:end)),rR(useN(1:22:end)),'bo'); axis square; hold on
% %     title('range and R Rad');
% %     xlabel('range to cricket (cm)'); ylabel('R eye');
% % %     x = linspace(-40,40);
% % %     y = linspace(-80,80);
% % %     plot(x,y);  xlim([-40 40]); ylim([-80 80]);
%  %   plot(r(use(1:22:end)),rR(use(1:22:end)),'og');
% %     if vid==(useData(end))
% %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %     end
% %     figure(3);
% % %     plot(r(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on
% % %     title('range and L Rad');
% % %     xlabel('range to cricket (cm)'); ylabel('L eye');
% % % %     x = linspace(-40,40);
% % % %     y = linspace(-40,40);  xlim([-40 40]); ylim([-40 40]);
% % % %     plot(-x,y);
% %      plot(r(use(1:22:end)),lR(use(1:22:end)),'og');
% % %     if vid==(useData(end))
% % %
% % %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% % %     end
% %      figure(4);
% % %     plot(mouseV(useN(1:22:end)),rR(useN(1:22:end)),'bo');axis square; hold on;
% % %     title('mouse speed & R rad');
% % %     xlabel('mouse Speed (cm/sec)'); ylabel('R Rad');
% % % %     x = linspace(-40,40);  xlim([-40 40]); ylim([-80 80]);
% % % %     y = linspace(-80,80);
% % % %     plot(x,y);
% %     plot(mouseV(use(1:22:end)),rR(use(1:22:end)),'og');
% % %     if vid==(useData(end))
% % %
% % %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% % %     end
% %      figure(5);%subplot(rownum,colnum,vid);
% % %     plot(mouseV(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on;
% % %     title('mouse speed & L rad');
% % %     xlabel('mouse Speed (cm/sec)'); ylabel('L Rad');
% % % %     x = linspace(-80,80);
% % % %     y = linspace(-50,50);
% % % %     plot(-x,y); xlim([-80 80]); ylim([-50 50]);
% %      plot(mouseV(use(1:22:end)),lR(use(1:22:end)),'og');
% % %     if vid==(useData(end))
% % %
% % %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% % %     end
% %      figure(6);
% % %     plot(crSp(useN(1:22:end)),rR(useN(1:22:end)),'bo');axis square; hold on;
% % % %     x = linspace(-80,80);
% % % %     y = linspace(-80,80);
% % % %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
% % %     title('cricket speed & R rad');
% % %     xlabel('cricket Speed'); ylabel('R Rad');
% %      plot(crSp(use(1:22:end)),rR(use(1:22:end)),'og');
% % %     if vid==(useData(end))
% % %
% % %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% % %     end
% % %
% %       figure(7);
% % %     plot(crSp(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on;
% % % %     x = linspace(-80,80);
% % % %     y = linspace(-80,80);
% % % %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
% % %   title('cricket speed & L rad');
% % %     xlabel('cricket Speed'); ylabel('l Rad');
% %     plot(crSp(use(1:22:end)),lR(use(1:22:end)),'og');
% % %        if vid==(useData(end))
% % %
% % %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% % %     end
% % %
% %     figure(8);
% % %     plot(azC(useN(1:22:end)),r(useN(1:22:end)),'b.');axis square; hold on;
% % % %     x = linspace(-80,80);
% % % %     y = linspace(-80,80);
% % % %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
% % %     title('az and range');
% %      plot(azC(use(1:22:end)),r(use(1:22:end)),'.g');
% % %        if vid==(useData(end))
% % %
% % %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% % %     end
%
% end




%%
% for vid=1:length(useData)
%
% nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
% nframe = min(nframe, length(az{useData(vid)}));
%
% tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe); azC=(az{useData(vid)}(1:end-4));
% pR=Rphi{useData(vid)}(1:nframe); pL=Lphi{useData(vid)}(1:nframe);
% useN= appEpoch{vid}==0%(1:nframe)%(1:end-4)==0;
% % figure(1)
% % plot(tR(useN),tL(useN),'bo'); axis square; hold on;%xlim([-pi pi]);
% % title('R & L theta');
% figure(2);
% plot(tR(useN),pR(useN),'bo');axis square; hold on;%xlim([-pi pi]);
% title('R theta & phi');
% figure(3);
% plot(tL(useN),pL(useN),'bo');axis square; hold on;%xlim([-pi pi]);
% title('L theta & phi');
% end
%%
% figure
% for vid=1:length(useData)
%
% nframe = min(length(dRtheta{useData(vid)}),length(dLtheta{useData(vid)}));
% nframe = min(nframe, length(az{useData(vid)}));
% nframe = min(nframe, length(appEpoch{vid}));
%
% tR=dRtheta{useData(vid)}(1:nframe); tL=dLtheta{useData(vid)}(1:nframe);
% ht=d_Theta{useData(vid)}(1:nframe); mnEye=.5*(tR+tL);
% azC=(az{useData(vid)})%(1:end-4));
% pR=Rphi{useData(vid)}(1:nframe); pL=Lphi{useData(vid)}(1:nframe);
% useN=appEpoch{vid}(1:nframe)==0;
% use= appEpoch{vid}(1:nframe)==1;
% plot(ht(useN),mnEye(useN),'b.'); hold on; axis square; axis([-frRate frRate -frRate frRate]);
% % plot(ht(use),mnEye(use),'.g')
% % figure(3)
% % plot(tR(use),azC(use),'b.'); axis square; hold on;%xlim([-pi pi]);
% % % plot(tR(useN),azC(useN),'b.'); axis square; hold on;%xlim([-pi pi]);
% % title('R th and az');
% % figure(4)
% % plot(tL(use),azC(use),'b.');axis square; hold on;%xlim([-pi pi]);
% % % plot(tL(useN),azC(useN),'b.');axis square; hold on;
% % title('L th and az');
%
% % figure(1)
% % plot(tR(use),tL(use),'go'); axis square; hold on;%xlim([-pi pi]);
% % title('R & L theta');
% % figure(2);
% % plot(tR(use),pR(use),'go');axis square; hold on;%xlim([-pi pi]);
% % title('R theta & phi');
% % figure(3);
% % plot(tL(use),pL(use),'go');axis square; hold on;%xlim([-pi pi]);
% % title('L theta & phi');
% end
%  if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%
% clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
% figure
% for vid=1:length(useData)
%     clear corrR lagsR corrRA lagsRA corrL corrLA lagsL lagsLA
% nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
% nframe = min(nframe, length(appEpoch{vid}));
% tR=Rtheta{useData(vid)}; tL=Ltheta{useData(vid)};hT=thetaHead{useData(vid)};
% useN= appEpoch{vid}==0; use= appEpoch{vid}==1; pre=find(use)-4%%appT=find(appEpoch{vid}==1);
% useN=find(pre)
% mnEye=(tR+tL)./(length(nframe*2))
% gaze=(hT+mnEye');
% % HGazeR=(hT+tR');HGazeL=(hT+tL');
% % cT=((cricketTh{useData(vid)}));
% cT=((cricketTh{useData(vid)}));
%
% % rGazeA=(cT(use)+tR(use)');lGazeA=(cT(use)+tL(use)')
% subplot(1,2,1)
% plot(gaze(pre),cT(pre));title('gaze, cricket Th');xlim([0 4frRate0]);
% hold on; plot(gaze(use),cT(use),'og');
% subplot(1,2,2)
% plot(hT(pre),mnEye(pre)); hold on; plot(tL);title('head, eye theta');xlim([0 4frRate0]);
% plot(hT(use),mnEye(use),'og')
% % subplot(3,1,3)
% % plot(HGazeR);hold on;plot(HGazeL); title('eye Gaze, head theta');xlim([0 4frRate0]);
% % plot(find(use),HGazeR(use),'og');hold on; plot(find(use),HGazeL(use),'og') %
% % plot(cT);
%
% % if sum(~isnan(cT(useN)))>20 &sum(~isnan(gaze(useN)))>20
% % [corrR lagsR]=nanxcorr(gaze(useN),cT(useN),frRate,'coeff');
% % % [corrL lagsL]=nanxcorr(cT(useN),HGazeL(useN),frRate,'coeff');
% % else
% %       lag=1:(2*(frRate)+1)
% %       corrR=NaN(size([lag],1),'like',lag);
% %     %  corrL=NaN(size([lag],1),'like',lag);
% %       %dbstop
% % end
% %
% % if sum(~isnan(cT(use)))>20 &sum(~isnan(gaze(use)))>20 %sum(use)>10 %& sum(~isnan(cT(use)))<1 & sum(~isnan(HGazeR(use)))<1 & sum(~isnan(HGazeL(use)))<1
% % [corrRA lagsRA]=nanxcorr(gaze(use),cT(use),frRate,'coeff');
% % %[corrLA lagsLA]=nanxcorr(HGazeL(use),cT(use),frRate,'coeff');
% % else
% %   lag=1:(2*(frRate)+1)
% %       corrRA=NaN(size([lag],1),'like',lag);
% %  %     corrLA=NaN(size([lag],1),'like',lag);
% %     end
% %  % if sum(isnan(cT(use)))>20 & sum(isnan(HGazeR(use)))>20 &sum(isnan(cT(useN)))>20 & sum(isnan(HGazeR(useN)))>20
% %        corrRAll(vid,:)=corrR;
% %        corrRAAll(vid,:)=corrRA;
% % %        corrLAll(vid,:)=corrL;
% % %        corrLAAll(vid,:)=corrLA;
% %
% % %     else
% %    %  end
%
% end
% %
% figure('units','normalized','outerposition',[0 0 1 1])
% errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errRA = nanstd(corrRAAll)/(sqrt(length(corrRAAll)));
% % errL= nanstd(corrLAll)/(sqrt(length(corrLAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));
%
% shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
% % shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
% shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
% % shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);
%
% plot([frRate+1,frRate+1],[1,-1],'--','Color', [.5 .5 .5]); %ylim([-.6 .6]);
% xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% %L(2) = plot(nan, nan, 'r-');
% L(2) = plot(nan, nan, 'g-');
% %L(4) = plot(nan, nan, 'c-');

%legend(L,{'R eye Non-app','L eye Non-app','R App','L App'}); title('mean gaze and cricket theta corr');%%
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
clear tNcounts tAcounts
figure
for vid=1:length(useData)
    nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
  %  pR=dRphi{useData(vid)}(1:nframe); pL=dLphi{useData(vid)}(1:nframe);
    tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
    tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
    useN= appEpoch{vid}==0;
    % subplot(1,2,1)
    plot(tR(useN(1:70:end)),tL(useN(1:70:end)),'bo'); hold on; axis square
    %plot(tR(use(1:70:end)),pR(use(1:70:end)),'go'); hold on
    %  subplot(1,2,2);
    %   plot(tL(useN(1:70:end)),pL(useN(1:70:end)),'bo'); hold on; axis square
    clear r l
    if sum(useN)>4
    nBins=-90:20:90
    r= (hist(tR(useN),nBins))/sum(useN);l= (hist(tL(useN),nBins))/sum(useN);
    else
        r=NaN;l=NaN;
    end
   tNcounts(vid,:,1)=r;
    tNcounts(vid,:,2)=l;
end

for vid=1:length(useData)
    nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
    tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
    tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;

    use = (appEpoch{vid});
    plot(tR(use(1:70:end)),tL(use(1:70:end)),'go'); hold on; axis square; %ylim([-frRate frRate]); xlim([-frRate frRate]);
    clear r l
    if sum(use)>4
        nBins=-90:20:90
        r= (hist(tR(use),nBins))/sum(use);l= (hist(tL(use),nBins))/sum(use);
    else
        r=NaN;l=NaN;
    end
    tAcounts(vid,:,1)=r;
    tAcounts(vid,:,2)=l;
end


allThN(:,1)=nansum(tNcounts(:,:,1),1)./(length(useData))
allThN(:,2)=nansum(tNcounts(:,:,2),1)./(length(useData))

allThA(:,1)=nansum(tAcounts(:,:,1),1)./(length(useData))
allThA(:,2)=nansum(tAcounts(:,:,2),1)./(length(useData))

figure;subplot(1,2,1); plot(nBins,allThN(:,1),'b'); hold on; plot(nBins,allThA(:,1),'g'); title('Right Eye Theta'); axis square
ylim([0 .6]);
subplot(1,2,2);plot(nBins,allThN(:,2),'b'); hold on; plot(nBins,allThA(:,2),'g');title('Left Eye Theta');axis square
ylim([0 .6]);
legend('non approach','approach')


%
clear tNcounts tAcounts allThN allThA
figure
for vid=1:length(useData)
    nframe = min(length(Rphi{useData(vid)}),length(Lphi{useData(vid)}));
  %  pR=Rphi{useData(vid)}(1:nframe); pL=Lphi{useData(vid)}(1:nframe);
    tR=Rphi{useData(vid)}(1:nframe); tL=Lphi{useData(vid)}(1:nframe);
    tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
    useN= appEpoch{vid}==0;
    % subplot(1,2,1)
    plot(tR(useN(1:70:end)),tL(useN(1:70:end)),'bo'); hold on; axis square
    %plot(tR(use(1:70:end)),pR(use(1:70:end)),'go'); hold on
    %  subplot(1,2,2);
    %   plot(tL(useN(1:70:end)),pL(useN(1:70:end)),'bo'); hold on; axis square
    clear r l
    if sum(useN)>4
    nBins=-50:5:50
    r= (hist(tR(useN),nBins))/sum(useN);l= (hist(tL(useN),nBins))/sum(useN);
    else
        r=NaN;l=NaN;
    end
   tNcounts(vid,:,1)=r;
    tNcounts(vid,:,2)=l;
end

for vid=1:length(useData)
    nframe = min(length(Rphi{useData(vid)}),length(Lphi{useData(vid)}));
    tR=Rphi{useData(vid)}(1:nframe); tL=Lphi{useData(vid)}(1:nframe);
   tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;

    use = (appEpoch{vid});
    plot(tR(use(1:70:end)),tL(use(1:70:end)),'go'); hold on; axis square; %ylim([-frRate frRate]); xlim([-frRate frRate]);
    clear r l
    if sum(use)>4
        nBins=-50:5:50
        r= (hist(tR(use),nBins))/sum(use);l= (hist(tL(use),nBins))/sum(use);
    else
        r=NaN;l=NaN;
    end
    tAcounts(vid,:,1)=r;
    tAcounts(vid,:,2)=l;
end


allThN(:,1)=nansum(tNcounts(:,:,1),1)./(length(useData))
allThN(:,2)=nansum(tNcounts(:,:,2),1)./(length(useData))

allThA(:,1)=nansum(tAcounts(:,:,1),1)./(length(useData))
allThA(:,2)=nansum(tAcounts(:,:,2),1)./(length(useData))

figure;subplot(1,2,1); plot(nBins,allThN(:,1),'b'); hold on; plot(nBins,allThA(:,1),'g'); title('Right Eye Phi'); axis square
ylim([0 .32]);
subplot(1,2,2);plot(nBins,allThN(:,2),'b'); hold on; plot(nBins,allThA(:,2),'g');title('Left Eye Phi');axis square
ylim([0 .32]);
legend('non approach','approach')




%%




if savePDF
    pSname='T:\PreyCaptureAnalysis\Data\';
    filen=sprintf('%s','DEINTERLACED_Analyzed_AllAnimals_052320_halfShift','.pdf')
    pdfilename=fullfile(pSname,filen);
    dos(['ps2pdf ' psfilename ' ' pdfilename]);
    delete(psfilename);
else
    pSname='T:\PreyCaptureAnalysis\Data\';
end


afilename=sprintf('%s','DEINTERLACED_Analyzed_AllAnimals_052320_halfShift','.mat')
save(fullfile(pSname, afilename))

