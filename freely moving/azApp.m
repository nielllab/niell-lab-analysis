close all 
clear all

% load('ACCAnalyzed_AllAnimals_121819_a.mat')
  load('ACCAnalyzed_AllAnimals_010820_noDLS.mat')
% load('ACC_deInter_Analyzed_AllAnimals_011520_a.mat')

clear L H R
skip = 0; %%% only shows figures at this interval
nthresh  = 60;  %%% threshold for number of approach points to be included in averages
set(groot,'defaultFigureVisible','on') %disable figure plotting

savePDF=0; dbstop if error
if savePDF
    psfilename = 'C:\analysisPS.ps';
    if exist(psfilename,'file')==2; delete(psfilename);end
end

vgDiff = [];
hdDiff = [];
gzDiff = [];
vergAll = [];
allHth = [];
allGz = [];
allVg = [];
gzApp = [];
reset=[];
stable=[];
headTurn=[];
phiVg=[];

hthAll=[];
dthAll = []; diffThAll = [];
gazeAll = [];
dgazeAll =[]
vergAll = [];
dvergAll=[];
lthAll = []; rthAll=[];

dHeadThApp = [];
dgazeApp=[];
d_vgApp=[];
mnEyeAll=[];
dEyeAll=[];

tiltAll = [];
rollAll = [];
acc_dthAll = [];

appAll=[];

maxlag = 30;
thbins = -60:5:60;

for i = 67%1:length(appEpoch)
    
    vid = useData(i)
    
    %%% get approaches
    app = appEpoch{i};
    appAll=[appAll app];
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
    hthnonan = hth;
    hthnonan(abs(diff(hth))>90)=NaN;
    
    % head, gaze, vergence not separated by approach/non-approach
    
    hthAll=[hthAll hth(1:end-1)];
    dthAll = [dthAll dth(1:end-1)'];
    
    diffThAll = [ diffThAll diff(hth)];
    gazeAll = [gazeAll gaze(1:end-1)];
    dgazeAll =[dgazeAll diff(gaze)];
    vergAll = [vergAll vergence(1:end-1)'];
    dvergAll=[dvergAll diff(vergence')];
    mnEyeAll=[mnEyeAll mnEyeTh(1:end-1)'];
    dEyeAll=[dEyeAll diff(mnEyeTh')];
    lthAll = [lthAll lth(1:end-1)'];
    rthAll = [rthAll rth(1:end-1)'];
    
    if exist('accelData','var')
        tiltAll = [tiltAll tilt(1:end-1)'];
        rollAll = [rollAll roll(1:end-1)'];
        acc_dthAll = [acc_dthAll acc_dth(1:end-1)'];
    end
    
    %     he=[diff(headTh);diff(mnEye)];
    %     % gm = fitgmdist(he',6);
    %     gm = fitgmdist(he',5);
    %     idx = cluster(gm,he')
    %     X=he;
    %
    %     figure;
    %     gscatter(X(1,:),X(2,:),idx);
    %     legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Location','best');
    %     if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    %
    %     resetPt=((diff(mnEye))>5 & (diff(fullgaze))>5) | ((diff(mnEye))<-5 & (diff(fullgaze))<-5);
    %     stablePt=idx==1|idx==2
    %     stablePt = ((diff(fullgaze))>-5 & (diff(fullgaze))<5) & ((diff(mnEye)<30) & ((diff(mnEye)>-30)));
    %     headTurnPt =((diff(mnEye)>-2 & diff(mnEye)<2) & ((diff(fullgaze)>10 |(diff(fullgaze)<-10))) & (diff(headTh)>10 | diff(headTh)<-10));
    %     resetPt=((diff(mnEye))>5 & (diff(fullgaze))>5) | ((diff(mnEye))<-5 & (diff(fullgaze))<-5);
    %
    %
    %     reset = [reset resetPt];
    %     stable = [stable stablePt];
    %     headTurn = [headTurn headTurnPt];
    
    
    
    
    %     dHeadThApp = [dHeadThApp diff(headTh(app(1:end-1)))];
    %     dgazeApp=[dgazeApp diff(fullgaze(app(1:end-1)))];
    %     d_vgApp=[d_vgApp diff(fullvg(app(1:end-1))')];
    
    
    hthApp = hth(appRange)-nanmedian(hth(mainApp));
    hthApp = mod(hthApp + 180,360)-180;
    
    gzApp = hthApp +0.5*( rth(appRange) +lth(appRange))';
    mnEyeApp =  0.5*(rth(appRange) +lth(appRange))';
    
    
    %%% calculate change in position at different lags, as measure of stability
    
    for lag = 1:20;
        d_gz(i,lag) = nanmean(abs(gaze(1:end-lag) - gaze((lag+1):end)));
        d_hd(i,lag) = nanmean(abs(hth(1:end-lag) - hth((lag+1):end)));
    end
    
    
    %%% draw figures
    
    if round(i/skip)==i/skip
        
        %%% eye/head correlations
        
        figure
        subplot(2,3,1)
        plot(dth(app),dlth(app),'.'); axis square; axis([-15 15 -15 15]); hold on; plot([-15 15],[15 -15])
        title(sprintf('ccL = %0.2f scale = %0.2f',ccL(vid),scaleL(vid))); xlabel('dHead th'); ylabel('dLeft theta');
        subplot(2,3,2)
        plot(dth(app),drth(app),'.'); axis square; axis([-15 15 -15 15]); hold on; plot([-15 15],[15 -15])
        title(sprintf('ccR = %0.2f scale = %0.2f',ccR(vid),scaleR(vid))); xlabel('dHead th'); ylabel('dRight theta');
        subplot(2,3,4)
        plot(dth(~app),dlth(~app),'.'); axis square; axis([-15 15 -15 15]); hold on; plot([-15 15],[15 -15])
        subplot(2,3,5)
        plot(dth(~app),drth(~app),'.'); axis square; axis([-15 15 -15 15]); hold on; plot([-15 15],[15 -15])
        
        subplot(2,3,3);
        plot(dth(app),0.5*(drth(app)+dlth(app)),'.');  axis square; axis([-15 15 -15 15]); hold on; plot([-15 15],[15 -15])
        title(sprintf('vid %d',vid));  xlabel('dHead th'); ylabel('dLeft+dRight theta');
        
        %         subplot(2,3,6);
        %         plot(dth(~app),0.5*(drth(~app)+dlth(~app)),'.');  axis square; axis([-15 15 -15 15]); hold on; plot([-15 15],[15 -15])
        %
        subplot(2,3,6);
        plot(d_gz(i,:));
        hold on
        plot(d_hd(i,:));
        title('gaze drift, head drift');
        axis([1 size(d_hd,2) 0 20])
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
        %%% pursuit eye/head/gaze positions
        
        figure
        subplot(6,1,1);
%         plot(hthnonan,'k'); hold on;
        plot(rth,'r');hold on; plot(lth,'b'); legend('right th','left th');
        plot(find(app),ones(sum(app),1)*90,'g.'); %ylim([-180 180])
        title(sprintf('vid %d',vid));
        
        roll = accelChannels{vid}(:,1); roll=roll-nanmean(roll); roll=medfilt1(roll,5);
        tilt=accelChannels{vid}(:,2); tilt=tilt-nanmean(tilt); tilt=medfilt1(tilt,8);
        g3=accelChannels{vid}(:,6); g3=g3-nanmean(g3);
        
        subplot(6,1,2)
        plot(roll); hold on;
        plot(tilt);
        plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        legend('roll','tilt');

        subplot(6,1,3);
        plot(0.5*( rth(appRange) +lth(appRange))','k','LineWidth',2); hold on; plot( rth(appRange)','r'); hold on; plot(lth(appRange)','b');
        ylim([-30 30]);  xlim([0 max(length(appRange),1)]);
        plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        legend('mean','right','left');
        
        subplot(6,1,4)
        plot(roll(appRange)); hold on;
        plot(tilt(appRange));
        plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        plot([1  max(length(appRange),1)],[0 0],'--')
        xlim([0 max(length(appRange),1)]);
        if ~isempty(appRange)
        ylim([((min(min(roll(appRange),tilt(appRange))))-5) ((max(max(roll(appRange),tilt(appRange))))+5)]);
        end
        legend('roll','tilt');
        
        subplot(6,1,5)
        vgPhi=rphi(appRange)-lphi(appRange);
        plot(vgPhi); hold on
        plot(roll(appRange),'k'); hold on
        plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        xlim([0 max(length(appRange),1)]);
        legend('R phi- L phi','acc roll')
            if ~isempty(appRange)
        ylim([((min(min(roll(appRange),vgPhi)))-5) ((max(max(roll(appRange),vgPhi)))+5)]);
        end
        
        gzApp = hthApp +0.5*( rth(appRange) +lth(appRange))';
        subplot(6,1,6);
        hold on; plot(hthApp +0.5*( rth(appRange) +lth(appRange))','k','LineWidth',2);
        plot(hthApp,'Color',[0 0.75 0],'LineWidth',2);
        if appOffset>0
            plot([appOffset appOffset],[-60 60],'g');end
        plot([endOffset endOffset],[-60 60],'r');
        if ~isempty(min(hthApp))
            ylim([min(hthApp)-20 max(hthApp)+20]);
        else  ylim([-60 60]);end
        xlim([0 max(length(appRange),1)]); xlabel('frames'); ylabel('deg');
        legend('gaze','head')
        
       if exist('gm','var')
           X = [dth(appRange)'; [diff(mnEyeApp) 0]; mnEyeApp]';
        idx = cluster(gm,X);
        sacc = find(idx==3);
        for i = 1:length(sacc)-1;
            plot(sacc(i):sacc(i)+1,gzApp(sacc(i):sacc(i)+1),'r')
        end
       end
       
        
        %         appGaze =(hthApp +0.5*( rth(appRange) +lth(appRange))')
        %         hold on;plot(find(resetPt(appRange)),appGaze(resetPt(appRange)),'ob');
        %         plot(find(stablePt(appRange)),appGaze(stablePt(appRange)),'og')
        %         plot(find(headTurnPt(appRange)),appGaze(headTurnPt(appRange)),'or')
        %
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
        
        %         legend('gaze','non-comp','eyes still');
        %   plot(hthApp(use)+mnEye(use)','c','LineWidth',2); %gaze
        %         plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        % ylim([-40 40]); xlim([0 max(length(appRange),1)]);
        %         xlabel('frames'); ylabel('gaze(deg)');
        %         set(gcf,'Position',[440 100 560 640])
        
        
        
    end
    
    drawnow
    
    %%% eye movements around the time of head saccades
    %%% work in progress!
    
    win = 3;
    thSum = conv(dth,ones(win,1));
    thrange = [2.5 5 10 20 40 80];
    for rep = 1:(length(thrange)-1)
        %         if rep==1  %%% positive and negative saccades
        %              saccEp = thSum>30;
        %         else
        %             saccEp = thSum<-30;
        %         end
        saccEp = thSum>thrange(rep) & thSum<thrange(rep+1);  %%% different size ranges
        
        sacc = find(diff(saccEp)>0);
        %sacc = find(saccEp);
        %%% maybe also find saccades by dtheta(n-5:n-1)<3 & dtheta(n)>15
        sacc= sacc(sacc>31 & sacc<length(dth)-30);
        range = -30:30;
        clear Hsacc dthSacc Lsacc Rsacc  AZsacc
        baserange = -5:-1;
        for j = 1:length(sacc);
            Lsacc(:,j) = lth(sacc(j)+range);% - nanmean(lth(sacc(j)+baserange));
            Rsacc(:,j) = rth(sacc(j)+range);%- nanmean(rth(sacc(j)+baserange));
            this_dth =dth(sacc(j)+range);
            this_dth(isnan(this_dth))=0;
            hd = cumsum(this_dth);
            Hsacc(:,j) = hd - nanmean(hd(25:29));
            dthSacc(:,j) = thSum(sacc(j)+range);
            AZsacc(:,j) = azdeg(sacc(j)+range);
        end
        
        
        if length(sacc)==0
            Lsacc=NaN;
            Rsacc= NaN;
            Hsacc=NaN;
            AZsacc= NaN;
            appSacc=NaN;
        else
            endAZ = nanmean(AZsacc(40:50,:),1);
            appSacc = abs(endAZ)<45;
            appSacc = max(abs(Hsacc(20:29,:)),[],1)<20;% & abs(endAZ)<45;
        end
        
        if sum(appSacc)>0 & sum(~appSacc)>0
            
            if round(i/skip)==i/skip+1
                %                 figure
                %                 subplot(2,2,1)
                %                 plot(range,Hsacc(:,~appSacc),'r');hold on;plot(range,nanmean(Hsacc(:,~appSacc),2),'r','LineWidth',4)
                %                 plot(range,Hsacc(:,appSacc),'b'); ylim([-90 90]); hold on; plot(range,nanmean(Hsacc(:,appSacc),2),'b','LineWidth',4)
                %                 title('head')
                %
                %                     subplot(2,2,2)
                %                     plot(range,dthSacc);ylim([-90 90]); hold on; plot(range,nanmean(dthSacc,2),'g','LineWidth',4)
                %                     title('dth')
                
                subplot(2,2,2);
                plot(range,nanmean(Hsacc(:,appSacc),2),'k'); hold on
                plot(range,nanmean(Lsacc(:,appSacc),2),'b');
                plot(range,nanmean(Rsacc(:,appSacc),2),'g');
                ylim([-90 90]); xlim([-20,20])
                title(sprintf('rep %d',rep))
                
                subplot(2,2,3)
                plot(range,Lsacc(:,~appSacc),'r'); hold on; plot(range,Lsacc(:,appSacc),'b');
                plot(range,nanmean(Lsacc(:,~appSacc),2),'r','linewidth',4);
                plot(range,nanmean(Lsacc(:,appSacc),2),'b','linewidth',4);
                ylim([-45 45])
                title('left')
                
                subplot(2,2,4)
                
                plot(range,Rsacc(:,~appSacc),'r'); hold on; plot(range,Rsacc(:,appSacc),'b');
                plot(range,nanmean(Rsacc(:,~appSacc),2),'r','linewidth',4); plot(range,nanmean(Rsacc(:,appSacc),2),'b','linewidth',4);
                ylim([-45 45])
                title('right')
                
                if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
                
            end
            
            
            H(:,rep,1,i) = nanmean(Hsacc(:,appSacc),2);
            R(:,rep,1,i) =nanmean(Rsacc(:,appSacc),2);
            L(:,rep,1,i)=nanmean(Lsacc(:,appSacc),2);
            
            H(:,rep,2,i) = nanmean(Hsacc(:,~appSacc),2);
            R(:,rep,2,i) =nanmean(Rsacc(:,~appSacc),2);
            L(:,rep,2,i)=nanmean(Lsacc(:,~appSacc),2);
        else
            H(:,rep,:,i) = NaN;
            R(:,rep,:,i)=NaN;
            L(:,rep,:,i)=NaN;
        end
        
        
    end
    %
    %         figure
    %     for j = 1:(min(24,length(sacc)))
    %         subplot(4,6,j);
    %         plot(Hsacc(:,j)); hold on;
    %         plot(Lsacc(:,j)); plot(Rsacc(:,j));
    %         ylim([-90 90])
    %     end
    
end


% figure
% subplot(2,3,1)
% plot(thbins,nanmean(lthHist,3)); ylim([0 0.5]); legend('app','non-app');
%
% subplot(2,3,2)
% plot(thbins,nanmean(rthHist,3)); ylim([0 0.5]); legend('app','non-app');
%
% subplot(2,3,4)
% plot(-maxlag:maxlag,nanmean(lthAz_xc,2));  ylim([-0.5 0.5])
%
% subplot(2,3,5)
% plot(-maxlag:maxlag,nanmean(rthAz_xc,2));  ylim([-0.5 0.5])
%
% subplot(2,3,6)
% plot(-maxlag:maxlag,(nanmean(dheadAz_xc,3))); ylim([-0.5 0.5])
%
% subplot(2,3,3);
% plot(thbins, nanmean(az_hist,3),'b')
% hold on
% plot(thbins, nanmean(azthL_hist,3),'r')
% plot(thbins, nanmean(azthR_hist,3),'g')

for rep =1:5
    clear mnEye
    figure
    subplot(2,2,1);
    plot(range,squeeze(H(:,rep,1,:))); axis([-20 20 -90 90])
    hold on;  plot(range,nanmean(squeeze(H(:,rep,1,:)),2),'k');
    subplot(2,2,3);
    plot(range,squeeze(L(:,rep,1,:))); axis([-20 20 -45 45])
    subplot(2,2,4);
    plot(range,squeeze(R(:,rep,1,:))); axis([-20 20 -45 45])
    
    mnEye=(squeeze(nanmean(L(:,rep,1,:),4))+squeeze(nanmean(R(:,rep,1,:),4)))*.5
    subplot(2,2,2);
    % plot(range,squeeze(nanmean(H(:,rep,1,:),4)),'k'); hold on
    plot(range,mnEye,'k'); hold on;
    plot(range,squeeze(nanmean(L(:,rep,1,:),4)),'b'); hold on
    plot(range,squeeze(nanmean(R(:,rep,1,:),4)),'r'); hold on
    % plot(range,squeeze(nanmean(H(:,rep,2,:),4)),'k:'); hold on
    % plot(range,squeeze(nanmean(L(:,rep,2,:),4)),'b:'); hold on
    % plot(range,squeeze(nanmean(R(:,rep,2,:),4)),'r:'); hold on
    axis([-30 30 -10 10])
    title(sprintf('rep %d',rep))
    if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end
%%

%%% clean up a couple values
diffThAll(diffThAll<-180) = diffThAll(diffThAll<-180)+360;
diffThAll(diffThAll>180) = diffThAll(diffThAll>180)-360;
dgazeAll = dEyeAll + dthAll;
appAll = appAll==1;

%%% large head movement seems to be DLC errors according to accelerometers
dthAll(abs(dthAll)>25)=NaN;
dgazeAll = dEyeAll + dthAll;
dgazeAll = dEyeAll + acc_dthAll;


figure
bins = -90:3:90;
h = hist(tiltAll(~appAll),bins)/sum(~appAll);
plot(bins,h,'r');
hold on
h = hist(tiltAll(appAll),bins)/sum(appAll);
plot(bins,h,'g');
xlabel('tilt') 

figure
bins = -90:3:90;
h = hist(rollAll(~appAll),bins)/sum(~appAll);
plot(bins,h,'r');
hold on
h = hist(rollAll(appAll),bins)/sum(appAll);
plot(bins,h,'g');
xlabel('roll')

figure
bins = -20:20;
h = hist(acc_dthAll(~appAll),bins)/sum(~appAll);
plot(bins,h,'r');
hold on
h = hist(acc_dthAll(appAll),bins)/sum(appAll);
plot(bins,h,'g');
xlabel('gyro 3')

figure
bins = -60:2:60;
h = hist(vergAll(~appAll),bins)/sum(~appAll);
plot(bins,h,'r');
hold on
h = hist(vergAll(appAll),bins)/sum(appAll);
plot(bins,h,'g');
xlabel('vergence')


figure
plot(dthAll(appAll),acc_dthAll(appAll),'.');
axis square; axis([-25 25 -25 25]);
xlabel('DLC dtheta'); ylabel('acc dtheta')

figure
plot(acc_dthAll(appAll),dEyeAll(appAll),'.')
axis square; axis([-25 25 -25 25])
xlabel('acc dtheta'); ylabel('dEye theta')

figure
plot(dthAll(appAll),dEyeAll(appAll),'.')
axis square; axis([-25 25 -25 25])
xlabel('DLC dtheta'); ylabel('dEye theta')


figure
plot(-30:30, nanxcorr(acc_dthAll(appAll),dEyeAll(appAll),30,'coeff'))
title('acc dtheta vs eye dtheta')

figure
plot(acc_dthAll(appAll),dgazeAll(appAll),'.');
axis equal;  hold on; plot([-10 10], [-10 10],'r')
xlabel('acc dtheta'); ylabel('gaze dtheta')

figure
plot(dEyeAll(appAll),dgazeAll(appAll),'.');
axis equal;  hold on; plot([-10 10], [-10 10],'r')
xlabel('eye dtheta'); ylabel('gaze dtheta')

figure
plot(mnEyeAll(appAll),dgazeAll(appAll),'.');
axis square; axis([-25 25 -25 25]); hold on;
xlabel('eye position'); ylabel('gaze dtheta')

figure
plot(mnEyeAll(appAll),dEyeAll(appAll),'.');
axis square; axis([-25 25 -25 25]); hold on; 

figure
plot(mnEyeAll(appAll),dthAll(appAll),'.');


figure
plot(max(lthAll(appAll), rthAll(appAll)),dgazeAll(appAll),'.');


he=[acc_dthAll(appAll);dEyeAll(appAll); mnEyeAll(appAll)];
% he=[dthAll(appAll);dgazeAll(appAll);mnEyeAll(appAll)];


% gm = fitgmdist(he',6);
gm = fitgmdist(he',3,'Replicates',10);
% he=[head+eyes;head - eyes];
idx = cluster(gm,he');
X=he;
% figure;
% gscatter(X(1,:),X(2,:),idx); axis equal; %axis([-25 25 -25 25]);
% figure;
% gscatter(X(3,:),X(2,:),idx); axis equal; %axis([-25 25 -25 25]);

figure
gscatter(mnEyeAll(appAll),dgazeAll(appAll),idx); axis equal

figure
gscatter(dEyeAll(appAll),dgazeAll(appAll),idx); axis equal

figure
gscatter(acc_dthAll(appAll),dgazeAll(appAll),idx); axis equal

figure
gscatter(acc_dthAll(appAll),dEyeAll(appAll),idx); axis equal

figure;
gscatter(X(1,:),X(3,:),idx); axis equal; %axis([-25 25 -25 25]);

figure
hist(acc_dthAll(appAll),[-20:20]); xlim([-20 20]);

figure
hist(dgazeAll(appAll),[-20:20]); xlim([-20 20]);





hbins = -20:0.5:20;
figure
plot(hbins,(hist(hdDiff +vgDiff,hbins)/length(vgDiff)))
hold on
plot(hbins,(hist(hdDiff,hbins)/length(hdDiff)))
legend('delta gaze','delta head')
xlabel('deg')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%
figure;
plot(dthAll,dvergAll,'.')
hold on; plot(dthAll(find(appAll)),dvergAll(find(appAll)),'.g')
xlabel('d head yaw'); ylabel('d vergence'); axis equal;axis([-40 40 -40 40])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

p=dthAll>0 & dEyeAll<0;

figure
plot(dthAll,dEyeAll,'.'); hold on
plot(dthAll(find(p)),dEyeAll(find(p)),'.'); hold on
plot(dthAll(find(appAll)),dEyeAll(find(appAll)),'.g');
xlabel('d head yaw'); ylabel('d mn eye th'); axis equal; axis([-40 40 -40 40])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end



figure
plot(dgazeAll,dvergAll,'.'); hold on
plot(dgazeAll(find(appAll)),dvergAll(find(appAll)),'.g');
xlabel('d gaze'); ylabel('d vergence'); axis equal; axis([-40 40 -40 40])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

figure
plot(dgazeAll,dEyeAll,'.'); hold on
plot(dgazeAll(find(appAll)),dEyeAll(find(appAll)),'.g');
xlabel('d gaze');  ylabel('d mn eye th'); axis equal; axis([-40 40 -40 40])
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

figure
plot(dEyeAll+dthAll,dEyeAll-dthAll,'.'); hold on; axis equal; axis([-90 90 -90 90])
plot(dEyeAll(find(appAll)) +dthAll(find(appAll)) ,dEyeAll(find(appAll))-dthAll(find(appAll)),'.g');
xlabel('eye + head');  ylabel('eye - head');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%


% figure
% plot(hdDiff,gzDiff,'.');
% axis equal; axis([-25 25 -25 25]);
% hold on; plot([-25 25],[-25 25],'r')
% xlabel('delta head'); ylabel('delta gaze');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

% figure
% plot(hdDiff,vgDiff,'.');
% axis equal; hold on
% % axis([-20 20 -20 20]);
% plot([-20 20],[20 -20],'r')
% xlabel('delta head'); ylabel('delta theta mean eye');
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
dvg=dvergAll;
dgz=dgazeAll;
dht =dthAll;


r=reset; s=stable; ht=headTurn;
figure
plot(dht,dEyeAll,'.'); hold on; axis square
plot(dht(find(r)),dEyeAll(find(r)),'.c')
plot(dht(find(s)),dEyeAll(find(s)),'.g')
plot(dht(find(ht)),dEyeAll(find(ht)),'.r')

figure
plot(dgz,dEyeAll,'.'); hold on; axis square
plot(dgz(find(r)),dEyeAll(find(r)),'.c')
plot(dgz(find(s)),dEyeAll(find(s)),'.g')
plot(dgz(find(ht)),dEyeAll(find(ht)),'.r')
% 
% figure;plot(randsample(dgz,20000),randsample(dvg,20000),'.'); hold on;
% axis equal; axis([-40 40 -40 40])
% plot([-90 90],[0 0]);
% plot([-0 0],[-50 50])
% plot([-100 100],[-50 50])
% 
% plot(dgz(find(r)),dvg(find(r)),'.c')
% plot(dgz(find(s)),dvg(find(s)),'.g')
% plot(dgz(find(ht)),dvg(find(ht)),'.r')
% axis equal; axis([-40 40 -40 40])
% xlabel('d gaze');ylabel ('d eye theta')
% title('non-approach')
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% 
% figure;plot(randsample(dgz(find(appAll)),20000),randsample(dvg(find(appAll)),20000),'.'); hold on;
% plot(dgz(find(appAll&r)),dvg(find(appAll&r)),'.c')
% plot(dgz(find(appAll&s)),dvg(find(appAll&s)),'.g')
% plot(dgz(find(appAll&ht)),dvg(find(appAll&ht)),'.r')
% axis equal; axis([-40 40 -40 40])
% plot([-90 90],[0 0]);
% plot([-0 0],[-50 50])
% plot([-100 100],[-50 50])
% xlabel('d gaze');ylabel ('d eye theta');
% title('approaches')
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%

% pts = randsample(1:length(dht),15000);
% head = dht(pts); eyes = dEyeAll(pts);


% mnEyeAll dthAll dEyeAll

% he=[dthAll(appAll);dEyeAll(appAll)];
% % he=[dthAll(appAll);dgazeAll(appAll);mnEyeAll(appAll)];
% 
% 
% % gm = fitgmdist(he',6);
% gm = fitgmdist(he',3);
% % he=[head+eyes;head - eyes];
% idx = cluster(gm,he')
% X=he;
% 
% figure;
% gscatter(X(1,:),X(2,:),idx); axis equal; %axis([-25 25 -25 25]);
% figure;
% gscatter(X(2,:),X(3,:),idx); axis equal; %axis([-25 25 -25 25]);
% 
% figure;
% gscatter(X(1,:),X(3,:),idx); axis equal; %axis([-25 25 -25 25]);
% 
% legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Location','best');
%  if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


%%
% he=[dht;dEye];
% [g,C]=kmeans(he',6);
%
% X=he;
% figure;
% plot(X(1,g==1),X(2,g==1),'r.','MarkerSize',12)
% hold on
% plot(X(1,g==2),X(2,g==2),'b.','MarkerSize',12)
% plot(X(1,g==3),X(2,g==3),'c.','MarkerSize',12)
% plot(X(1,g==4),X(2,g==4),'g.','MarkerSize',12)
% plot(X(1,g==5),X(2,g==5),'y.','MarkerSize',12)
% plot(X(1,g==6),X(2,g==6),'m.','MarkerSize',12)
%
% plot(C(1,:),C(2,:),'kx','MarkerSize',15,'LineWidth',3)
% plot(C(3,:),C(4,:),'kx','MarkerSize',15,'LineWidth',3)
%%
% figure;
% subplot(1,2,1)
% plot(randsample(dht,15000),randsample(dEyeAll,15000),'.'); hold on;
% plot(dht(find(r)),dvg(find(r)),'.c')
% plot(dht(find(s)),dvg(find(s)),'.g')
% plot(dht(find(ht)),dvg(find(ht)),'.r')
% axis equal; axis([-40 40 -40 40])
% % plot([-90 90],[0 0]);
% % plot([-0 0],[-50 50])
% % plot([-100 100],[-50 50])
% xlabel('d head th');ylabel ('d eye theta')
% title('non-approach')
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

% for approaches
% figure;
% subplot(1,2,2);
% plot(randsample(dht(find(appAll)),15000),randsample(dEyeAll(find(appAll)),15000),'.'); hold on;
% % plot(dht(find(appT&r)),dvg(find(appT&r)),'.c')
% % plot(dht(find(appT&s)),dvg(find(appT&s)),'.g')
% % plot(dht(find(appT&ht)),dvg(find(appT&ht)),'.r')
% axis equal; axis([-40 40 -40 40])
% % plot([-90 90],[0 0]);
% % plot([-0 0],[-50 50])
% % plot([-100 100],[-50 50])
% xlabel('d head th');ylabel ('d eye theta');
% title('approaches')
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%
figure;plot3(dgz,dvg,dht,'o'); hold on
plot3(dgz(find(r)),dvg(find(r)),dht(find(r)),'oc');
plot3(dgz(find(ht)),dvg(find(ht)),dht(find(ht)),'or')
plot3(dgz(find(s)),dvg(find(s)),dht(find(s)),'og');
axis([-90 90 -40 40 -90 90])
xlabel('d gaze'); ylabel('d eye'); zlabel('d head')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

% figure;plot(dgz,dht,'o'); axis equal; hold on; axis([-90 90 -40 40])
% figure;subplot(1,2,1); plot(dgz,dvg,'o'); hold on; axis equal;axis([-90 90 -40 40])
%
% subplot(1,2,2);plot(dht,dvg,'o'); axis equal;  axis([-40 40 -90 90])
%
% hold on; plot(dvg,dgz,'o')


%%

% figure
% subplot(1,2,1)
% plot(diff(allGz),allVg(1:end-1),'o'); hold on;
% axis square; axis([-90 90 -90 90])
% xlabel('delta gaze'); ylabel('mean eye theta')
% plot(gzDiff,vergAll,'go');
% title('delta gaze, eye theta')
%
% subplot(1,2,2)
% % plot(diff(allGz),diff(allVg),'o'); hold on;
% plot(dHeadTh,dvergence,'o'); hold on;
%
% % plot(diff(allGz(r)'),diff(allVg(r)'),'co')
% % plot(diff(allGz(s)),diff(allVg(s)),'ro')
% % plot(diff(allGz(ht)),diff(allVg(ht)),'go')
%
% axis equal; axis([-90 90 -40 40])
% xlabel('delta gaze'); ylabel('delta mn eye theta')
% plot(gzDiff,vgDiff,'og');
% title('delta gaze, delta eye theta')
%
%%
% clear gaze
% gaze= diff(allGz); verg=diff(allVg)
% figure
% subplot(3,2,1)
% scatplot(gaze(1:10:end),verg(1:10:end),'circles',5,[],[],1,4); %smaller radius elongates along horizontal axis
% axis xy
% xlabel('delta gaze'); ylabel('delta eye theta')
% axis equal;
% axis([-40 40 -40 40]);
% title('rad = 10')
%
% subplot(3,2,2)
% scatplot(gzDiff(1:10:end),vgDiff(1:10:end),'circles',5,[],[],1,4);
% title('Approach: delta gaze, delta eye theta')
% axis equal; axis([-90 90 -90 90])
% xlabel('delta gaze'); ylabel('delta eye theta')
%
% subplot(3,2,3)
% scatplot(gaze(1:10:end),verg(1:10:end),'circles',10,30,[],1,4);
% axis xy
% xlabel('delta gaze'); ylabel('delta eye theta')
% axis equal; axis([-90 90 -90 90]);
% title('rad = 10, mesh=30')
%
%
% subplot(3,2,4)
% scatplot(gzDiff(1:10:end),vgDiff(1:10:end),'circles',10,30,[],1,4);
% title('Approach: delta gaze, delta eye theta')
% axis equal; axis([-90 90 -90 90])
% xlabel('delta gaze'); ylabel('delta eye theta')
%
%
% subplot(3,2,5)
% scatplot(diff(gaze),diff(verg),'circles',10,30,10,1,4);
% axis xy
% xlabel('delta gaze'); ylabel('delta eye theta')
% axis equal; axis([-90 90 -90 90]);
% title('rad = 10, mesh=30 r_mean =10')
%
%
% subplot(3,2,6)
% scatplot(gzDiff(1:10:end),vgDiff(1:10:end),'circles',10,30,10,1,4);
% title('Approach: delta gaze, delta eye theta')
% axis equal; axis([-90 90 -90 90])
% xlabel('delta gaze'); ylabel('delta eye theta')

%%

figure
plot(vergAll,dvergAll,'.');
axis equal; axis([-20 20 -20 20])
xlabel('mean eye theta'); ylabel('delta mean eye theta')

figure
plot(nanmean(d_gz(d_hd(:,1)<5,:),1)); %%% select to get rid of crazy shit
hold on
plot(nanmean(d_hd(d_hd(:,1)<5,:),1));
legend('gaze drift','head drift');
xlabel('lag (frames)');
title('mean change over time - metric for gaze stabilization?')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
if savePDF
    pSname='T:\PreyCaptureAnalysis\Data\azApp\';
    filen=sprintf('%s','azAppAll_012820_medfilt','.pdf')
    pdfilename=fullfile(pSname,filen)
    dos(['ps2pdf ' psfilename ' ' pdfilename]);
    delete(psfilename);
else
    pFile='T:\PreyCaptureAnalysis\Data\azApp\';
end


% afilename=sprintf('AzApp_AllAnimals_121219_a','.mat');
% save(fullfile(pSname, afilename));