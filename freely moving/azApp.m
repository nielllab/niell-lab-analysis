close all
clear all

load('Analyzed_AllAnimals_120219_a.mat')
clear L H R
skip = 5; %%% only shows figures at this interval
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

appAll=[];

maxlag = 30;
thbins = -60:5:60;

for i = 1:length(appEpoch)
    
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
    nr(i) =sum(~isnan(lth(app))); %%% # good eye approach points
    rthHist(:,1,i) = hist(rth(app),thbins)/nr(i);
    rthHist(:,2,i) = hist(rth(~app),thbins)/sum(~isnan(rth(~app)));
    if nr(i)<nthresh
        rthHist(:,:,i) = NaN;
    end
    
    %%% get head positions
    hth = thetaHead{vid};
    dth = d_Theta{vid};
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
        subplot(3,1,1);
        plot(hthnonan,'k'); hold on; plot(rth,'r'); plot(lth,'b'); legend('head th','right th','left th');
        plot(find(app),ones(sum(app),1)*90,'g.'); ylim([-180 180])
        title(sprintf('vid %d',vid));
        
        subplot(3,1,2);
        plot(0.5*( rth(appRange) +lth(appRange))','k','LineWidth',2); hold on; plot( rth(appRange)','r'); hold on; plot(lth(appRange)','b');
        ylim([-30 30]);  xlim([0 max(length(appRange),1)]);
        plot([appOffset appOffset],[-60 60],'g'); plot([endOffset endOffset],[-60 60],'r');
        legend('mean','right','left');
        
        subplot(3,1,3);
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

%%% clean up a couple values
diffThAll(diffThAll<-180) = diffThAll(diffThAll<-180)+360;
diffThAll(diffThAll>180) = diffThAll(diffThAll>180)-360;
dgazeAll = dEyeAll + dthAll;
appAll = appAll==1;

%%% large head movement seems to be DLC errors according to accelerometers
dthAll(abs(dthAll)>25)=NaN;

%%% plot eye, gaze, head versus each other
%%% hopefully can be used for clustering
figure
plot(dthAll(appAll),dEyeAll(appAll),'.');
axis square; axis([-25 25 -25 25]); hold on; plot([-10 10], [10 -10],'r'); plot([-10 10], [-10 10],'r')

figure
plot(nanxcorr(dthAll(appAll),dEyeAll(appAll),30,'coeff'))

figure
plot(dthAll(appAll),dgazeAll(appAll),'.');
axis equal; axis equal; hold on; plot([-10 10], [-10 10],'r')


figure
plot(dthAll(appAll)-dEyeAll(appAll),dthAll(appAll)+dEyeAll(appAll),'.')


figure
plot(mnEyeAll(appAll),dgazeAll(appAll),'.');

figure
plot(mnEyeAll(appAll),dEyeAll(appAll),'.');

figure
plot(mnEyeAll(appAll),dthAll(appAll),'.');


figure
plot(max(lthAll(appAll), rthAll(appAll)),dgazeAll(appAll),'.');



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

figure;plot(randsample(dgz,20000),randsample(dvg,20000),'.'); hold on;
axis equal; axis([-40 40 -40 40])
plot([-90 90],[0 0]);
plot([-0 0],[-50 50])
plot([-100 100],[-50 50])

plot(dgz(find(r)),dvg(find(r)),'.c')
plot(dgz(find(s)),dvg(find(s)),'.g')
plot(dgz(find(ht)),dvg(find(ht)),'.r')
axis equal; axis([-40 40 -40 40])
xlabel('d gaze');ylabel ('d eye theta')
title('non-approach')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

figure;plot(randsample(dgz(find(appAll)),20000),randsample(dvg(find(appAll)),20000),'.'); hold on;
plot(dgz(find(appAll&r)),dvg(find(appAll&r)),'.c')
plot(dgz(find(appAll&s)),dvg(find(appAll&s)),'.g')
plot(dgz(find(appAll&ht)),dvg(find(appAll&ht)),'.r')
axis equal; axis([-40 40 -40 40])
plot([-90 90],[0 0]);
plot([-0 0],[-50 50])
plot([-100 100],[-50 50])
xlabel('d gaze');ylabel ('d eye theta');
title('approaches')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%

pts = randsample(1:length(dht),20000);
head = dht(pts); eyes = dEyeAll(pts);

he=[dht+dEyeAll;dht-dEyeAll];


% gm = fitgmdist(he',6);
gm = fitgmdist(he',5);
he=[head+eyes;head - eyes];
idx = cluster(gm,he')
X=he;

figure;
gscatter(X(1,:),X(2,:),idx);
legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Location','best');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end


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
figure;
subplot(1,2,1)
plot(randsample(dht,15000),randsample(dEyeAll,15000),'.'); hold on;
plot(dht(find(r)),dvg(find(r)),'.c')
plot(dht(find(s)),dvg(find(s)),'.g')
plot(dht(find(ht)),dvg(find(ht)),'.r')
axis equal; axis([-40 40 -40 40])
% plot([-90 90],[0 0]);
% plot([-0 0],[-50 50])
% plot([-100 100],[-50 50])
xlabel('d head th');ylabel ('d eye theta')
title('non-approach')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

% for approaches
% figure;
subplot(1,2,2);
plot(randsample(dht(find(appAll)),15000),randsample(dEyeAll(find(appAll)),15000),'.'); hold on;
% plot(dht(find(appT&r)),dvg(find(appT&r)),'.c')
% plot(dht(find(appT&s)),dvg(find(appT&s)),'.g')
% plot(dht(find(appT&ht)),dvg(find(appT&ht)),'.r')
axis equal; axis([-40 40 -40 40])
% plot([-90 90],[0 0]);
% plot([-0 0],[-50 50])
% plot([-100 100],[-50 50])
xlabel('d head th');ylabel ('d eye theta');
title('approaches')
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
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
    filen=sprintf('%s','azAppAll_121219_a','.pdf')
    pdfilename=fullfile(pSname,filen)
    dos(['ps2pdf ' psfilename ' ' pdfilename]);
    delete(psfilename);
else
    pFile='T:\PreyCaptureAnalysis\Data\azApp\';
end


afilename=sprintf('AzApp_AllAnimals_121219_a','.mat');
save(fullfile(pSname, afilename));