close all
clear L H R
skip = 5; %%% only shows figures at this interval
nthresh  = 60;  %%% threshold for number of approach points to be included in averages

vgDiff = [];
hdDiff = [];
gzDiff = [];
vergAll = [];
for i = 1:length(appEpoch)
    vid = useData(i)
    
    azdeg = az{vid}*180/pi;
    app = appEpoch{i};
    
    thbins = -60:5:60;
    
    %figure
    
    lth = Ltheta{vid} - nanmedian(Ltheta{vid});
    nl(i) = sum(~isnan(lth(app))); %%% # good eye approach points
    lthHist(:,1,i) = hist(lth(app),thbins)/nl(i);
    lthHist(:,2,i) = hist(lth(~app),thbins)/sum(~isnan(lth(~app)));
    %     subplot(2,3,1)
    %     plot(thbins,lthHist(:,:,i)); ylim([0 0.5]); legend('app','non-app'); title('left')
    %     title(sprintf('left %d pts app',sum(~isnan(lth(app)))))
    if nl(i)<nthresh
        lthHist(:,:,i) = NaN;
    end
    
    maxlag = 30;
    
    rth = Rtheta{vid} - nanmedian(Rtheta{vid});
    nr(i) =sum(~isnan(lth(app))); %%% # good eye approach points
    rthHist(:,1,i) = hist(rth(app),thbins)/nr(i);
    rthHist(:,2,i) = hist(rth(~app),thbins)/sum(~isnan(rth(~app)));
    %     subplot(2,3,2)
    %     plot(thbins,rthHist(:,:,i)); ylim([0 0.5]); legend('app','non-app');
    %     title(sprintf('right %d pts app',sum(~isnan(rth(app)))))
    %
    if nr(i)<nthresh
        rthHist(:,:,i) = NaN;
    end
    
    if  nl(i)>nthresh
    lthAz_xc(:,i) = nanxcorr(azdeg(app),lth(app),maxlag,'zero');
    %     subplot(2,3,4)
    %     plot(-maxlag:maxlag,lthAz_xc(:,i)); ylim([-0.5 0.5])
    else
        lthAz_xc(:,i) = NaN;
    end
    
    if nr(i)>nthresh
    rthAz_xc(:,i) = nanxcorr(azdeg(app),rth(app),maxlag,'zero');
    %     subplot(2,3,5)
    %     plot(-maxlag:maxlag,rthAz_xc(:,i));  ylim([-0.5 0.5])
    else
        rthAz_xc(:,i) = NaN;
    end
    
    %%% close all figs except on interval of skip
    
    
    hth = thetaHead{vid};
    dth = d_Theta{vid};
    n(i) = sum(~isnan(azdeg(app)) & ~isnan(dth(app))');
    if n(i)>0
    dheadAz_xc(:,1,i) =nanxcorr(azdeg(app),dth(app),maxlag,'zero');
    dheadAz_xc(:,2,i) =nanxcorr(azdeg(~app),dth(~app),maxlag,'zero');
    else
        dheadAz_xc(:,1,i)=NaN;
        dheadAz_xc(:,2,i)=NaN;
    end
    %      subplot(2,3,6)
%          plot(-maxlag:maxlag,dheadAz_xc(:,:,i)); axis([-maxlag maxlag -0.5 0.5]);
%         title(sprintf('az vs dhead %d pts',n(i)))
    
    
    az_hist(:,1,i) = hist(-azdeg(app),thbins)/sum(~isnan(azdeg(app)));
    
    
    n= sum(~isnan(azdeg(app)) & ~isnan(lth(app)'));
    azthL_hist(:,1,i) = hist(-azdeg(app)-lth(app)',thbins)/n;
    
    n= sum(~isnan(azdeg(app)) & ~isnan(rth(app)'));
    azthR_hist(:,1,i) = hist(-azdeg(app)-rth(app)',thbins)/n;
    %
    %      subplot(2,3,3)
    %     plot(thbins,az_hist(:,1,i),'b'); hold on
    %      plot(thbins,azthL_hist(:,1,i),'r');
    %      plot(thbins,azthR_hist(:,1,i),'g');
    
    
    %%% close all figs except on interval of skip
    % if round(i/skip)~=i/skip, close(gcf),  end
    
    %%% alignment of eyes during approaches
    %%% vergence is cool! it gets very tight around 0 during approaches
    drth = dRtheta{vid};
    dlth = dLtheta{vid};
    
    vergence = rth-lth;
    n= sum(~isnan(vergence(app)));
    vergeHist(:,1,i) = hist(vergence(app),thbins)/n;
    vergeHist(:,2,i) = hist(vergence(~app),thbins)/sum(~isnan(vergence(~app)));
    %         subplot(2,3,6)
    %         plot(thbins,vergeHist(:,:,i))
    
    az_hist(:,1,i) = hist(-azdeg(app),thbins)/sum(~isnan(azdeg(app)));
    gz = 0.5*(rth+lth);
    
    n= sum(~isnan(azdeg(app)) & ~isnan(gz(app)'));
    azthRL_hist(:,1,i) = hist(-azdeg(app)-gz(app)',thbins)/n;
    %
    %        figure
    %        subplot(2,2,1)
    %
    %     plot(-azdeg(app),gz(app),'.')
    %     axis square ; axis([-45 45 -45 45])
    %
    %     subplot(2,2,2)
    %             plot(thbins,az_hist(:,1,i),'k'); hold on
    %          plot(thbins,azthRL_hist(:,1,i),'g');
    
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
    
    hthApp = hth(appRange)-nanmedian(hth(mainApp));
 hthApp = mod(hthApp + 180,360)-180;
    
    gz = hthApp +0.5*( rth(appRange) +lth(appRange))';
    vg =  0.5*( rth(appRange) +lth(appRange))';
    vgDiff = [vgDiff diff(vg)];
    hdDiff = [hdDiff diff(hthApp)];
    gzDiff = [gzDiff diff(gz)];
    vergAll =[vergAll vg(1:end-1)];
    
    
    %%% calculate change in position at different lags, as measure of stability
    for lag = 1:20;
        d_gz(i,lag) = nanmean(abs(gz(1:end-lag) - gz((lag+1):end)));
        d_hd(i,lag) = nanmean(abs(hthApp(1:end-lag) - hthApp((lag+1):end)));
    end
    
    
  
    if round(i/skip)==i/skip
       
%             win = 15;
%     xc = slidingXC(dth(1:end-1)', drth+dlth,10,-win:win);
%     figure
%     subplot(7,5,16:35)
%    imagesc(xc',[-1 1]); colormap jet; title('xcorr dhead dEye theta'); xlabel('frames')
%   subplot(7,5,1:5);
%     plot(mouseSp{vid}); xlim([1 size(xc,1)]); ylim([0 50]); 
%     title('speed')   
%     subplot(7,5,6:10);
%     plot(azdeg(win:(end-win)));xlim([1 size(xc,1)]); ylim([-180 180])
%     title('azimuth')  
%     subplot(7,5,11:15);
%     plot(app(win:(end-win)),'g','Linewidth',2); ylim([0.5 1.2]); xlim([1 size(xc,1)])
%     title('approach')
%     
%     
%     
%     win = 15;
%      xc = slidingXC(dth', (azdeg),10,-win:win);
%     figure
%       subplot(7,5,16:35)
%     imagesc(xc',[-1 1]); colormap jet; title('xcorr dHead vs Azimuth'); xlabel('frames')
%     
%     subplot(7,5,1:5);
%     plot(mouseSp{vid}); xlim([1 size(xc,1)]); ylim([0 50]); 
%     title('speed')   
%     subplot(7,5,6:10);
%     plot(azdeg(win:(end-win)));xlim([1 size(xc,1)]); ylim([-180 180])
%     title('azimuth')  
%     subplot(7,5,11:15);
%     plot(app(win:(end-win)),'g','Linewidth',2); ylim([0.5 1.2]); xlim([1 size(xc,1)])
%     title('approach')
        
        
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

        %%% pursuit eye/head/gaze positions
        
        figure
        subplot(4,1,1);
        plot(hthnonan,'k'); hold on; plot(rth,'r'); plot(lth,'b'); legend('head th','right th','left th');
        plot(find(app),ones(sum(app),1)*90,'g.'); ylim([-180 180])
           title(sprintf('vid %d',vid));
           
        subplot(4,1,2);
        plot(0.5*( rth(appRange) +lth(appRange))','k','LineWidth',2); hold on; plot( rth(appRange)','r'); hold on; plot(lth(appRange)','b');
        ylim([-30 30]);  xlim([0 max(length(appRange),1)]); plot(appOffset,0,'g*'); plot(endOffset,0,'r*');
        legend('mean','right','left');
        
        subplot(4,1,3);
        hold on; plot(hthApp +0.5*( rth(appRange) +lth(appRange))','k','LineWidth',2);   plot(hthApp,'g');       
        plot(appOffset,0,'g*'); plot(endOffset,0,'r*');
        ylim([-60 60]); xlim([0 max(length(appRange),1)]); xlabel('frames'); ylabel('deg');
        legend('gaze','head')
        
        subplot(4,1,4); hold on 
        plot(azdeg(appRange)); plot(mouseSp{vid}(appRange)*2); ylim([-90 90]); 
        plot(appOffset,0,'g*'); xlim([0 max(length(appRange),1)]); plot(endOffset,0,'r*');
        legend('azimuth','speed');
        
        set(gcf,'Position',[440 100 560 640])
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
                figure
                subplot(2,2,1)
                plot(range,Hsacc(:,~appSacc),'r');hold on;plot(range,nanmean(Hsacc(:,~appSacc),2),'r','LineWidth',4)
                plot(range,Hsacc(:,appSacc),'b'); ylim([-90 90]); hold on; plot(range,nanmean(Hsacc(:,appSacc),2),'b','LineWidth',4)
                title('head')
                
                %     subplot(2,2,2)
                %     plot(range,dthSacc);ylim([-90 90]); hold on; plot(range,nanmean(dthSacc,2),'g','LineWidth',4)
                %     title('dth')
                
                subplot(2,2,2);
                plot(range,nanmean(Hsacc(:,appSacc),2),'k'); hold on
                plot(range,nanmean(Lsacc(:,appSacc),2),'b');
                plot(range,nanmean(Rsacc(:,appSacc),2),'g');
                ylim([-90 90]); xlim([-20,20])
                title(sprintf('rep %d',rep))
                
                subplot(2,2,3)
                plot(range,Lsacc(:,~appSacc),'r'); hold on; plot(range,Lsacc(:,appSacc),'b');
                plot(range,nanmean(Lsacc(:,~appSacc),2),'r','linewidth',4); plot(range,nanmean(Lsacc(:,appSacc),2),'b','linewidth',4);
                ylim([-45 45])
                title('left')
                
                subplot(2,2,4)
                
                plot(range,Rsacc(:,~appSacc),'r'); hold on; plot(range,Rsacc(:,appSacc),'b');
                plot(range,nanmean(Rsacc(:,~appSacc),2),'r','linewidth',4); plot(range,nanmean(Rsacc(:,appSacc),2),'b','linewidth',4);
                ylim([-45 45])
                title('right')
                
                
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
    %     figure
    % for j = 1:(min(24,length(sacc)))
    %     subplot(4,6,j);
    %     plot(Hsacc(:,j)); hold on;
    %     plot(Lsacc(:,j)); plot(Rsacc(:,j));
    %     ylim([-90 90])
    % end
    
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
    figure
    subplot(2,2,1);
    plot(range,squeeze(H(:,rep,1,:))); axis([-20 20 -90 90])
    subplot(2,2,3);
    plot(range,squeeze(L(:,rep,1,:))); axis([-20 20 -45 45])
    subplot(2,2,4);
    plot(range,squeeze(R(:,rep,1,:))); axis([-20 20 -45 45])
    
    subplot(2,2,2);
    %plot(range,squeeze(nanmean(H(:,rep,1,:),4)),'k'); hold on
    plot(range,squeeze(nanmean(L(:,rep,1,:),4)),'b'); hold on
    plot(range,squeeze(nanmean(R(:,rep,1,:),4)),'r'); hold on
    % plot(range,squeeze(nanmean(H(:,rep,2,:),4)),'k:'); hold on
    % plot(range,squeeze(nanmean(L(:,rep,2,:),4)),'b:'); hold on
    % plot(range,squeeze(nanmean(R(:,rep,2,:),4)),'r:'); hold on
    axis([-30 30 -20 20])
    title(sprintf('rep %d',rep))
    
end

hbins = -20:0.5:20;
figure
plot(hbins,hist(hdDiff +vgDiff,hbins))
hold on
plot(hbins,hist(hdDiff,hbins))
legend('delta gaze','delta head')
xlabel('deg')


figure
plot(hdDiff,gzDiff,'.');
axis equal; axis([-25 25 -25 25]);
hold on; plot([-25 25],[-25 25],'r')
xlabel('delta head'); ylabel('delta gaze');

figure
plot(hdDiff,vgDiff,'.');
axis equal; hold on
axis([-20 20 -20 20]);
plot([-20 20],[20 -20],'r')
xlabel('delta head'); ylabel('delta theta mean eye');


figure
plot(vergAll,gzDiff,'.');
axis equal; axis([-20 20 -20 20])
xlabel('mean eye theta'); ylabel('delta gaze')

figure
plot(vergAll,vgDiff,'.');
axis equal; axis([-20 20 -20 20])
xlabel('mean eye theta'); ylabel('delta mean eye theta')

figure
plot(nanmean(d_gz(d_hd(:,1)<5,:),1)); %%% select to get rid of crazy shit
hold on
plot(nanmean(d_hd(d_hd(:,1)<5,:),1));
legend('gaze drift','head drift');
xlabel('lag (frames)');
title('mean change over time - metric for gaze stabilization?')