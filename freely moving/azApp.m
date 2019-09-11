close all
clear L H R
skip = 20;
nthresh  = 60;
for i = 1:length(appEpoch)
    vid = useData(i)
    
    azdeg = az{vid}*180/pi;
    app = appEpoch{i};
    
    thbins = -60:5:60;
    
    figure
    
    lth = Ltheta{vid} - nanmedian(Ltheta{vid});
    nl(i) = sum(~isnan(lth(app))); %%% # good eye approach points
    lthHist(:,1,i) = hist(lth(app),thbins)/nl(i);
    lthHist(:,2,i) = hist(lth(~app),thbins)/sum(~isnan(lth(~app)));
    subplot(2,3,1)
    plot(thbins,lthHist(:,:,i)); ylim([0 0.5]); legend('app','non-app'); title('left')
    title(sprintf('left %d pts app',sum(~isnan(lth(app)))))
    if nl(i)<nthresh      
        lthHist(:,:,i) = NaN;
    end
    
    maxlag = 30;
    
    rth = Rtheta{vid} - nanmedian(Rtheta{vid});
    nr(i) =sum(~isnan(lth(app))); %%% # good eye approach points
    rthHist(:,1,i) = hist(rth(app),thbins)/nr(i);
    rthHist(:,2,i) = hist(rth(~app),thbins)/sum(~isnan(rth(~app)));
    subplot(2,3,2)
    plot(thbins,rthHist(:,:,i)); ylim([0 0.5]); legend('app','non-app');
    title(sprintf('right %d pts app',sum(~isnan(rth(app)))))
    
    if nr(i)<nthresh      
        rthHist(:,:,i) = NaN;
    end
    
    lthAz_xc(:,i) = nanxcorr(azdeg(app),lth(app),maxlag,'zero');
    subplot(2,3,4)
    plot(-maxlag:maxlag,lthAz_xc(:,i)); ylim([-0.5 0.5])
    if nl(i)<nthresh
        lthAz_xc(:,i) = NaN;
    end
    
    rthAz_xc(:,i) = nanxcorr(azdeg(app),rth(app),maxlag,'zero');
    subplot(2,3,5)
    plot(-maxlag:maxlag,rthAz_xc(:,i));  ylim([-0.5 0.5])
    if nr(i)<nthresh
        rthAz_xc(:,i) = NaN;
    end
   
    %%% close all figs except on interval of skip
   
    
    hth = thetaHead{vid};
    dth = d_Theta{vid};
    n(i) = sum(~isnan(azdeg(app)) & ~isnan(dth(app))');
    dheadAz_xc(:,1,i) =nanxcorr(azdeg(app),dth(app),maxlag,'zero');
     dheadAz_xc(:,2,i) =nanxcorr(azdeg(~app),dth(~app),maxlag,'zero');
     subplot(2,3,6)
     plot(-maxlag:maxlag,dheadAz_xc(:,:,i)); axis([-maxlag maxlag -0.5 0.5]);
    title(sprintf('az vs dhead %d pts',n(i)))
     
    subplot(2,3,3)
    az_hist(:,1,i) = hist(-azdeg(app),thbins)/sum(~isnan(azdeg(app)));
    plot(thbins,az_hist(:,1,i),'b'); hold on
    
    n= sum(~isnan(azdeg(app)) & ~isnan(lth(app)'));
    azthL_hist(:,1,i) = hist(-azdeg(app)-lth(app)',thbins)/n;
    plot(thbins,azthL_hist(:,1,i),'r');
    n= sum(~isnan(azdeg(app)) & ~isnan(rth(app)'));
    azthR_hist(:,1,i) = hist(-azdeg(app)-rth(app)',thbins)/n;
     plot(thbins,azthR_hist(:,1,i),'g');
     
     
    %%% close all figs except on interval of skip
    if round(i/skip)~=i/skip, close(gcf),  end
    
    %%% alignment of eyes during approaches
    %%% vergence is cool! it gets very tight around 0 during approaches
%     drth = dRtheta{vid};
%     dlth = dLtheta{vid};
%     figure
%     subplot(2,2,1)
%     plot(dth(app),dlth(app),'.'); axis square; axis([-15 15 -15 15]); hold on; plot([-15 15],[15 -15])
%    subplot(2,2,2)
%         plot(dth(app),drth(app),'.'); axis square; axis([-15 15 -15 15]); hold on; plot([-15 15],[15 -15])
%     subplot(2,2,3)
%     plot(dth(~app),dlth(~app),'.'); axis square; axis([-15 15 -15 15]); hold on; plot([-15 15],[15 -15])
%    subplot(2,2,4)
%         plot(dth(~app),drth(~app),'.'); axis square; axis([-15 15 -15 15]); hold on; plot([-15 15],[15 -15])
% 
%         vergence = rth-lth;
%         n= sum(~isnan(vergence(app)));
%         vergeHist(:,1,i) = hist(vergence(app),thbins)/n;
%          vergeHist(:,2,i) = hist(vergence(~app),thbins)/sum(~isnan(vergence(~app)));
%         figure
%         plot(thbins,vergeHist(:,:,i))
%         
%     figure
%     subplot(2,1,1);
%     plot(dth); hold on
%    % plot(lth + hth');
%     subplot(2,1,2);
%     plot(hth); hold on
%     plot(rth+hth');
%     
    win = 2;
    thSum = conv(dth,ones(win,1));
    for rep = 1:2  %%% positive and negative saccades
        if rep==1
             saccEp = thSum>30;
        else
            saccEp = thSum<-30;
        end

    sacc = find(diff(saccEp)>0);
    sacc= sacc(sacc>31 & sacc<length(dth)-30);
    range = -30:30;
    clear Hsacc dthSacc Lsacc Rsacc  AZsacc
    baserange = -5:-1;
    for j = 1:length(sacc);
        Lsacc(:,j) = lth(sacc(j)+range) - nanmean(lth(sacc(j)+baserange));
         Rsacc(:,j) = rth(sacc(j)+range)- nanmean(rth(sacc(j)+baserange));
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
    end
    
   if sum(appSacc)>0 & sum(~appSacc)>0
    figure
    subplot(2,2,1)
    plot(range,Hsacc(:,~appSacc),'r');hold on;plot(range,nanmean(Hsacc(:,~appSacc),2),'r','LineWidth',4) 
    plot(range,Hsacc(:,appSacc),'b'); ylim([-90 90]); hold on; plot(range,nanmean(Hsacc(:,appSacc),2),'b','LineWidth',4)
    title('head')
    
%     subplot(2,2,2)
%     plot(range,dthSacc);ylim([-90 90]); hold on; plot(range,nanmean(dthSacc,2),'g','LineWidth',4)
%     title('dth')
    
subplot(2,2,2);
plot(range,nanmean(Hsacc(:,~appSacc),2),'k'); hold on
plot(range,nanmean(Lsacc(:,~appSacc),2),'b');
plot(range,nanmean(Rsacc(:,~appSacc),2),'g');
ylim([-90 90]); xlim([-20,20])
title(sprintf('rep %d',rep))

    subplot(2,2,3)
    plot(range,Lsacc); hold on; plot(range,nanmean(Lsacc,2),'g','linewidth',4); ylim([-45 45])
    title('left')
    
   subplot(2,2,4)
    plot(range,Rsacc); hold on; plot(range,nanmean(Rsacc,2),'g','linewidth',4); ylim([-45 45])
    title('right')
    
  
    %%% close all figs except on interval of skip
    if round(i/skip)~=i/skip, close(gcf),  end
    
    
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


figure
subplot(2,3,1)
plot(thbins,nanmean(lthHist,3)); ylim([0 0.5]); legend('app','non-app');

subplot(2,3,2)
plot(thbins,nanmean(rthHist,3)); ylim([0 0.5]); legend('app','non-app');

subplot(2,3,4)
plot(-maxlag:maxlag,nanmean(lthAz_xc,2));  ylim([-0.5 0.5])

subplot(2,3,5)
plot(-maxlag:maxlag,nanmean(rthAz_xc,2));  ylim([-0.5 0.5])

subplot(2,3,6)
plot(-maxlag:maxlag,(nanmean(dheadAz_xc,3))); ylim([-0.5 0.5])

subplot(2,3,3);
plot(thbins, nanmean(az_hist,3),'b')
hold on
plot(thbins, nanmean(azthL_hist,3),'r')
plot(thbins, nanmean(azthR_hist,3),'g')

for rep =1:2
figure
subplot(2,2,1);
plot(range,squeeze(H(:,rep,:))); axis([-20 20 -180 180])
subplot(2,2,3);
plot(range,squeeze(L(:,rep,:))); axis([-20 20 -45 45])
subplot(2,2,4);
plot(range,squeeze(R(:,rep,:))); axis([-20 20 -45 45])

subplot(2,2,2);
plot(range,squeeze(nanmean(H(:,rep,:,:),4)),'k'); hold on
plot(range,squeeze(nanmean(L(:,rep,:,:),4)),'b'); hold on
plot(range,squeeze(nanmean(R(:,rep,:,:),4)),'r'); hold on
axis([-20 20 -45 45])

end


