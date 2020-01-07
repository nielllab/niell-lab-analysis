
%% FIGURE 1: a freely moving eye and head tracking system
% Figure 1G: gyro & DLC traces

% Figure 1H: measure of similarity of gyro & DLC (scatter or hist of
% difference)

figure
clear h
bins=-20:2:20
h=hist(gyro3All-dlcDhth,bins)
plot(bins,h/length(gyro3All)); hold on


% Figure 1I (where is the right place to put this?): head angle is directed
% towards cricket - corr of az and change in head yaw

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
for vid=1:length(useData)
    clear hT azC use useN corrR corrRA corrL corrLA
    nframe = min(length(az{useData(vid)}),length(d_Theta{useData(vid)}));
    nframe = min(nframe,length(appEpoch{vid}));
    useN=appEpoch{vid}(1:nframe)==0; use=appEpoch{vid}(1:nframe)==1;

    azC=-(az{useData(vid)}); hT=d_Theta{useData(vid)}; hT=hT-nanmedian(hT);
    
    if sum(useN)>3 & sum(~isnan(hT(useN)))>20 & sum(~isnan(azC(useN)))>20
        [corrR lagsR]= nanxcorr(azC(useN),hT(useN),30,'coeff');
        hold on;
    else
            corrR=NaN;
    end
    
%     subplot(4,5,vid)
%     plot(rad2deg(azC(use))); hold on; plot(hT(use));
   if sum(use)>4 & sum(~isnan(hT(use)))>20
        [corrRA lagsRA]= nanxcorr(azC(use),hT(use),30,'coeff');
   else
              corrRA=NaN;
   end
        corrRAll(vid,:)=corrR; 
        corrRAAll(vid,:)=corrRA; 
end

figure%('units','normalized','outerposition',[0 0 1 1])
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); 
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); xlim([1 61]); axis square; ylim([-.2 .5]);

% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('az, eye theta both eyes');

title('az and dHead Theta')
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

%% FIGURE 2: Coordination of eyes during free movement

% Figure 2A: scatter plots of R & L eyes yaw vs pitch, approach and non-approach
%%
figure
for vid=1:length(useData)
    nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
    pR=Rphi{useData(vid)}(1:nframe); pL=Lphi{useData(vid)}(1:nframe);
    pR=nanmean(pR)-pR; pL=nanmean(pL)-pL;
    tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
    tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
    for c=0:1
        clear use
        use= appEpoch{vid}==c;
        subplot(1,2,1)
        plot(tR(use(1:20:length(use))),pR(use(1:20:length(use))),'.'); hold on; xlim([-60 60]);ylim([-60 60]); axis square
        title('R eye')
        subplot(1,2,2)
        plot(tL(use(1:20:length(use))),pL(use(1:20:length(use))),'.'); hold on; xlim([-60 60]);ylim([-60 60]); axis square
        title('L eye')
    end
    xlabel('yaw'); ylabel('pitch')
end
%%
% Figure 2B: histograms of yaw and pitch from panel above
clear tNcounts tAcounts allThN allThA
 figure
for vid=1:length(useData)
    nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
    pR=Rphi{useData(vid)}(1:nframe); pL=Lphi{useData(vid)}(1:nframe);
    pR=nanmean(pR)-pR; pL=nanmean(pL)-pL;
    tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
    tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
    useN= appEpoch{vid}==0; use= appEpoch{vid}==1;
    subplot(1,2,1)
    plot(tR(useN(1:20:end)),pR(useN(1:20:end)),'b.'); hold on; axis square
    plot(tR(use(1:20:end)),pR(use(1:20:end)),'g.'); hold on
     subplot(1,2,2);
      plot(tL(useN(1:20:end)),pL(useN(1:20:end)),'b.'); hold on; axis square
      plot(tL(use(1:20:end)),pL(use(1:20:end)),'g.'); hold on; axis square

    clear r l
    if sum(useN)>4
    nBins=-90:20:90; nBinsP=-50:5:50
    r= (hist(tR(useN),nBins))/sum(useN);l= (hist(tL(useN),nBins))/sum(useN);
     rP= (hist(pR(useN),nBinsP))/sum(useN);lP= (hist(pL(useN),nBinsP))/sum(useN);
    else
        r=NaN;l=NaN;
        rP=NaN;lP=NaN;
    end
   tNcounts(vid,:,1)=r;
   tNcounts(vid,:,2)=l;
   tNcountsP(vid,:,1)=rP;
   tNcountsP(vid,:,2)=lP;
    
end

for vid=1:length(useData)
    nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
    tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
    tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
   pR=Rphi{useData(vid)}(1:nframe); pL=Lphi{useData(vid)}(1:nframe);
    pR=nanmean(pR)-pR; pL=nanmean(pL)-pL;
    use = (appEpoch{vid});
%     plot(tR(use(1:20:end)),tL(use(1:20:end)),'go'); hold on; axis square; %ylim([-30 30]); xlim([-30 30]);
    clear r l
    if sum(use)>4
        nBins=-90:20:90; nBinsP=-50:5:50
        r= (hist(tR(use),nBins))/sum(use);l= (hist(tL(use),nBins))/sum(use);
        rP= (hist(pR(use),nBinsP))/sum(use);lP= (hist(pL(use),nBinsP))/sum(use);

    else
        r=NaN;l=NaN;
        rP=NaN;lP=NaN;
    end
    tAcounts(vid,:,1)=r;
    tAcounts(vid,:,2)=l;
    tAcountsP(vid,:,1)=rP;
    tAcountsP(vid,:,2)=lP;
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


allThNP(:,1)=nansum(tNcountsP(:,:,1),1)./(length(useData))
allThNP(:,2)=nansum(tNcountsP(:,:,2),1)./(length(useData))

allThAP(:,1)=nansum(tAcountsP(:,:,1),1)./(length(useData))
allThAP(:,2)=nansum(tAcountsP(:,:,2),1)./(length(useData))

figure;subplot(1,2,1); plot(nBinsP,allThNP(:,1),'b'); hold on; plot(nBinsP,allThAP(:,1),'g'); title('Right Eye Phi'); axis square
ylim([0 .32]);
subplot(1,2,2);plot(nBinsP,allThNP(:,2),'b'); hold on; plot(nBinsP,allThAP(:,2),'g');title('Left Eye Phi');axis square
ylim([0 .32]);
legend('non approach','approach')

%%
% Figure 2C: overlaid trace of two eyes converging and diverging

% Figure 2D: scatter plot of R vs L eye yaw to show conv/divergence 
clear numPts congMv convMv divMv Rmv Lmv
figure
    for vid=1:length(useData)
        clear nframe
        nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
        nframe=min(nframe, length(appEpoch{vid}));
        tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
        tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
        useN = appEpoch{vid}(1:nframe)==0; use = appEpoch{vid}(1:nframe)==1;
        plot(tR(useN(1:20:end)),tL(useN(1:20:end)),'bo'); hold on; axis square; xlim([-50 50]);ylim([-50 50]);
        plot(tR(use(1:20:end)),tL(use(1:20:end)),'go');
%         numPts(vid,1)=length(tR(1:20:end));
%         numPts(vid,2)=length(tR);

    Rmv=nansum(tR(useN)>=0& tL(useN)>=0);
    Lmv=nansum(tR(useN)<=0& tL(useN)<=0);
    congMv(vid,:)=(Rmv+Lmv)./(nansum(useN & ~isnan(tR)'& ~isnan(tL)'));
    convMv(vid,:)=(nansum(tR(useN)<=0& tL(useN)>=0))./sum(useN&~isnan(tR)'&~isnan(tL)');
    divMv(vid,:)=(nansum(tR(useN)>=0 & tL(useN)<=0))./sum(useN&~isnan(tR)'&~isnan(tL)');
    end

% Figure 2E: quantification of 2D, bar plots


% Figure 2F: hist of difference in eye yaw (aka vergence)

clear vergCounts vgDiff nBins vgHist

for c=0:1
    for vid=1:length(useData)
        nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
        nframe=min(nframe, length(appEpoch{vid}))
        tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
        tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
        vg=tR-tL;
        useN= appEpoch{vid}(1:nframe)==c;
        nBins=-90:5:90
        if sum(useN)>4
            vgHist = (hist(vg(useN),nBins))/sum(useN&~isnan(vg)');
        else
            vgHist=NaN;
        end
        vgDiff(vid,:,c+1)=vgHist;
    end
end

vergCounts(:,1)=nansum(vgDiff(:,:,1),1)./(length(useData))
vergCounts(:,2)=nansum(vgDiff(:,:,2),1)./(length(useData))

figure;plot(nBins,vergCounts(:,1),'b');
hold on; 
plot(nBins,vergCounts(:,2),'g'); title('diff in vergence'); axis square




% Figure 2G: correlation of eye thetas & eye phis 
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
corrRAAll=[]
for vid=1:length(useData)
    nframe = min(length(dRtheta{useData(vid)}),length(dLtheta{useData(vid)}));
    dtR=dRtheta{useData(vid)}(1:nframe); dtL=dLtheta{useData(vid)}(1:nframe);
    nonapp=appEpoch{vid}==0;
    use = (nonapp==1)';%& ~isnan(dtR(1:nframe));
    if sum(use)>3
        [corrR lagsR]= nanxcorr(dtR(use),dtL(use),30,'coeff');
        uselagsR=(lagsR>=-30& lagsR<=30);

        clear nframe
        nframe = min(length(dRphi{useData(vid)}),length(dLphi{useData(vid)}));
        dpR=dRphi{useData(vid)}(1:nframe); dpL=dLphi{useData(vid)}(1:nframe);
        [corrL lagsL]= nanxcorr(dpR(use),dpL(use),30,'coeff');
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    use=appEpoch{vid}==1;
    if sum(use)>3 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dtR(use),dtL(use),30,'coeff');
        uselagsRA=(lagsRA>=-30& lagsRA<=30);

        nframe = min(length(dRphi{useData(vid)}),length(dLphi{useData(vid)}));
        dpR=dRphi{useData(vid)}(1:nframe); dpL=dLphi{useData(vid)}(1:nframe);
        [corrLA lagsLA]= nanxcorr(dpR(use),dpL(use),30,'coeff');
        uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    else
    end

    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
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

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.7 .7]);
xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dTheta non-app','dPhi non-app','dTh Approach','dPhi app'}); title('mean between eye corr, dtheta & dphi');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%% FIGURE 3: during approaches, eyes are more centered (smaller vergence), due to less head movements in pitch and roll

% Figure 3A: overlaid traces of pitch & eye vergence

% Figure 3B: scatter & corr of pitch (tilt) and vergence
for c=0:1
    clear use
use = find(appAll==c);

figure(3)
plot(tiltAll(use),vergDlc(use),'.'); axis equal; hold on; lsline
xlabel('acc tilt'); ylabel('eye vergence');
ylim([-90 90]); xlim([-90 90]);
R=corrcoef(tiltAll(use), vergDlc(use),'Rows','pairwise')
text(-80,(80-c*5), ['corrcoef = ' num2str(R(1,2),'%.2f')],'FontSize',10)
title('acc ch 2');
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end


figure(4)
[corr lags]=(nanxcorr(tiltAll(use),vergDlc(use),30,'coeff'));
plot(lags,corr); hold on
axis square;ylim([-1 1]);
title('acc tilt, vergence');
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end

end

% Figure 3C: example trace, less change in pitch during approach

% Figure 3D: hist or corr, less change in pitch during approach

% Figure 3E: overlaid traces of roll & eye phi

% Figure 3F: scatter & corr of roll & eye phi
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
[corr lags]=(nanxcorr(rollAll(use),dlcPhi(use),30,'coeff'));
plot(lags,corr);axis square; hold on
title('acc roll, phi');
ylim([-1 1]);
if c==1
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
end

end

% Figure 3G: example trace, less change in roll during approach

% Figure 3H: hist or corr, less change in roll during approach

%% FIGURE 4: Coordination of eyes and head

% Figure 4A: hist of head yaw app and non approach (not different)

% Figure 4B: hist of mean eye yaw, app and non app 

% Figure 4C: ex trace when dTh is 0, eyes are also 0, then hist

% Figure 4D: scatter of change in head yaw and change in eye yaw, shows
% mostly congruent but not all

% Figure 4E: corr of change in head and eye yaw, congruence
% (non-compensation) at short time scales

%% FIGURE 5: Non-compensatory eye movements

% Figure 5A: same as 4D (change in eye and head yaw), but separated by eye mvmt types from gmm
  % maybe bar plot showing prop of time non-app and app that each occurs

% Figure 5B: example of target selection saccades 

% Figure 5C: quantification of target selection sacc

% Figure 5D: example of resetting saccades

% Figure 5E: quantification of resetting sacc











