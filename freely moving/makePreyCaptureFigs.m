
%% FIGURE 1: a freely moving eye and head tracking system
% Figure 1G: gyro & DLC traces

% Figure 1H: measure of similarity of gyro & DLC (scatter or hist of
% difference)

% Figure 1I (where is the right place to put this?): head angle is directed
% towards cricket - corr of az and change in head yaw

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
%% FIGURE 2: Coordination of eyes during free movement

% Figure 2A: scatter plots of R & L eyes yaw vs pitch, approach and non-approach

% Figure 2B: histograms of yaw and pitch from panel above
clear tNcounts tAcounts allThN allThA
% figure
for vid=1:length(useData)
    nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
    pR=Rphi{useData(vid)}(1:nframe); pL=Lphi{useData(vid)}(1:nframe);
    pR=nanmean(pR)-pR; pL=nanmean(pL)-pL;
    tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe);
    tR=nanmean(tR)-tR; tL=nanmean(tL)-tL;
    useN= appEpoch{vid}==0;
    % subplot(1,2,1)
%     plot(tR(useN(1:70:end)),tL(useN(1:70:end)),'bo'); hold on; axis square
    %plot(tR(use(1:70:end)),pR(use(1:70:end)),'go'); hold on
    %  subplot(1,2,2);
    %   plot(tL(useN(1:70:end)),pL(useN(1:70:end)),'bo'); hold on; axis square
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
%     plot(tR(use(1:70:end)),tL(use(1:70:end)),'go'); hold on; axis square; %ylim([-30 30]); xlim([-30 30]);
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

% Figure 2C: overlaid trace of two eyes converging and diverging

% Figure 2D: scatter plot of R vs L eye yaw to show conv/divergence 

% Figure 2E: quantification of 2D, bar plots

% Figure 2F: hist of difference in eye yaw (aka vergence)

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

% Figure 3C: example trace, less change in pitch during approach

% Figure 3D: hist or corr, less change in pitch during approach

% Figure 3E: overlaid traces of roll & eye phi

% Figure 3F: scatter & corr of roll & eye phi

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











