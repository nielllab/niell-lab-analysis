clear all; close all;
% set(groot,'defaultFigureVisible','off') %disable figure plotting

savePDF=1;
if savePDF
    psfilename = 'C:\analysisPS_all.ps';
    if exist(psfilename,'file')==2;delete(psfilename);end
end


pname = 'T:\PreyCaptureAnalysis\Data\';
files={'J462aAnalyzed_All_090319_new.mat','J462bAnalyzed_All_090319_new.mat','J462cAnalyzed_All_090319_new.mat',...
    'J463bAnalyzed_All_090319_new.mat','J463cAnalyzed_All_090319_new.mat','J465dAnalyzed_All_090319_new.mat',...
    'J470bAnalyzed_All_090319_new.mat','J475cAnalyzed_All_090319_new.mat'}%,'J470bAnalyzed_All_090319_new.mat'}; %'J462bAnalyzed_All_090319_1.mat', 'J463dAnalyzed_All_090319_1.mat',
fname=[];cricketPos=[];mousePos=[];cricketSp=[];mouseSp=[];thetaHead=[];
d_Theta=[];az=[];dist=[];Rtheta=[];Ltheta=[];Rphi=[];Lphi=[];Rcent=[];Lcent=[];n_vid=[]; mousePos=[];theta_good=[];approach=[];useData=[];
fitParamsR=[];fitParamsL=[]; dRtheta=[];dLtheta=[]; dRphi=[];dLphi=[]; dxR=[];dyR=[];Rcent=[];Lcent=[];radR=[];radL=[];Rgood=[];ccR=[];slopeR=[];
scaleR=[]; Lgood=[];ccL=[];slopeL=[];scaleL=[]; ani=[]; sess=[]; date=[];clipnum=[]; cricketPosHead=[]; cricketSpHead=[];
n=0;
for i = 1:length(files)
    
%     clear 'mouse_xy','cricketV','theta','slip','dTheta','azT', 'mouseV', 'appEpoch','useData','EllipseParamsR','EllipseParamsL',...
%         'dthetaR','dthetaL','dphiR','dphiL','dXRcent';
    fname= fullfile(pname,files{i});
    load(fname,'cricket_xy','mouse_xy','cricketV','cricketH','theta', 'dTheta','azT','mouseV','goodTheta','range','theta','dTheta','appEpoch','useData','EllipseParamsR',...
        'EllipseParamsL','thetaR','thetaL','phiR','phiL','dthetaR','dthetaL','dphiR','dphiL','dXRcent','XRcent','YRcent','YLcent','XLcent','RRad','LRad',...
        'goodR','Rcc','Rslope','Rscale','Rngood','goodL','Lcc','Lslope','Lscale','Lngood','animal','sessionN','clip','expdate');
    nc = length(cricket_xy(:,1))'; vidrange= n+1:n+nc;
    nvid{i} = length(cricket_xy(:,1))';
    %  nvidTot = nvid(i+1)
    %   delayFull{end+1,1} = cell2mat(slip);
    
    for j = 1:length(vidrange)
        n_vid(end+1) = j;
        ani{end+1} = animal{j};
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
        
        thetaHead{end+1,1} = theta{j,:};
        d_Theta{end+1,1,1} = dTheta{j,:};
        az{end+1,1}=azT{j,:};
        dist{end+1,1}=range{j,:};
        theta_good(end+1)=goodTheta(j);
        %         useData(end+1)=useData(j);
        %         approach{end+1}=appEpoch{j,:};
        %         Rcent{end+1,1}=x_centR{j,:}; Rcent{end,2}=y_centR{j,:};
        %         Lcent{end+1,1}=x_centL{j,:}; Lcent{end,2}=y_centL{j,:};
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
        
        Lgood{end+1,1}= Lngood{j,:};
        ccL(end+1)= Lcc(j);
        slopeL(end+1)= Lslope(j);
        scaleL(end+1)= Lscale(j);
        
        
    end
    
    
end
clear goodR goodL
goodR=ccR>.3;
goodL=ccL>.3;
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
    deltaR = diff(dist{useData(vid)})*30;
    badDist=(isnan(deltaR));
    vsmooth = conv(mouseSp{useData(vid)},ones(5,1)/5,'same');
    dRThresh=-10; %%%cm/sec
    vThresh=10;
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
clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    %subplot(rownum,colnum,vid);
    nframe = min(length(d_Theta{useData(vid)}),length(dRtheta{useData(vid)}));
    nframe = min(nframe, length(dLtheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe)==0
    dT=d_Theta{useData(vid)}; dtR=dRtheta{useData(vid)}; dtL=dLtheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    %  if sum(use)>3
    [corrR lagsR]= nanxcorr(dT(use),dtR(use),30,'coeff');
    %  plot(lagsR/30,corrR,'b');%xlim([-.3 .3])
    hold on;
    uselagsR=(lagsR>=-30& lagsR<=30);
    
    [corrL lagsL]= nanxcorr(dT(use),dtL(use),30,'coeff');
    % plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
    uselagsL=(lagsL>=-30 & lagsL<=30);
    %     else
    %     end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),30,'coeff');
        %  plot(lagsLA/30,corrLA,'c');%xlim([-.3 .3]);
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
subplot(1,2,1)
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]);
ylim([-.5 .5]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('d_head theta, both eyes');

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll
for vid=1:length(useData)
    %     nframe = min(length(dTheta{useData(vid)}),length(dthetaR{useData(vid)}));
    %     nframe = min(nframe, length(dthetaL{useData(vid)}));
    nonapp=appEpoch{vid}==0;
    dT=d_Theta{useData(vid)}; dpR=dRphi{useData(vid)}; dpL=dLphi{useData(vid)};
    clear use
    use = (nonapp==1)'% & ~isnan(dT(1:length(dpR)));
    if sum(use)>3
        [corrR lagsR]= nanxcorr(dT(use),dpR(use),30,'coeff');
        uselagsR=(lagsR>=-30& lagsR<=30);
        
        [corrL lagsL]= nanxcorr(dT(use),dpL(use),30,'coeff');
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1
    if sum(use)>3 &sum(~isnan(dpR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dpR(use),30,'coeff');
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dpL(use),30,'coeff');
        uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    else
    end
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        
    else
    end
end
% if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

subplot(1,2,2);
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

legend(L,{'dPhi R non-app','dPhi L non-app','dPhi R Approach','dPhi L Approach'}); title('d_head Theta and Eye Phi, both eyes');
 if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%
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

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.65 .65]);
xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dTheta non-app','dPhi non-app','dTh Approach','dPhi app'}); title('mean between eye corr, dtheta & dphi');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    nframe = min(length(az{useData(vid)}),length(dRtheta{useData(vid)}));
    nframe = min(nframe, length(dLtheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe)==0
    dT=az{useData(vid)}; dtR=dRtheta{useData(vid)}; dtL=dLtheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
        [corrR lagsR]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsR/30,corrR,'b');%xlim([-.3 .3])
        hold on;
        uselagsR=(lagsR>=-30& lagsR<=30);
        
        [corrL lagsL]= nanxcorr(dT(use),dtL(use),30,'coeff');
        % plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),30,'coeff');
        %  plot(lagsLA/30,corrLA,'c');%xlim([-.3 .3]);
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
subplot(1,2,1)
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.1 .1]); xlim([21 41]); axis square
% L(1) = plot(nan, nan, 'b-');
% L(2) = plot(nan, nan, 'r-');
% L(3) = plot(nan, nan, 'g-');
% L(4) = plot(nan, nan, 'c-');
%
% legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('az, eye theta both eyes');

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    nframe = min(length(diff(az{useData(vid)})),length(dRtheta{useData(vid)}));
    nframe = min(nframe, length(dLtheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe)==0
    dT=diff(az{useData(vid)}); dtR=dRtheta{useData(vid)}; dtL=dLtheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
        [corrR lagsR]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsR/30,corrR,'b');%xlim([-.3 .3])
        hold on;
        uselagsR=(lagsR>=-30& lagsR<=30);
        
        [corrL lagsL]= nanxcorr(dT(use),dtL(use),30,'coeff');
        % plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),30,'coeff');
        %  plot(lagsLA/30,corrLA,'c');%xlim([-.3 .3]);
        uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    else
    end
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        
    else
    end
end
subplot(1,2,2)
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.3 .3]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('d_az, eye theta both eyes');
 if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    nframe = min(length(az{useData(vid)}),length(Rtheta{useData(vid)}));
    nframe = min(nframe, length(Ltheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe-1)==0
    dT=az{useData(vid)}; dtR=Rtheta{useData(vid)}; dtL=Ltheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
        [corrR lagsR]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsR/30,corrR,'b');%xlim([-.3 .3])
        hold on;
        uselagsR=(lagsR>=-30& lagsR<=30);
        
        [corrL lagsL]= nanxcorr(dT(use),dtL(use),30,'coeff');
        % plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),30,'coeff');
        %  plot(lagsLA/30,corrLA,'c');%xlim([-.3 .3]);
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
subplot(1,2,1)
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

legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'}); title('az, eye theta pos both eyes');

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    nframe = min(length(diff(az{useData(vid)})),length(Rtheta{useData(vid)}));
    nframe = min(nframe, length(Ltheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe)==0
    dT=diff(az{useData(vid)}); dtR=Rtheta{useData(vid)}; dtL=Ltheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
        [corrR lagsR]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsR/30,corrR,'b');%xlim([-.3 .3])
        hold on;
        uselagsR=(lagsR>=-30& lagsR<=30);
        
        [corrL lagsL]= nanxcorr(dT(use),dtL(use),30,'coeff');
        % plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),30,'coeff');
        %  plot(lagsLA/30,corrLA,'c');%xlim([-.3 .3]);
        uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    else
    end
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        
    else
    end
end
subplot(1,2,2)
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.3 .3]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'}); title('d_az, eye theta position both eyes');
 if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

%%
for vid=1:length(useData)

nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
nframe = min(nframe, length(az{useData(vid)}));

% tR=Rtheta{useData(vid)}(5:end); tL=Ltheta{useData(vid)}(5:end); azC=(az{useData(vid)}(1:end-4));
tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe); azC=(az{useData(vid)}(1:nframe));
useN= appEpoch{vid}==0;

%useN= appEpoch{vid}(1:end-4)==0;
figure(1)
plot(azC(useN(1:15:end)),tR(useN(1:15:end)),'bo'); axis square; hold on;xlim([-pi pi]);
title('az, R eye Theta');
figure(2);
plot(azC(useN(1:15:end)),tL(useN(1:15:end)),'bo');axis square; hold on;xlim([-pi pi]);
title('az, L eye Theta');
end
for vid=1:length(useData)

nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
nframe = min(nframe, length(az{useData(vid)}));
tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe); azC=(az{useData(vid)}(1:nframe));
% tR=Rtheta{useData(vid)}(5:end); tL=Ltheta{useData(vid)}(5:end); azC=(az{useData(vid)}(1:end-4));
use= appEpoch{vid}==1;
% use= appEpoch{vid}(1:end-4)==1;

figure(1);
plot(azC(use(1:15:end)),tR(use(1:15:end)),'go');axis square; hold on; xlim([-pi pi]);
title('az, R eye Theta');
figure(2);
plot(azC(use(1:15:end)),tL(use(1:15:end)),'go');axis square; hold on;xlim([-pi pi]);
title('az, L eye Theta');
end





%%

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    nframe = min(length(cricketSp{useData(vid)}),length(dRtheta{useData(vid)}));
    nframe = min(nframe, length(dLtheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe)==0;
    dT=cricketSp{useData(vid)}; dtR=dRtheta{useData(vid)}; dtL=dLtheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
        [corrR lagsR]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsR/30,corrR,'b');%xlim([-.3 .3])
        hold on;
        uselagsR=(lagsR>=-30& lagsR<=30);
        
        [corrL lagsL]= nanxcorr(dT(use),dtL(use),30,'coeff');
        % plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),30,'coeff');
        %  plot(lagsLA/30,corrLA,'c');%xlim([-.3 .3]);
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
subplot(1,2,1)
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.1 .1]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('cricketSp, eye theta both eyes');

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    nframe = min(length(diff(cricketSp{useData(vid)})),length(dRtheta{useData(vid)}));
    nframe = min(nframe, length(dLtheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe)==0;
    dT=diff(cricketSp{useData(vid)}); dtR=dRtheta{useData(vid)}; dtL=dLtheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
        [corrR lagsR]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsR/30,corrR,'b');%xlim([-.3 .3])
        hold on;
        uselagsR=(lagsR>=-30& lagsR<=30);
        
        [corrL lagsL]= nanxcorr(dT(use),dtL(use),30,'coeff');
        % plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),30,'coeff');
        %  plot(lagsLA/30,corrLA,'c');%xlim([-.3 .3]);
        uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    else
    end
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        
    else
    end
end
subplot(1,2,2)
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.1 .1]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('d_cricketSp, eye theta both eyes');
 if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    nframe = min(length(cricketSp{useData(vid)}),length(Rtheta{useData(vid)}));
    nframe = min(nframe, length(Ltheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe-1)==0
    dT=cricketSp{useData(vid)}; dtR=Rtheta{useData(vid)}; dtL=Ltheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
        [corrR lagsR]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsR/30,corrR,'b');%xlim([-.3 .3])
        hold on;
        uselagsR=(lagsR>=-30& lagsR<=30);
        
        [corrL lagsL]= nanxcorr(dT(use),dtL(use),30,'coeff');
        % plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),30,'coeff');
        %  plot(lagsLA/30,corrLA,'c');%xlim([-.3 .3]);
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
subplot(1,2,1)
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.3 .3]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'}); title('cricketSp, eye theta pos both eyes');

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    nframe = min(length(diff(cricketSp{useData(vid)})),length(Rtheta{useData(vid)}));
    nframe = min(nframe, length(Ltheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe)==0
    dT=diff(cricketSp{useData(vid)}); dtR=Rtheta{useData(vid)}; dtL=Ltheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
        [corrR lagsR]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsR/30,corrR,'b');%xlim([-.3 .3])
        hold on;
        uselagsR=(lagsR>=-30& lagsR<=30);
        
        [corrL lagsL]= nanxcorr(dT(use),dtL(use),30,'coeff');
        % plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),30,'coeff');
        %  plot(lagsLA/30,corrLA,'c');%xlim([-.3 .3]);
        uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    else
    end
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        
    else
    end
end
subplot(1,2,2)
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.1 .1]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'}); title('d_cricketSp, eye theta position both eyes');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    nframe = min(length(mouseSp{useData(vid)}),length(dRtheta{useData(vid)}));
    nframe = min(nframe, length(dLtheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe)==0;
    dT=mouseSp{useData(vid)}; dtR=dRtheta{useData(vid)}; dtL=dLtheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
        [corrR lagsR]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsR/30,corrR,'b');%xlim([-.3 .3])
        hold on;
        uselagsR=(lagsR>=-30& lagsR<=30);
        
        [corrL lagsL]= nanxcorr(dT(use),dtL(use),30,'coeff');
        % plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),30,'coeff');
        %  plot(lagsLA/30,corrLA,'c');%xlim([-.3 .3]);
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
subplot(1,2,1)
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.1 .1]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');
%
legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('mouseSp, eye theta both eyes');

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    nframe = min(length(diff(mouseSp{useData(vid)})),length(dRtheta{useData(vid)}));
    nframe = min(nframe, length(dLtheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe)==0;
    dT=diff(mouseSp{useData(vid)}); dtR=dRtheta{useData(vid)}; dtL=dLtheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
        [corrR lagsR]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsR/30,corrR,'b');%xlim([-.3 .3])
        hold on;
        uselagsR=(lagsR>=-30& lagsR<=30);
        
        [corrL lagsL]= nanxcorr(dT(use),dtL(use),30,'coeff');
        % plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),30,'coeff');
        %  plot(lagsLA/30,corrLA,'c');%xlim([-.3 .3]);
        uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    else
    end
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        
    else
    end
end
subplot(1,2,2)
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.2 .2]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'dtheta R non-app','dtheta L non-app','dtheta R Approach','dtheta L Approach'}); title('d_mouseSp, eye theta both eyes');
 if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%%

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    nframe = min(length(mouseSp{useData(vid)}),length(Rtheta{useData(vid)}));
    nframe = min(nframe, length(Ltheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe-1)==0;
    dT=mouseSp{useData(vid)}; dtR=Rtheta{useData(vid)}; dtL=Ltheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
        [corrR lagsR]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsR/30,corrR,'b');%xlim([-.3 .3])
        hold on;
        uselagsR=(lagsR>=-30& lagsR<=30);
        
        [corrL lagsL]= nanxcorr(dT(use),dtL(use),30,'coeff');
        % plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),30,'coeff');
        %  plot(lagsLA/30,corrLA,'c');%xlim([-.3 .3]);
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
subplot(1,2,1)
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

legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'}); title('mouseSp, eye theta pos both eyes');

clear corrR lagsR corrRAll corrLAll uselagsL uselagsR corrL lagsL corrRAAll corrLAAll

for vid=1:length(useData)
    nframe = min(length(diff(mouseSp{useData(vid)})),length(Rtheta{useData(vid)}));
    nframe = min(nframe, length(Ltheta{useData(vid)}));
    nonapp=appEpoch{vid}(1:nframe)==0;
    dT=diff(mouseSp{useData(vid)}); dtR=Rtheta{useData(vid)}; dtL=Ltheta{useData(vid)};
    clear use
    use = (nonapp==1)'; %~isnan(dT(1:nframe)) &
    if sum(use)>3 & sum(~isnan(dtR(use)))>20 & sum(~isnan(dT(use)))>20
        [corrR lagsR]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsR/30,corrR,'b');%xlim([-.3 .3])
        hold on;
        uselagsR=(lagsR>=-30& lagsR<=30);
        
        [corrL lagsL]= nanxcorr(dT(use),dtL(use),30,'coeff');
        % plot(lagsL/30,corrL,'r');xlim([-.3 .3]);
        uselagsL=(lagsL>=-30 & lagsL<=30);
    else
    end
    clear use;
    use=appEpoch{vid}==1;
    if sum(use)>4 & sum(~isnan(dtR(use)))>20
        [corrRA lagsRA]= nanxcorr(dT(use),dtR(use),30,'coeff');
        %  plot(lagsRA/30,corrRA,'g');xlim([-.3 .3]);
        uselagsRA=(lagsRA>=-30& lagsRA<=30);
        [corrLA lagsLA]= nanxcorr(dT(use),dtL(use),30,'coeff');
        %  plot(lagsLA/30,corrLA,'c');%xlim([-.3 .3]);
        uselagsLA=(lagsLA>=-30 & lagsLA<=30);
    else
    end
    if sum(uselagsR)==61 & sum(uselagsL)==61 &sum(uselagsRA)==61 & sum(uselagsLA)==61
        corrRAll(vid,:)=corrR(uselagsR); corrLAll(vid,:)=corrL(uselagsL);
        corrRAAll(vid,:)=corrRA(uselagsRA); corrLAAll(vid,:)=corrLA(uselagsLA);
        
    else
    end
end
subplot(1,2,2)
errR= nanstd(corrRAll)/(sqrt(length(corrRAll))); errL = nanstd(corrLAll)/(sqrt(length(corrLAll)));
errRA= nanstd(corrRAAll)/(sqrt(length(corrRAAll))); errLA = nanstd(corrLAAll)/(sqrt(length(corrLAAll)));

shadedErrorBar(1:size(corrRAll,2),nanmean(corrRAll,1),errR,'-b',1); hold on
shadedErrorBar(1:size(corrLAll,2),nanmean(corrLAll,1),errL,'-r',1);
shadedErrorBar(1:size(corrRAAll,2),nanmean(corrRAAll,1),errRA,'-g',1); hold on
shadedErrorBar(1:size(corrLAAll,2),nanmean(corrLAAll,1),errLA,'-c',1);

plot([31,31],[1,-1],'--','Color', [.5 .5 .5]); ylim([-.2 .2]); xlim([21 41]); axis square
L(1) = plot(nan, nan, 'b-');
L(2) = plot(nan, nan, 'r-');
L(3) = plot(nan, nan, 'g-');
L(4) = plot(nan, nan, 'c-');

legend(L,{'theta R non-app','theta L non-app','theta R Approach','theta L Approach'}); title('d_mouseSp, eye theta position both eyes');
if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end

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

dTRAll=[];
close all
% figure%('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useData)
    clear dT dpR dpL dtR dtL
    %     subplot(rownum,colnum,vid);
    nframe = min(length(d_Theta{useData(vid)}),length(dRphi{useData(vid)}));
    nframe = min(nframe, length(dLphi{useData(vid)}));
    nframe = min(nframe, length(dRtheta{useData(vid)}));
    nframe = min(nframe,length(dLtheta{useData(vid)}));
    
    dT=d_Theta{useData(vid)}(1:nframe); pR=dRphi{useData(vid)}(1:nframe); pL=dLphi{useData(vid)}(1:nframe);
    tR=dRtheta{useData(vid)}(1:nframe); tL=dLtheta{useData(vid)}(1:nframe);
    useN= appEpoch{vid}==0; azC=az{useData(vid)}(1:nframe);
    use = (appEpoch{vid});

    
    figure(1);
    plot(dT(useN(1:15:end)),pR(useN(1:15:end)),'bo');axis square; hold on
    title('head d_Theta, phi R, nonapp & app');
    x = linspace(-40,40);
    y = linspace(-60,60);
    plot(-x,y); xlim([-40 40]); ylim([-60 60]);
    plot(dT(use(1:15:end)),pR(use(1:15:end)),'og');
    if vid==(useData(end))
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(2);
    plot(dT(useN(1:15:end)),pL(useN(1:15:end)),'bo'); axis square; hold on
    title('head d_Theta, phi L, nonapp & app');
    x = linspace(-40,40);
    y = linspace(-80,80);
    plot(x,y);  xlim([-40 40]); ylim([-80 80]);
    plot(dT(use(1:15:end)),pL(use(1:15:end)),'og');
    if vid==(useData(end))
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(3);
    plot(dT(useN(1:15:end)),tR(useN(1:15:end)),'bo');axis square; hold on
    title('head d_Theta, theta R, nonapp & app');
    x = linspace(-40,40);
    y = linspace(-40,40);  xlim([-40 40]); ylim([-40 40]);
    plot(-x,y);
    plot(dT(use(1:15:end)),tR(use(1:15:end)),'og');
    if vid==(useData(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(4);
    plot(dT(useN(1:15:end)),tL(useN(1:15:end)),'bo');axis square; hold on;
    title('head d_Theta, theta L, nonapp & app');
    x = linspace(-40,40);  xlim([-40 40]); ylim([-80 80]);
    y = linspace(-80,80);
    plot(x,y);
    plot(dT(use(1:15:end)),tL(use(1:15:end)),'og');
    if vid==(useData(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(5);%subplot(rownum,colnum,vid);
    plot(tR(useN(1:15:end)),tL(useN(1:15:end)),'bo');axis square; hold on;
    title('r theta, l theta, nonapp & app');
    x = linspace(-80,80);
    y = linspace(-50,50);
    plot(-x,y); xlim([-80 80]); ylim([-50 50]);
    plot(tR(use(1:15:end)),tL(use(1:15:end)),'og');
    if vid==(useData(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    figure(6);
    plot(pR(useN(1:15:end)),pL(useN(1:15:end)),'bo');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('right phi, left phi, nonapp & app');
    plot(pR(use(1:15:end)),pL(use(1:15:end)),'og');
    if vid==(useData(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
     figure(7);
    plot(tR(useN(1:15:end)),pR(useN(1:15:end)),'bo');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('right eye - theta/phi pos - nonapp & app');
    plot(tR(use(1:15:end)),pR(use(1:15:end)),'og');
    if vid==(useData(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
         figure(8);
    plot(tL(useN(1:15:end)),pL(useN(1:15:end)),'bo');axis square; hold on;
    x = linspace(-80,80);
    y = linspace(-80,80);
    plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
    title('left eye - theta/phi pos - nonapp & app');
    plot(tL(use(1:15:end)),pL(use(1:15:end)),'og');
    if vid==(useData(end))
        
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
    
end

%%

close all

% figure%('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useData)
    clear dT dpR dpL dtR dtL
    %     subplot(rownum,colnum,vid);
    nframe = min(length(dist{useData(vid)}),length(radR{useData(vid)}));
    nframe = min(nframe, length(radL{useData(vid)}));
    nframe = min(nframe, length(mouseSp{useData(vid)}));
    nframe = min(nframe,length(cricketSp{useData(vid)}));
    nframe = min(nframe,length(az{useData(vid)}));

    r=dist{useData(vid)}(1:nframe); 
    rR=radR{useData(vid)}(1:nframe);%rR = rR - min(rR(:)); rR = rR ./ max(rR(:)) 
    lR=radL{useData(vid)}(1:nframe);%lR = lR - min(lR(:)); lR = lR ./ max(lR(:)) 
    
    mouseV=mouseSp{useData(vid)}(1:nframe); crSp=cricketSp{useData(vid)}(1:nframe);
    azC=az{useData(vid)}(1:nframe);
    useN= appEpoch{vid}==0;
    use = (appEpoch{vid});
    figure(1);
    plot(rR(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on
    xlabel('r eye'); ylabel('l eye');
    title('two eyes, rad');
%     x = linspace(-40,40);
%     y = linspace(-60,60);
%     plot(-x,y); 
  %  xlim([0 1]); ylim([0 1]);
   % plot(rR(use(1:22:end)),lR(use(1:22:end)),'og');
   % plot(rR(use),lR(use),'og');


%     figure(2);
%     plot(r(useN(1:22:end)),rR(useN(1:22:end)),'bo'); axis square; hold on
%     title('range and R Rad');
%     xlabel('range to cricket (cm)'); ylabel('R eye');
% %     x = linspace(-40,40);
% %     y = linspace(-80,80);
% %     plot(x,y);  xlim([-40 40]); ylim([-80 80]);
% %     figure(3);
%     plot(r(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on
%     title('range and L Rad');
%     xlabel('range to cricket (cm)'); ylabel('L eye');
% %     x = linspace(-40,40);
% %     y = linspace(-40,40);  xlim([-40 40]); ylim([-40 40]);
% %     plot(-x,y);
%     figure(4);
%     plot(mouseV(useN(1:22:end)),rR(useN(1:22:end)),'bo');axis square; hold on;
%     title('mouse speed & R rad');
%     xlabel('mouse Speed (cm/sec)'); ylabel('R Rad');
% %     x = linspace(-40,40);  xlim([-40 40]); ylim([-80 80]);
% %     y = linspace(-80,80);
% %     plot(x,y);
%     plot(mouseV(use(1:22:end)),rR(use(1:22:end)),'og');
%    figure(5);%subplot(rownum,colnum,vid);
%     plot(mouseV(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on;
%     title('mouse speed & L rad');
%     xlabel('mouse Speed (cm/sec)'); ylabel('L Rad');
% %     x = linspace(-80,80);
% %     y = linspace(-50,50);
% %     plot(-x,y); xlim([-80 80]); ylim([-50 50]);
%     plot(mouseV(use(1:22:end)),lR(use(1:22:end)),'og');
%       figure(6);
%     plot(crSp(useN(1:22:end)),rR(useN(1:22:end)),'bo');axis square; hold on;
% %     x = linspace(-80,80);
% %     y = linspace(-80,80);
% %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
%     title('cricket speed & R rad');
%      xlabel('cricket Speed'); ylabel('R Rad');    
%     plot(crSp(use(1:22:end)),rR(use(1:22:end)),'og');
%     if vid==(useData(end))
%         
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     
%      figure(7);
%     plot(crSp(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on;
% % %     x = linspace(-80,80);
% % %     y = linspace(-80,80);
% % %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
% %   title('cricket speed & L rad');
%      xlabel('cricket Speed'); ylabel('l Rad');      
% %     plot(crSp(use(1:22:end)),lR(use(1:22:end)),'og');
% %     
%     figure(8);
%     plot(azC(useN(1:22:end)),r(useN(1:22:end)),'b.');axis square; hold on;
% % %     x = linspace(-80,80);
% % %     y = linspace(-80,80);
% % %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
%     title('az and range');
% %     plot(azC(use(1:22:end)),r(use(1:22:end)),'.g');
%        if vid==(useData(end))
%         
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
    
end



% figure%('units','normalized','outerposition',[0 0 1 1])
for vid=1:length(useData)
    clear dT dpR dpL dtR dtL
    %     subplot(rownum,colnum,vid);
    nframe = min(length(dist{useData(vid)}),length(radR{useData(vid)}));
    nframe = min(nframe, length(radL{useData(vid)}));
    nframe = min(nframe, length(mouseSp{useData(vid)}));
    nframe = min(nframe,length(cricketSp{useData(vid)}));
    nframe = min(nframe,length(az{useData(vid)}));

    r=dist{useData(vid)}(1:nframe); 
    rR=radR{useData(vid)}(1:nframe);%rR = rR - min(rR(:)); rR = rR ./ max(rR(:)) 
    lR=radL{useData(vid)}(1:nframe);%lR = lR - min(lR(:)); lR = lR ./ max(lR(:)) 
    
    mouseV=mouseSp{useData(vid)}(1:nframe); crSp=cricketSp{useData(vid)}(1:nframe);
    azC=az{useData(vid)}(1:nframe);
  %  useN= appEpoch{vid}==0;
    use = (appEpoch{vid});
    figure(1);
   % plot(rR(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on
    xlabel('r eye'); ylabel('l eye');
    title('two eyes, rad');
%     x = linspace(-40,40);
%     y = linspace(-60,60);
%     plot(-x,y); 
    xlim([15 46]); ylim([15 46]);
    plot(rR(use(1:22:end)),lR(use(1:22:end)),'og');
    if vid==(useData(end))
        if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
    end
   % figure(2);
%     plot(r(useN(1:22:end)),rR(useN(1:22:end)),'bo'); axis square; hold on
%     title('range and R Rad');
%     xlabel('range to cricket (cm)'); ylabel('R eye');
% %     x = linspace(-40,40);
% %     y = linspace(-80,80);
% %     plot(x,y);  xlim([-40 40]); ylim([-80 80]);
 %   plot(r(use(1:22:end)),rR(use(1:22:end)),'og');
%     if vid==(useData(end))
%         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
%     end
%     figure(3);
% %     plot(r(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on
% %     title('range and L Rad');
% %     xlabel('range to cricket (cm)'); ylabel('L eye');
% % %     x = linspace(-40,40);
% % %     y = linspace(-40,40);  xlim([-40 40]); ylim([-40 40]);
% % %     plot(-x,y);
%      plot(r(use(1:22:end)),lR(use(1:22:end)),'og');
% %     if vid==(useData(end))
% %         
% %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %     end
%      figure(4);
% %     plot(mouseV(useN(1:22:end)),rR(useN(1:22:end)),'bo');axis square; hold on;
% %     title('mouse speed & R rad');
% %     xlabel('mouse Speed (cm/sec)'); ylabel('R Rad');
% % %     x = linspace(-40,40);  xlim([-40 40]); ylim([-80 80]);
% % %     y = linspace(-80,80);
% % %     plot(x,y);
%     plot(mouseV(use(1:22:end)),rR(use(1:22:end)),'og');
% %     if vid==(useData(end))
% %         
% %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %     end
%      figure(5);%subplot(rownum,colnum,vid);
% %     plot(mouseV(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on;
% %     title('mouse speed & L rad');
% %     xlabel('mouse Speed (cm/sec)'); ylabel('L Rad');
% % %     x = linspace(-80,80);
% % %     y = linspace(-50,50);
% % %     plot(-x,y); xlim([-80 80]); ylim([-50 50]);
%      plot(mouseV(use(1:22:end)),lR(use(1:22:end)),'og');
% %     if vid==(useData(end))
% %         
% %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %     end
%      figure(6);
% %     plot(crSp(useN(1:22:end)),rR(useN(1:22:end)),'bo');axis square; hold on;
% % %     x = linspace(-80,80);
% % %     y = linspace(-80,80);
% % %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
% %     title('cricket speed & R rad');
% %     xlabel('cricket Speed'); ylabel('R Rad');    
%      plot(crSp(use(1:22:end)),rR(use(1:22:end)),'og');
% %     if vid==(useData(end))
% %         
% %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %     end
% %     
%       figure(7);
% %     plot(crSp(useN(1:22:end)),lR(useN(1:22:end)),'bo');axis square; hold on;
% % %     x = linspace(-80,80);
% % %     y = linspace(-80,80);
% % %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
% %   title('cricket speed & L rad');
% %     xlabel('cricket Speed'); ylabel('l Rad');      
%     plot(crSp(use(1:22:end)),lR(use(1:22:end)),'og');
% %        if vid==(useData(end))
% %         
% %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %     end
% %     
%     figure(8);
% %     plot(azC(useN(1:22:end)),r(useN(1:22:end)),'b.');axis square; hold on;
% % %     x = linspace(-80,80);
% % %     y = linspace(-80,80);
% % %     plot(-x,y);  xlim([-80 80]); ylim([-80 80]);
% %     title('az and range');
%      plot(azC(use(1:22:end)),r(use(1:22:end)),'.g');
% %        if vid==(useData(end))
% %         
% %         if savePDF, set(gcf, 'PaperPositionMode', 'auto');print('-bestfit','-dpsc',psfilename,'-append'); close(gcf); end
% %     end
    
end




%%
for vid=1:length(useData)

nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
nframe = min(nframe, length(az{useData(vid)}));

tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe); azC=(az{useData(vid)}(1:end-4));
pR=Rphi{useData(vid)}(1:nframe); pL=Lphi{useData(vid)}(1:nframe);
useN= appEpoch{vid}%(1:nframe)%(1:end-4)==0;
% figure(1)
% plot(tR(useN),tL(useN),'bo'); axis square; hold on;%xlim([-pi pi]);
% title('R & L theta');
figure(2);
plot(tR(useN),pR(useN),'bo');axis square; hold on;%xlim([-pi pi]);
title('R theta & phi');
figure(3);
plot(tL(useN),pL(useN),'bo');axis square; hold on;%xlim([-pi pi]);
title('L theta & phi');
end

for vid=1:length(useData)

nframe = min(length(Rtheta{useData(vid)}),length(Ltheta{useData(vid)}));
nframe = min(nframe, length(az{useData(vid)}));

tR=Rtheta{useData(vid)}(1:nframe); tL=Ltheta{useData(vid)}(1:nframe); azC=(az{useData(vid)}(1:end-4));
pR=Rphi{useData(vid)}(1:nframe); pL=Lphi{useData(vid)}(1:nframe);

use= appEpoch{vid}%(1:nframe)==1;
% figure(1)
% plot(tR(use),tL(use),'go'); axis square; hold on;%xlim([-pi pi]);
% title('R & L theta');
figure(2);
plot(tR(use),pR(use),'go');axis square; hold on;%xlim([-pi pi]);
title('R theta & phi');
figure(3);
plot(tL(use),pL(use),'go');axis square; hold on;%xlim([-pi pi]);
title('L theta & phi');
end





%%


if savePDF
    pSname='T:\PreyCaptureAnalysis\Data\';
    filen=sprintf('%s','Analyzed_090319_AllAnimals_CricketBody_c','.pdf')
    pdfilename=fullfile(pSname,filen);
    dos(['ps2pdf ' psfilename ' ' pdfilename]);
    delete(psfilename);
else
     pSname='T:\PreyCaptureAnalysis\Data\';
end


afilename=sprintf('%s','Analyzed_AllAnimals_090319_New_cricketBody_c','.mat')
save(fullfile(pSname, afilename))

