%%% quick analysis of jumping/eye data
%%% mostly just cross-correlations
%%% cmn 2019

close all
clear all
load('CK2-CK2-7P1-RT_AllSessions_112719.mat')
allTheta = []; allPhi = []; allDiv = [];

%%% setup for pdf
savePDF = 1;
psfilename = 'C:\analysisPS.ps';
if exist(psfilename,'file')==2; delete(psfilename);end

%%% loop over videas
for i= 1:length(Data);
    
    clear R L
    %%% zero-center eye data (can also edit this to select pixels vs angles from fit    
    R(:,1) = Data(i).Rtheta - nanmedian(Data(i).Rtheta);
    R(:,2) = Data(i).Rphi - nanmedian(Data(i).Rphi);
    L(:,1) = Data(i).Ltheta - nanmedian(Data(i).Ltheta);
    L(:,2) = Data(i).Lphi - nanmedian(Data(i).Lphi);
    
    R = R*Data(i).scaleR/50;
    L = L*Data(i).scaleL/50;
    
    %%% zero-center head theta, and get rid of wrap-around effect (mod 360)
    th = Data(i).theta;
    th = th*180/pi; th = mod(th+360,360); th = th- nanmean(th); th = -th;
    
    %%% div = eye divergence (theta)
    div = 0.5*(R(:,1)-L(:,1));
    %%% gaze th = mean theta of eyes
    gaze_th = (R(:,1) + L(:,1))*0.5;
    %%% gaze phi = mean phi of eyes
    gaze_phi = (R(:,2) + L(:,2))*0.5;
    
    figure
    
    %%% plot traces
    subplot(2,3,1:3);
    hold on
    plot(th);
    % plot(R(:,1)); plot(L(:,1));
    plot(gaze_th); plot(div); plot(gaze_phi);
    title(sprintf('%d %s %s %s R=%0.2f L=%0.2f',i,Data(i).ani{1},Data(i).date{1},Data(i).clipnum{1}, Data(i).scaleR,Data(i).scaleL));
    legend('head','eye th','eye div','eye phi');
    ylabel('deg'); xlabel('frames');
    
    %%% calculate and plot xcorrs
    subplot(2,3,4);
    
    [th_gz(:,i) lags] = nanxcorr(th,gaze_th,30,'coeff');
    hold on
    plot(lags,th_gz(:,i)); ylim([-1 1]);
    
    [th_div(:,i) lags] = nanxcorr(th,div,30,'coeff');
    plot(lags,th_div(:,i)); ylim([-1 1]);
    
    [th_phi(:,i) lags] = nanxcorr(th, gaze_phi,30,'coeff');
    plot(lags,th_phi(:,i)); ylim([-1 1]); title('head theta vs phi');
    title('head theta xcorr');
   if i ==1, legend('gaze','div','phi');end
    
    %%% scatter plots
    subplot(2,3,5);
    plot(th, div,'.'); axis square; axis([-40 40 -40 40]); hold on; plot([-40 40],[40 -40],'r:');
    xlabel('head th deg'); ylabel('eye div deg');
    
    subplot(2,3,6);
    plot(th,gaze_phi,'.'); axis square; axis([-40 40 -40 40]); hold on; plot([-40 40],[-40 40],'r:');
    xlabel('head th deg'); ylabel('eye phi deg');
    
    allTheta = [allTheta th];
    allPhi = [allPhi gaze_phi'];
    allDiv = [allDiv div'];
    
    if savePDF,  set(gcf, 'PaperPositionMode', 'auto'); print('-dpsc',psfilename,'-append');  end
    
end

%%% plot pooled data
figure
plot(allTheta,allPhi,'.'); xlabel('head theta'); ylabel('phi');  axis square; axis([-60 60 -60 60])
if savePDF,  set(gcf, 'PaperPositionMode', 'auto'); print('-dpsc',psfilename,'-append');  end

figure
plot(allTheta,allDiv,'.'); xlabel('head theta'); ylabel('eye theta div'); axis  square; axis([-60 60 -60 60]);
if savePDF,  set(gcf, 'PaperPositionMode', 'auto'); print('-dpsc',psfilename,'-append');  end


figure
hold on
errorbar(lags,nanmean(th_gz,2),std(th_gz,[],2)/sqrt(size(th_gz,2)));
errorbar(lags,nanmean(th_div,2),std(th_div,[],2)/sqrt(size(th_gz,2)));
errorbar(lags,nanmean(th_phi,2),std(th_phi,[],2)/sqrt(size(th_gz,2)));
ylim([-1 1]); ylabel('correlation');
title('xcorr with head angle'); legend('mean theta','theta divergence','mean phi');
if savePDF,  set(gcf, 'PaperPositionMode', 'auto'); print('-dpsc',psfilename,'-append');  end

if savePDF
    [f p] = uiputfile('*.pdf','pdf file');
    pdfilename = fullfile(p,f);
    dos(['ps2pdf ' psfilename ' ' pdfilename]);
end
