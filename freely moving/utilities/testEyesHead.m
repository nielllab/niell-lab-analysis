clear all

j= 1;


[f p] = uigetfile('*.csv','left eye csv file');
Lfname = fullfile(p,f)


Data(j).DataL = (csvread(Lfname,3,0))

Data(j).xL= Data(j).DataL(:,2:3:end);   %%% flip left eye
Data(j).yL=Data(j).DataL(:,3:3:end);
Data(j).LLikelihood=Data(j).DataL(:,4:3:end);
[Data(j).Lthetaraw,Data(j).Lphiraw,Data(j).EllipseParamsL,Data(j).ExtraParamsL,Data(j).goodLeye] = EyeCameraCalc1(length(Data(j).xL(:,1)),Data(j).xL,Data(j).yL, Data(j).LLikelihood)
Data(j).XLcentraw=Data(j).EllipseParamsL(:,1);
Data(j).YLcentraw=Data(j).EllipseParamsL(:,2);
Data(j).Lphiraw = -Data(j).Lphiraw;   %%% reverse because images have 1 at top left corner.

[f p] = uigetfile('*.csv','right eye csv file');
Rfname = fullfile(p,f)

Data(j).DataR = (csvread(Rfname,3,0))
Data(j).xR=640 - Data(j).DataR(:,2:3:end);
Data(j).yR=Data(j).DataR(:,3:3:end);
Data(j).RLikelihood=Data(j).DataR(:,4:3:end);
[Data(j).Rthetaraw,Data(j).Rphiraw,Data(j).EllipseParamsR,Data(j).ExtraParamsR,Data(j).goodReye] = EyeCameraCalc1(length(Data(j).xR(:,1)), Data(j).xR,Data(j).yR, Data(j).RLikelihood)
Data(j).XRcentraw=Data(j).EllipseParamsR(:,1);  Data(j).YRcentraw=Data(j).EllipseParamsR(:,2);
Data(j).Rphiraw = -Data(j).Rphiraw;   %%% reverse because images have 1 at top left corner.


%
% [f p] = uigetfile('*.mat','data file');
% load(fullfile(p,f));

[f p]= uigetfile('*.avi','left eye movie file');
topMov = VideoReader(fullfile(p,f));
nframe = round(topMov.duration*topMov.framerate);
nframe = 300;

for frm = 1:nframe
    if round((frm-1)/25) == (frm-1)/25
        display(sprintf('done %d/%d frames',frm-1,nframe))
    end
    leftmov(:,:,:,frm)= topMov.readFrame;
end

[f p]= uigetfile('*.avi','right eye movie file');
topMov = VideoReader(fullfile(p,f));
nframe = round(topMov.duration*topMov.framerate)
nframe = 300
for frm = 1:nframe
    if round((frm-1)/25) == (frm-1)/25
        display(sprintf('done %d/%d frames',frm-1,nframe))
    end
    rightmov(:,:,:,frm)= topMov.readFrame;
end


[f p] = uigetfile('*.csv','top csv file');
fname = fullfile(p,f);

    psfilename = 'C:\analysisPS.ps';

aligned = alignHead(fname,8,0,psfilename,.90, .95)
Data(j).mouse_xyRaw=aligned.mouse_xy;
Data(j).mouseVRaw=aligned.mouseSp;
Data(j).thetaRaw=aligned.theta;
Data(j).dThetaRaw=aligned.dTheta;

Lth = Data(j).Lthetaraw- nanmean(Data(j).Lthetaraw);
Rth = Data(j).Rthetaraw- nanmean(Data(j).Rthetaraw);

[f p] = uigetfile('*.csv','right eye timestamps');
rTSfile = fullfile(p,f);
RTS = dlmread(rTSfile);
RTS= RTS(:,1)*60*60 + RTS(:,2)*60 + RTS(:,3);

[f p] = uigetfile('*.csv','left eye timestamps');
lTSfile = fullfile(p,f);
LTS = dlmread(lTSfile);
LTS= LTS(:,1)*60*60 + LTS(:,2)*60 + LTS(:,3);

[f p] = uigetfile('*.csv','top timestamps');
topTSfile = fullfile(p,f);
TopTs = dlmread(topTSfile);
TopTs= TopTs(:,1)*60*60 + TopTs(:,2)*60 + TopTs(:,3);

hth = circInterp(TopTs,Data(j).thetaRaw,RTS);
Lth = interp1(LTS,Lth,RTS);

eyeTh = 0.5*(Lth + Rth);

figure
nframe = 300;
clear mov

hthDeg = (hth-nanmean(hth(1:nframe)))*180/pi;

for i =1:nframe
    subplot(3,2,1)
    imagesc((leftmov(:,:,:,i))); axis equal;
    
    subplot(3,2,2)
    imagesc(fliplr(rightmov(:,:,:,i))); axis equal;
    
    % subplot(3,2,3); hold off
    % plot(Data(j).Lthetaraw- nanmean(Data(j).Lthetaraw), Data(j).Lphiraw- nanmean(Data(j).Lphiraw));
    % hold on
    % plot(Data(j).Lthetaraw(i)- nanmean(Data(j).Lthetaraw), Data(j).Lphiraw(i)- nanmean(Data(j).Lphiraw),'ro');
    % axis([-60 60 -60 60])
    
    subplot(3,2,3:4); hold off
    plot(eyeTh,'k','Linewidth',2); hold on; plot(Lth,'b');  plot(Rth,'r'); 
    plot(i,Lth(i),'bo'); plot(i,Rth(i),'ro');
    ylim([-30 30]); xlim([1 nframe])
  %  xlim([70 175]);
    
    subplot(3,2,5:6);hold off
    plot(hthDeg + eyeTh,'k','LineWidth',2); hold on; plot(hthDeg,'g'); plot(i,hthDeg(i),'go');
    ylim([-250 250]); legend('gaze','head');
    xlim([1 nframe]);
    mov(i) = getframe(gcf);
end

[f p] = uiputfile('*.avi','video out')
    movObj = VideoWriter(fullfile(p,f));
    open(movObj);
    writeVideo(movObj,mov);
    close(movObj);
