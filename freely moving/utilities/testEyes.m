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
nframe = 600;

for frm = 1:nframe
    if round((frm-1)/25) == (frm-1)/25
        display(sprintf('done %d/%d frames',frm-1,nframe))
    end
    leftmov(:,:,:,frm)= topMov.readFrame;
end

[f p]= uigetfile('*.avi','right eye movie file');
topMov = VideoReader(fullfile(p,f));
nframe = round(topMov.duration*topMov.framerate)
nframe = 600
for frm = 1:nframe
    if round((frm-1)/25) == (frm-1)/25
        display(sprintf('done %d/%d frames',frm-1,nframe))
    end
    rightmov(:,:,:,frm)= topMov.readFrame;
end


figure
for i = 1:nframe
subplot(2,2,1)
imagesc((leftmov(:,:,:,i)))

subplot(2,2,2)
imagesc(fliplr(rightmov(:,:,:,i)))

subplot(2,2,3); hold off
plot(Data(j).Lthetaraw- nanmean(Data(j).Lthetaraw), Data(j).Lphiraw- nanmean(Data(j).Lphiraw));
hold on
plot(Data(j).Lthetaraw(i)- nanmean(Data(j).Lthetaraw), Data(j).Lphiraw(i)- nanmean(Data(j).Lphiraw),'ro');
axis([-60 60 -60 60])

subplot(2,2,4); hold off
plot(Data(j).Rthetaraw- nanmean(Data(j).Rthetaraw), Data(j).Rphiraw- nanmean(Data(j).Rphiraw));
hold on
plot(Data(j).Rthetaraw(i)- nanmean(Data(j).Rthetaraw), Data(j).Rphiraw(i)- nanmean(Data(j).Rphiraw),'ro');
axis([-60 60 -60 60])

mov(i) = getframe(gcf);
end

