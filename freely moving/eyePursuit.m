%%% create movies of pursuit with eye positions
clear all; close all

[f p] = uigetfile('*.mat','data file');
load(fullfile(p,f));

[f,p] = uigetfile({'*.avi';},'choose Left eye Vid');
% Read in Eye Tracking Video
if f ~=0
    TempVid = VideoReader(fullfile(p,f));
    frame=1; k=1;
    while hasFrame(TempVid) 
        EyeVidL(:,:,frame) = rgb2gray(readFrame(TempVid));
        if mod(frame,100)==0
            fprintf('frame = %d\n',frame)
        end
        frame=frame+1;
        
    end

    end
fprintf('Done Reading in Left Video \n');
  [f,p] = uigetfile({'*.avi';},'choose Right eye vid');
  
  
% Read in Eye Tracking Video
if f ~=0
    TempVidR = VideoReader(fullfile(p,f));   
    frameR=1; k=1;
    while hasFrame(TempVidR) 
        EyeVidR(:,:,frameR) = rgb2gray(readFrame(TempVidR));
        if mod(frameR,100)==0
            fprintf('frame = %d\n',frameR)
        end
        frameR=frameR+1;
        
    end  
end
fprintf('Done Reading in Right Video \n');


[fT,pT] = uigetfile({'*.avi';},'choose Top video');
videoScale = 1%0.25;
% Read in Eye Tracking Video
if f ~=0
    TempVidT = VideoReader(fullfile(pT,fT));
    frame=1; k=1;
    
    while hasFrame(TempVidT)
%         TopVid(:,:,frame) = imresize(rgb2gray(readFrame(TempVidT)),videoScale);
                 TopVid(:,:,frame) = (rgb2gray(readFrame(TempVidT)));
        if mod(frame,100)==0
            fprintf('frame = %d\n',frame)
        end
        frame=frame+1;
    end
end
fprintf('Done Reading in Top Video \n');
figure
firstFrame = double(TopVid(:,:,1));
hist(firstFrame(:))
videoGain = floor(255/prctile(firstFrame(:),99))

%23
vid = input('which trial? ');

% Data(vid).theta =interpNan(Data(vid).thetaRaw-nanmedian(Data(vid).thetaRaw),15,'linear')
rthRaw = Data(vid).Rthetaraw  - nanmedian(Data(vid).Rthetaraw);
lthRaw = Data(vid).Lthetaraw - nanmedian(Data(vid).Lthetaraw);
rthRaw(abs(rthRaw)>30) = NaN;  %%% remove extreme values
lthRaw(abs(lthRaw)>30) = NaN;

Data(vid).Rtheta = interpNan(rthRaw,15,'linear');
Data(vid).Ltheta = interpNan(lthRaw ,15,'linear');

figure

subplot(2,1,1); plot(Data(vid).Rtheta); hold on; plot(rthRaw); title('rth interp'); legend('nan')
subplot(2,1,2); plot(Data(vid).Ltheta); hold on; plot(lthRaw); title('lth interp'); legend('nan')

% Data(vid).Rphi = Data(vid).Rphiraw - nanmedian(Data(vid).Rphiraw);
% Data(vid).Lphi = Data(vid).Lphiraw - nanmedian(Data(vid).Lphiraw);

mousexy = Data(vid).mouse_xyRaw * videoScale; %resize
cricketxy =Data(vid).cricketxyRaw * videoScale;

RTime=Data(vid).RTS;
LTime=Data(vid).LTS;
TopTime=Data(vid).TopTs;
% Data(vid).usedTime=Data(vid).usedTS %desired TS for interpolation, used in loadAllCsv, not needed here

azCorrection = -6 *pi/180; %%% factor to correct head angles,based on DLC fit to head (this is actually calculated in makeFigs, but value is -6deg)
headtheta=Data(vid).thetaRaw - azCorrection;

figure
plot(headtheta)


%%% don't subtract meadian of head (since it's real-world coordinate)
headTh =interpNan(headtheta,15,'linear');
bad = find(abs(diff(headTh))>pi/4);
for i = 1:length(bad)
    headTh(bad(i))= headTh(bad(i)-1);
end
figure
plot(headTh,'.'); hold on; plot(headtheta); title('head th'); legend('nan')

Rth =interp1(RTime,Data(vid).Rtheta,TopTime)*pi/180;
Lth =interp1(LTime,Data(vid).Ltheta,TopTime)*pi/180;

%%% extrapolate eye videos
LVid=im2double(EyeVidL);
leftVidShift = shiftdim(LVid,2); %%% moves 2 dimensions to left
leftVidShiftInterp= interp1(LTime,leftVidShift,TopTime,'linear');
LVidInterp = shiftdim(leftVidShiftInterp,1); %%% move it one more to the left to get back in place.

RVid=im2double(EyeVidR);
rightVidShift = shiftdim(RVid,2); %%% moves 2 dimensions to left
rightVidShiftInterp= interp1(RTime,rightVidShift,TopTime,'linear');
RVidInterp = shiftdim(rightVidShiftInterp,1); %%% move it one more to the left to get back in place.

figure
subplot(3,1,1); plot(headTh); title('head th')
subplot(3,1,2); plot(Rth); title('rth');
subplot(3,1,3); plot(Lth); title('lth');

angles = 0:0.01:2*pi; %R= 5; %%% substitute actual eye radii
R=((Data(vid).RRad)/4);
RL=((Data(vid).LRad)/4);
startLag = 10; %%% offset from beginning (needed for trails behind mouse)

% if resized to .25, 1 cm = 6.75 px, 5 cm = 33.7500 pixels
% if resized to .25, 1 cm = 6.75 px, 5 cm = 33.7500 pixels

%%
figure
vidObj = VideoWriter('J462a_111119_1_4_Supp2_033120.avi');
VidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
eyeSz = 0;
headSz = 165;
gazeSz = 165;
showCricket = 0;
lineW = 2;
%for i = (1+startLag)+3150:3550%length(headTh) 
%for i= (300+startLag):length(headTh) %for supplemental movie 1, clip at 10
%seconds
for i=(90+startLag):240 %for supplemental movie 2, use 3-8 seconds
    %%% plot cricket and mouse tracks
    imshow(videoGain*flipud((TopVid(:,:,i)))); hold on;
    %%% calculate head vector
    hx = headSz*cos(headTh(:,i));
    hy = headSz*sin(headTh(:,i));
    plot((mousexy(1,i)+ [0 hx']),(mousexy(2,i)+ [0 hy'])-120,'Color', [.07 .226 .59], 'Linewidth',lineW+.3); hold on

    %%% calculate gaze direction (head + eyes); assume each eye is centered
    %%% at 45deg (pi/4)

    rth = headTh(i) + pi/4 + Rth(i);
    lth = headTh(i) - pi/4 + Lth(i); 

%     mn=.5*((headTh(i) + pi/4 + Rth(i)*pi/180)+(headTh(i) + pi/4 + Lth(i)*pi/180))
%     gaze=mn+headTh(i);
    gaze = headTh(i) + 0.5*(Rth(i) + Lth(i));

    %%% plotleft and right eyes
    plot(((mousexy(1,i))) + [0 eyeSz*cos(lth)' ],((mousexy(2,i))) + [0 eyeSz*sin(lth)']-120,'Color', [.43 .62 .9], 'Linewidth',lineW+.3)%.363 .586 .808
    plot(((mousexy(1,i))) + [0 eyeSz*cos(rth)' ],((mousexy(2,i))) + [0 eyeSz*sin(rth)']-120,'Color',[.5 .238 .52], 'Linewidth',lineW+.3)% .339 .238 .4336

    %%% plot gaze
    hxg = gazeSz*cos(gaze);
    hyg = gazeSz*sin(gaze);
    plot((mousexy(1,i)+ [0 hxg]),(mousexy(2,i)+ [0 hyg])-120,'c', 'Linewidth',lineW+.3); hold on
    %%% plot cricket
   if showCricket
       plot(cricketxy(1,i),cricketxy(2,i)-120,'b.');
   end
   
    axis equal; hold off
    %     imshow(TopVid(:,:,i)); hold off
    drawnow
    %     axis([0 1600 0 1200]);
    %   axis([0 400 0 350]);
%     subplot(3,2,5);
%     imshow(LVidInt(:,:,i));hold off; drawnow
%     subplot(3,2,6);
%     imshow(RVidInt(:,:,i)); hold off; drawnow
%     
    mov(i-startLag) = getframe(gcf);
    if mod(i,100)==0
        fprintf('frame = %d\n',i)
    end
    currFrame = getframe(gcf);
   writeVideo(vidObj,currFrame);
end
close(vidObj);

%%

figure
vidObj = VideoWriter('J462c_111119_2_6_Supp1_noHead_LEye.avi');
VidObj.Quality = 100;
vidObj.FrameRate = 10;
open(vidObj);
eyeSz = 0;
lineW = 1.75;
for i= (90+startLag):240%length(headTh)
    imshow(LVidInterp(:,:,i));hold off; drawnow
    mov(i-startLag) = getframe(gcf);
    if mod(i,100)==0
        fprintf('frame = %d\n',i)
    end
    currFrame = getframe(gcf);
   writeVideo(vidObj,currFrame);
end
close(vidObj)


figure
vidObj = VideoWriter('J462c_111119_2_6_Supp1_noHead_REye.avi');
VidObj.Quality = 100;
vidObj.FrameRate = 30;
open(vidObj);
eyeSz = 10;
lineW = 1.75;
for i= (90+startLag):240%length(headTh)
    imshow(RVidInterp(:,:,i));hold off; drawnow
    mov(i-startLag) = getframe(gcf);
    if mod(i,100)==0
        fprintf('frame = %d\n',i)
    end
    currFrame = getframe(gcf);
   writeVideo(vidObj,currFrame);
end
close(vidObj)



%%
%%% save video


%
% sprintf('TEST subj %s session %s date %s clipnum %s',Data(vid).ani{1},Data(vid).sessionnum{1}, Data(vid).date{1},Data(vid).clipnum{1})
% [f p] = uiputfile('*.avi','output video file');
% movObj = VideoWriter(fullfile(p,f));
% open(movObj);
% writeVideo(movObj,mov);
% close(movObj);