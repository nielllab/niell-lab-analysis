%%% create movies of pursuit with eye positions
clear all

% [f p] = uigetfile('*.mat','data file');
% load(fullfile(p,f));
load('CK2-CK2-7P1-RT_AllSessions_112719.mat')


vid = input('which trial? ');

sprintf('%s %s %s',Data(vid).ani{1},Data(vid).date{i},Data(vid).clipnum{i});

[fT,pT] = uigetfile({'*.avi';},'choose side video');
% Read in Eye Tracking Video

TempVidT = VideoReader(fullfile(pT,fT));
frame=1; k=1;

while hasFrame(TempVidT)
    %         TopVid(:,:,frame) = imresize(rgb2gray(readFrame(TempVidT)),.25);
    TopVid(:,:,frame) = (rgb2gray(readFrame(TempVidT)));
    
    if mod(frame,100)==0
        fprintf('frame = %d\n',frame)
    end
    frame=frame+1;
    
end
fprintf('Done Reading in Top Video \n');

R(:,1) = Data(vid).Rtheta - nanmedian(Data(vid).Rtheta);
R(:,2) = Data(vid).Rphi - nanmedian(Data(vid).Rphi);
L(:,1) = Data(vid).Ltheta - nanmedian(Data(vid).Ltheta);
L(:,2) = Data(vid).Lphi - nanmedian(Data(vid).Lphi);
R = R*Data(vid).scaleR/50;
L = L*Data(vid).scaleL/50;


%%% zero-center head theta, and get rid of wrap-around effect (mod 360)
%%% add pi/8 since this is roughly head tilt in movies relative to mean theta
th = Data(vid).theta- nanmedian(Data(vid).theta) + pi + pi/8;


Data(vid).mouse_xy(2,:) = Data(vid).mouse_xy(2,:)-640;

R = interp1(Data(vid).usedTS,R,Data(vid).TopTs);
L = interp1(Data(vid).usedTS,L,Data(vid).TopTs);
Data(vid).mouse_xy = interp1(Data(vid).usedTS,Data(vid).mouse_xy',Data(vid).TopTs)';
th = interp1(Data(vid).usedTS,th,Data(vid).TopTs);
%%% div = eye divergence (theta)
div = 0.5*(R(:,1)-L(:,1));
%%% gaze th = mean theta of eyes
gaze_th = (R(:,1) + L(:,1))*0.5;
%%% gaze phi = mean phi of eyes
gaze_phi = (R(:,2) + L(:,2))*0.5;


startLag = 10; %%% offset from beginning (neededfor trails behind mouse)

lag = 3; %%% fix slight offset between eye and top cams
figure
for i = (1+startLag):length(R)-lag
   %%% display image (flipped into xy coords0   
    imshow(flipud(TopVid(:,:,i))); hold on;
    
    %%% plot mouse head position with "tracers"
    plot(Data(vid).mouse_xy(1,i + (-startLag:0)), Data(vid).mouse_xy(2,i + (-10:0)),'b', 'Linewidth',2)
    plot(Data(vid).mouse_xy(1,i), Data(vid).mouse_xy(2,i) ,'bo', 'Markersize',8)
    
    %%% calculate and plot head vector
    hx = 200*cos(th(i));
    hy = 200*sin(th(i));
    plot((Data(vid).mouse_xy(1,i) + [0 hx ]) ,(Data(vid).mouse_xy(2,i) + [0 hy]),'k', 'Linewidth',2)
    
    %%% calculate gaze direction (head + eyes);  
    %%%subtract off the pi/8 that was added above
    rth = th(i) - div(i+lag)*pi/180 - pi/8;
    plot(Data(vid).mouse_xy(1,i) + [0 200*cos(rth) ] ,Data(vid).mouse_xy(2,i) + [0 200*sin(rth)],'c', 'Linewidth',1)
   
    %%% set axes and grab movie frame
    axis equal; hold off; axis xy
    mov(i-startLag) = getframe(gcf);
    if mod(i,100)==0
        fprintf('frame = %d\n',i)
    end
end


%%% save video
sprintf('subj %s session %s date %s clipnum %s',Data(vid).ani{1}, Data(vid).date{1},Data(vid).clipnum{1})
[f p] = uiputfile('*.avi','output video file');
movObj = VideoWriter(fullfile(p,f));
movObj.FrameRate = 50;
open(movObj);
writeVideo(movObj,mov);
close(movObj);
