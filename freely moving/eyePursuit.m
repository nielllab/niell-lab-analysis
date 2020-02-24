%%% create movies of pursuit with eye positions
clear all
%814 16
[f p] = uigetfile('*.mat','data file');
load(fullfile(p,f));

[fT,pT] = uigetfile({'*.avi';},'choose Top video');
% Read in Eye Tracking Video
if f ~=0
    TempVidT = VideoReader(fullfile(pT,fT));
    frame=1; k=1;
 
    while hasFrame(TempVidT)
         TopVid(:,:,frame) = imresize(rgb2gray(readFrame(TempVidT)),.25);
%         TopVid(:,:,frame) = (rgb2gray(readFrame(TempVidT)));
 if mod(frame,100)==0
            fprintf('frame = %d\n',frame)
        end
        frame=frame+1;
    end
end
fprintf('Done Reading in Top Video \n');

vid = input('which trial? ');

%%% zero-center eye positions
% Data(vid).Rtheta = Data(vid).Rtheta - nanmedian(Data(vid).Rtheta);
% Data(vid).Ltheta = Data(vid).Ltheta - nanmedian(Data(vid).Ltheta);
% Data(vid).Rphi = Data(vid).Rphi - nanmedian(Data(vid).Rphi);
% Data(vid).Lphi = Data(vid).Lphi - nanmedian(Data(vid).Lphi);

Data(vid).Rtheta = interpNan(Data(vid).Rtheta - nanmedian(Data(vid).Rtheta),15,'linear');
Data(vid).Ltheta = interpNan(Data(vid).Ltheta - nanmedian(Data(vid).Ltheta),15,'linear');
Data(vid).Rphi = Data(vid).Rphi - nanmedian(Data(vid).Rphi);
Data(vid).Lphi = Data(vid).Lphi - nanmedian(Data(vid).Lphi);
Data(vid).mouse_xyN = (Data(vid).mouse_xy)*.25; %resize
Data(vid).cricketxy =(Data(vid).cricketxy);

angles = 0:0.01:2*pi; %R= 5; %%% substitute actual eye radii
R=((Data(vid).RRad)/4);
RL=((Data(vid).LRad)/4);
startLag = 10; %%% offset from beginning (needed for trails behind mouse)

% if resized to .25, 1 cm = 6.75 px, 5 cm = 33.7500 pixels
%%
figure
vidObj = VideoWriter('testVid_J470c_111519_1_2_b'); 
open(vidObj);
for i = (1+startLag):length(Data(vid).Rtheta)-79
    %%% plot left eye position (with circle)
%     subplot(2,2,3);
%     plot(Data(vid).Ltheta,Data(vid).Lphi); hold on
%     plot(Data(vid).Ltheta(i)+RL(i)*cos(angles),Data(vid).Lphi(i)+RL(i)*sin(angles),'Linewidth',2.5);
%     axis equal; axis([-60 60 -60 60]); hold off; title('left eye'); xlabel('theta'); ylabel('phi')
    %%% plot right eye position (with circle)
%     subplot(2,2,4);
%     plot(Data(vid).Rtheta,Data(vid).Rphi); hold on
%     plot(Data(vid).Rtheta(i)+R(i)*cos(angles),Data(vid).Rphi(i)+R(i)*sin(angles),'Linewidth',2.5);
%     axis equal; axis([-60 60 -60 60]); hold off; title('right eye'); xlabel('theta'); ylabel('phi')
    %%% plot cricket and mouse tracks
%     subplot(2,2,[1:2]);
     imshow(flipud((TopVid(:,:,i)))); hold on;
%       imshow(flipud(TopVid(:,:,i))); hold on;
   % subplot(2,2,[1:2]);
%        plot(Data(vid).cricketxy(1,i),Data(vid).cricketxy(2,i),'g*');
%     plot(Data(vid).cricketxy(1,i)-20,Data(vid).cricketxy(2,i)-135,'g*');
%      plot(Data(vid).mouse_xy(1,i + (-startLag:0)), Data(vid).mouse_xy(2,i + (-10:0)),'b', 'Linewidth',2)
 
    %%% calculate head vector

    %     hx = 20*cos(cumsum(accelData{vid}(i,6)));
    %     hy = 20*sin(cumsum(accelData{vid}(i,6)));
    hx = 20*cos(Data(vid).theta(:,i));
    hy = 20*sin(Data(vid).theta(:,i));
    plot((Data(vid).mouse_xyN(1,i)+ [0 hx']),(Data(vid).mouse_xyN(2,i)+ [0 hy']),'w', 'Linewidth',.75); hold on
    %%% calculate gaze direction (head + eyes); assume each eye is centered
    %%% at 45deg (pi/4)
%     rth = cumsum(accelData{vid}(i,6)) + pi/4 + Data(vid).Rtheta(i)*pi/180;
%     lth = cumsum(accelData{vid}(i,6)) - pi/4 + Data(vid).Ltheta(i)*pi/180;
        rth = Data(vid).theta(i) + pi/4 + Data(vid).Rtheta(i)*pi/180;
        lth = Data(vid).theta(i) - pi/4 + Data(vid).Ltheta(i)*pi/180;

    
  plot(((Data(vid).mouse_xyN(1,i))) + [0 50*cos(lth)' ],((Data(vid).mouse_xyN(2,i))) + [0 50*sin(lth)'],'m', 'Linewidth',1)
  plot(((Data(vid).mouse_xyN(1,i))) + [0 50*cos(rth)' ],((Data(vid).mouse_xyN(2,i))) + [0 50*sin(rth)'],'c', 'Linewidth',1)

    axis equal; hold off;
%     imshow(TopVid(:,:,i)); hold off
   drawnow
%     axis([0 1600 0 1200]);
%   axis([0 400 0 350]);
    mov(i-startLag) = getframe(gcf);
    if mod(i,100)==0
            fprintf('frame = %d\n',i)
    end
          currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end
close(vidObj);
%%
%%% save video


% 
% sprintf('TEST subj %s session %s date %s clipnum %s',Data(vid).ani{1},Data(vid).sessionnum{1}, Data(vid).date{1},Data(vid).clipnum{1})
% [f p] = uiputfile('*.avi','output video file');
% movObj = VideoWriter(fullfile(p,f));
% open(movObj);
% writeVideo(movObj,mov);
% close(movObj);
