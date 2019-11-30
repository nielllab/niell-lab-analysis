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
%         TopVid(:,:,frame) = imresize(rgb2gray(readFrame(TempVidT)),.25);
        TopVid(:,:,frame) = (rgb2gray(readFrame(TempVidT)));

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
Data(vid).mouse_xy = (Data(vid).mouse_xy);
Data(vid).cricketxy =(Data(vid).cricketxy+20);


angles = 0:0.01:2*pi; %R= 5; %%% substitute actual eye radii
R=((Data(vid).RRad)/4);
RL=((Data(vid).LRad)/4);
startLag = 10; %%% offset from beginning (neededfor trails behind mouse)


% if resized to .25, 1 cm = 6.75 px, 5 cm = 33.7500 pixels
figure
for i = (1+startLag):length(Data(vid).Rtheta)

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
     imshow(flipud(TopVid(:,:,i))); hold on;
   % subplot(2,2,[1:2]);
    plot(Data(vid).cricketxy(1,i)-20,Data(vid).cricketxy(2,i)-135,'g*');
    hold on
    plot(Data(vid).mouse_xy(1,i + (-startLag:0))-20, Data(vid).mouse_xy(2,i + (-10:0))-115,'b', 'Linewidth',2)
 
    %%% calculate head vector
    hx = 20*cos(Data(vid).theta(i));
    hy = 20*sin(Data(vid).theta(i));
    plot((Data(vid).mouse_xy(1,i) + [0 hx ])-20 ,(Data(vid).mouse_xy(2,i) + [0 hy])-115,'k', 'Linewidth',2)
 
    %%% calculate gaze direction (head + eyes); assume each eye is centered
    %%% at 45deg (pi/4)
    rth = Data(vid).theta(i) + pi/4 + Data(vid).Rtheta(i)*pi/180;
    lth = Data(vid).theta(i) - pi/4 + Data(vid).Ltheta(i)*pi/180;
%    rth = Data(vid).theta(i) + Data(vid).Rtheta(i)*pi/180;
%     lth = Data(vid).theta(i) + Data(vid).Ltheta(i)*pi/180;
    
    plot(Data(vid).mouse_xy(1,i) + [0 200*cos(rth) ]-20 ,Data(vid).mouse_xy(2,i) + [0 200*sin(rth)]-115,'c', 'Linewidth',1)
    plot(Data(vid).mouse_xy(1,i) + [0 200*cos(lth) ] -20,Data(vid).mouse_xy(2,i) + [0 200*sin(lth)]-115,'m', 'Linewidth',1)
    axis equal; hold off;
 %   imshow(TopVid(:,:,i)); hold off
   % drawnow
    axis([0 1600 0 1200]);
%   axis([0 400 0 350]);
    mov(i-startLag) = getframe(gcf);
    if mod(i,100)==0
            fprintf('frame = %d\n',i)
        end
end


%%% save video
sprintf('subj %s session %s date %s clipnum %s',Data(vid).ani{1},Data(vid).sessionnum{1}, Data(vid).date{1},Data(vid).clipnum{1})
[f p] = uiputfile('*.avi','output video file');
movObj = VideoWriter(fullfile(p,f));
open(movObj);
writeVideo(movObj,mov);
close(movObj);
