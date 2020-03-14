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

Data(vid).theta =interpNan(Data(vid).thetaRaw-nanmedian(Data(vid).thetaRaw),15,'linear')
Data(vid).Rtheta = interpNan(Data(vid).Rthetaraw - nanmedian(Data(vid).Rthetaraw),15,'linear');
Data(vid).Ltheta = interpNan(Data(vid).Lthetaraw - nanmedian(Data(vid).Lthetaraw),15,'linear');
Data(vid).Rphi = Data(vid).Rphiraw - nanmedian(Data(vid).Rphiraw);
Data(vid).Lphi = Data(vid).Lphiraw - nanmedian(Data(vid).Lphiraw);
Data(vid).mouse_xyN = (Data(vid).mouse_xyRaw)*.25; %resize
Data(vid).cricketxy =(Data(vid).cricketxyRaw);

Data(vid).RTime =Data(vid).RTS;
Data(vid).LTime=Data(vid).LTS;
Data(vid).TopTime=Data(vid).TopTs;
Data(vid).usedTime=Data(vid).usedTS %desired TS for interpolation

headTh =Data(vid).theta;
Rth =interp1(Data(vid).RTime,Data(vid).Rtheta,Data(vid).TopTime);
Lth =interp1(Data(vid).LTime,Data(vid).Ltheta,Data(vid).TopTime);
mousexy=Data(vid).mouse_xyN;
% Top=interp3(Data(vid).TopTime,TopVid,Data(vid).usedTime);

angles = 0:0.01:2*pi; %R= 5; %%% substitute actual eye radii
R=((Data(vid).RRad)/4);
RL=((Data(vid).LRad)/4);
startLag = 10; %%% offset from beginning (needed for trails behind mouse)

% if resized to .25, 1 cm = 6.75 px, 5 cm = 33.7500 pixels
%%
figure
vidObj = VideoWriter('test_TRIAL19.avi');
open(vidObj);
for i = (1+startLag)+3150:3550%length(headTh)
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
    hx = 20*cos(headTh(:,i));
    hy = 20*sin(headTh(:,i));
    plot((mousexy(1,i)+ [0 hx']),(mousexy(2,i)+ [0 hy'])-30,'w', 'Linewidth',.75); hold on
    %%% calculate gaze direction (head + eyes); assume each eye is centered
    %%% at 45deg (pi/4)
    %     rth = cumsum(accelData{vid}(i,6)) + pi/4 + Data(vid).Rtheta(i)*pi/180;
    %     lth = cumsum(accelData{vid}(i,6)) - pi/4 + Data(vid).Ltheta(i)*pi/180;
    rth = headTh(i) + pi/4 + Rth(i)*pi/180;
    lth = headTh(i) - pi/4 + Lth(i)*pi/180;
    
    
    plot(((mousexy(1,i))) + [0 50*cos(lth)' ],((mousexy(2,i))) + [0 50*sin(lth)']-30,'m', 'Linewidth',1)
    plot(((mousexy(1,i))) + [0 50*cos(rth)' ],((mousexy(2,i))) + [0 50*sin(rth)']-30,'c', 'Linewidth',1)
    
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
