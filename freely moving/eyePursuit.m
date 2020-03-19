%%% create movies of pursuit with eye positions
clear all
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

% Data(vid).theta =interpNan(Data(vid).thetaRaw-nanmedian(Data(vid).thetaRaw),15,'linear')
Data(vid).Rtheta = interpNan(Data(vid).Rthetaraw - nanmedian(Data(vid).Rthetaraw),15,'linear');
Data(vid).Ltheta = interpNan(Data(vid).Lthetaraw - nanmedian(Data(vid).Lthetaraw),15,'linear');

Data(vid).Rphi = Data(vid).Rphiraw - nanmedian(Data(vid).Rphiraw);
Data(vid).Lphi = Data(vid).Lphiraw - nanmedian(Data(vid).Lphiraw);
mousexy = (Data(vid).mouse_xyRaw)*.25; %resize
Data(vid).cricketxy =(Data(vid).cricketxyRaw);

Data(vid).RTime=Data(vid).RTS;
Data(vid).LTime=Data(vid).LTS;
Data(vid).TopTime=Data(vid).TopTs;
% Data(vid).usedTime=Data(vid).usedTS %desired TS for interpolation, used in loadAllCsv, not needed here

headtheta=Data(vid).thetaRaw;
headtheta=circInterp(Data(vid).TopTime, headtheta, Data(vid).TopTime);
headTh =interpNan(headtheta-nanmedian(headtheta),15,'linear');
Rth =circInterp(Data(vid).RTime,Data(vid).Rtheta,Data(vid).TopTime);
Lth =circInterp(Data(vid).LTime,Data(vid).Ltheta,Data(vid).TopTime);

angles = 0:0.01:2*pi; %R= 5; %%% substitute actual eye radii
R=((Data(vid).RRad)/4);
RL=((Data(vid).LRad)/4);
startLag = 10; %%% offset from beginning (needed for trails behind mouse)

% if resized to .25, 1 cm = 6.75 px, 5 cm = 33.7500 pixels
% if resized to .25, 1 cm = 6.75 px, 5 cm = 33.7500 pixels
%%
figure
% vidObj = VideoWriter('J462a_111119_2_5_laserEyesGaze.avi');
% open(vidObj);
for i = (1+startLag)+3150:3550%length(headTh)
    %%% plot cricket and mouse tracks
    %     subplot(2,2,[1:2]);
    imshow(flipud((TopVid(:,:,i)))); hold on;
    %       imshow(flipud(TopVid(:,:,i))); hold on;
    % subplot(2,2,[1:2]);
    %        plot(Data(vid).cricketxy(1,i),Data(vid).cricketxy(2,i),'g*');
    %     plot(Data(vid).cricketxy(1,i)-20,Data(vid).cricketxy(2,i)-135,'g*');
    %      plot(Data(vid).mouse_xy(1,i + (-startLag:0)), Data(vid).mouse_xy(2,i + (-10:0)),'b', 'Linewidth',2)
    
    %%% calculate head vector
    hx = 20*cos(headTh(:,i));
    hy = 20*sin(headTh(:,i));
    plot((mousexy(1,i)+ [0 hx']),(mousexy(2,i)+ [0 hy'])-30,'w', 'Linewidth',.75); hold on
    %%% calculate gaze direction (head + eyes); assume each eye is centered
    %%% at 45deg (pi/4)

    rth = headTh(i) + pi/4 + Rth(i)*pi/180;
    lth = headTh(i) - pi/4 + Lth(i)*pi/180; 

    mn=.5*((headTh(i) + pi/4 + Rth(i)*pi/180)+(headTh(i) + pi/4 + Lth(i)*pi/180))
    gaze=mn+headTh(i);
    
    plot(((mousexy(1,i))) + [0 50*cos(lth)' ],((mousexy(2,i))) + [0 50*sin(lth)']-30,'Color', [.363 .586 .808], 'Linewidth',1)
    plot(((mousexy(1,i))) + [0 50*cos(rth)' ],((mousexy(2,i))) + [0 50*sin(rth)']-30,'Color',[.339 .238 .4336], 'Linewidth',1)

    hxg = 20*cos(gaze);
    hyg = 20*sin(gaze);
    plot((mousexy(1,i)+ [0 hxg]),(mousexy(2,i)+ [0 hyg])-30,'c', 'Linewidth',1); hold on
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
%     writeVideo(vidObj,currFrame);
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