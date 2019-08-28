%%% create movies of pursuit with eye positions
clear all

[f p] = uigetfile('*.mat','data file');
load(fullfile(p,f));
vid = input('which trial? ');

%%% zero-center eye positions
% Data(vid).Rtheta = Data(vid).Rtheta - nanmedian(Data(vid).Rtheta);
% Data(vid).Ltheta = Data(vid).Ltheta - nanmedian(Data(vid).Ltheta);
% Data(vid).Rphi = Data(vid).Rphi - nanmedian(Data(vid).Rphi);
% Data(vid).Lphi = Data(vid).Lphi - nanmedian(Data(vid).Lphi);


Data(vid).Rtheta = Data(vid).Rtheta - nanmedian(Data(vid).Rtheta);
Data(vid).Ltheta = Data(vid).Ltheta - nanmedian(Data(vid).Ltheta);
Data(vid).Rphi = Data(vid).Rphi - nanmedian(Data(vid).Rphi);
Data(vid).Lphi = Data(vid).Lphi - nanmedian(Data(vid).Lphi);

angles = 0:0.01:2*pi; R= 5; %%% substitute actual eye radii

startLag = 10; %%% offset from beginning (neededfor trails behind mouse)

figure
for i = (1+startLag):length(Data(vid).Rtheta)
    
    %%% plot left eye position (with circle)
    subplot(2,2,3);
    plot(Data(vid).Ltheta,Data(vid).Lphi); hold on
    plot(Data(vid).Ltheta(i)+R*cos(angles),Data(vid).Lphi(i)+R*sin(angles),'Linewidth',4);
    axis equal; axis([-60 60 -60 60]); hold off; title('left eye'); xlabel('theta'); ylabel('phi')
    
    %%% plot right eye position (with circle)
    subplot(2,2,4);
    plot(Data(vid).Rtheta,Data(vid).Rphi); hold on
    plot(Data(vid).Rtheta(i)+R*cos(angles),Data(vid).Rphi(i)+R*sin(angles),'Linewidth',4);
    axis equal; axis([-60 60 -60 60]); hold off; title('right eye'); xlabel('theta'); ylabel('phi')
    
    %%% plot cricket and mouse tracks
    subplot(2,2,[1:2]);
    plot(Data(vid).cricketxy(1,i),Data(vid).cricketxy(2,i),'g*');
    hold on
    plot(Data(vid).mouse_xy(1,i + (-startLag:0)), Data(vid).mouse_xy(2,i + (-10:0)),'b', 'Linewidth',2)
    
    %%% calculate head vector
    hx = 100*cos(Data(vid).theta(i));
    hy = 100*sin(Data(vid).theta(i));
    plot(Data(vid).mouse_xy(1,i) + [0 hx ] ,Data(vid).mouse_xy(2,i) + [0 hy],'k', 'Linewidth',2)
    
    %%% calculate gaze direction (head + eyes); assume each eye is centered
    %%% at 45deg (pi/4)
    rth = Data(vid).theta(i) + pi/4 + Data(vid).Rtheta(i)*pi/180;
    lth = Data(vid).theta(i) - pi/4 + Data(vid).Ltheta(i)*pi/180;
%     rth = Data(vid).theta(i) + Data(vid).Rtheta(i)*pi/180;
%     lth = Data(vid).theta(i) + Data(vid).Ltheta(i)*pi/180;
    
    plot(Data(vid).mouse_xy(1,i) + [0 200*cos(rth) ] ,Data(vid).mouse_xy(2,i) + [0 200*sin(rth)],'c', 'Linewidth',2)
    plot(Data(vid).mouse_xy(1,i) + [0 200*cos(lth) ] ,Data(vid).mouse_xy(2,i) + [0 200*sin(lth)],'m', 'Linewidth',2)
    
    axis equal; axis([0 1600 0 1200]);  hold off
    
    drawnow
    mov(i-startLag) = getframe(gcf);
end


%%% save video
sprintf('subj %s session %s date %s clipnum %s',Data(vid).ani{1},Data(vid).sessionnum{1}, Data(vid).date{1},Data(vid).clipnum{1})
[f p] = uiputfile('*.avi','output video file');
movObj = VideoWriter(fullfile(p,f));
open(movObj);
writeVideo(movObj,mov);
close(movObj);
