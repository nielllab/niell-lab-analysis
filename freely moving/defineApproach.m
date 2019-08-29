close all
framerate =30;

vid = 1
figure
plot(mouse_xy{vid}(1,:),mouse_xy{vid}(2,:),'b');
hold on
plot(cricket_xy{vid}(1,:),cricket_xy{vid}(2,:),'g');

figure
bar([mean(isnan(mouse_xy{vid}(1,:))) mean(isnan(cricket_xy{vid}(1,:))) ]);
ylim([0 1]); ylabel('fraction bad')
set(gca,'Xtick',[1 2]); set(gca,'Xticklabel',{'mouse','cricket'});

figure
subplot(2,1,1);
plot(mouse_xy{vid}(1,:)); hold on; plot(mouse_xy{vid}(2,:)); title('mouse')

subplot(2,1,2);
plot(cricket_xy{vid}(1,:)); hold on; plot(cricket_xy{vid}(2,:)); title('crick')

deltaR = diff(range{vid})*30;
figure
plot(deltaR);
vsmooth = conv(mouseV{vid},ones(5,1)/5,'same');

figure
subplot(4,1,1);
plot(range{vid}); title('range');
subplot(4,1,2);
plot(deltaR); title('deltaR');
subplot(4,1,3);
plot(azT{vid}); title('az');
subplot(4,1,4);
plot(vsmooth); title('speed')

%%% identify approach!!!

dRThresh=-10; %%%cm/sec
vThresh=10;
azThresh = pi/4;  %%% pi/4 = 45 deg
approach = deltaR<dRThresh & vsmooth(1:end-1)>vThresh & abs(azT{vid}(1:end-1))<azThresh;
approach(1)=0; approach(end)=0; %%% boundary conditions

starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% find start;stop

for j= 1:length(ends)-1;  %%% stitch together approachs with gap of less than 5 frames
    if (starts(j+1)-ends(j))<5 & (range{vid}(starts(j+1))- range{vid}(ends(j)))<3
        approach(ends(j) : starts(j+1))=1;
    end
end
starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% update start;stop

for j = 1:length(starts);   %%% remove short approaches (less than 10 frames)
    if ends(j)-starts(j)<10
        approach(starts(j):ends(j))=0;
    end
end
starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% update start;stop

figure
plot(approach), ylim([0 1.5])

%%% plot with eye data!!!
Data(vid).Rtheta = Data(vid).Rtheta - nanmedian(Data(vid).Rtheta);
Data(vid).Ltheta = Data(vid).Ltheta - nanmedian(Data(vid).Ltheta);
Data(vid).Rphi = Data(vid).Rphi - nanmedian(Data(vid).Rphi);
Data(vid).Lphi = Data(vid).Lphi - nanmedian(Data(vid).Lphi);

angles = 0:0.01:2*pi; R= 5; %%% substitute actual eye radii

startLag = 10; %%% offset from beginning (neededfor trails behind mouse)

figure
set(gcf,'Position',[100 0 560 840] )
for i = (1+startLag):length(Data(vid).Rtheta)
    
    %%% plot left eye position (with circle)
    subplot(5,2,3);
    plot(Data(vid).Ltheta,Data(vid).Lphi); hold on
    plot(Data(vid).Ltheta(i)+R*cos(angles),Data(vid).Lphi(i)+R*sin(angles),'Linewidth',4);
    axis equal; axis([-60 60 -60 60]); hold off; title('left eye'); xlabel('theta'); ylabel('phi')
    
    %%% plot right eye position (with circle)
    subplot(5,2,4);
    plot(Data(vid).Rtheta,Data(vid).Rphi); hold on
    plot(Data(vid).Rtheta(i)+R*cos(angles),Data(vid).Rphi(i)+R*sin(angles),'Linewidth',4);
    axis equal; axis([-60 60 -60 60]); hold off; title('right eye'); xlabel('theta'); ylabel('phi')
    
    if approach(i)
        marker = 'go';
    else
        marker = 'ro';
    end
    
    subplot(5,2,5:6);
    plot(deltaR); hold on; plot(i,deltaR(i),marker,'LineWidth',2); xlim([0 length(deltaR)]); ylim([-50 50])
    plot(find(approach),zeros(sum(approach),1),'g.') ; plot([1 length(deltaR)],[dRThresh,dRThresh],'r:')
    title('deltaR'); hold off
    
    subplot(5,2,7:8);
    plot(vsmooth); hold on; plot(i,vsmooth(i),marker,'LineWidth',2); xlim([0 length(deltaR)]); ylim([0 50])
    plot(find(approach),zeros(sum(approach),1),'g.'); plot([1 length(deltaR)],[vThresh,vThresh],'r:')
    title('speed');hold off
    
    subplot(5,2,9:10);
    plot(azT{vid},'.'); hold on; plot(i,azT{vid}(i),marker,'LineWidth',2); xlim([0 length(deltaR)]); ylim([-pi pi])
    plot(find(approach),zeros(sum(approach),1),'g.'); plot([1 length(deltaR)],[azThresh azThresh],'r:'); plot([1 length(deltaR)],[-azThresh -azThresh],'r:');
    title('azimuth');hold off
    
    %%% plot cricket and mouse tracks
    subplot(5,2,[1:2]);
    plot(Data(vid).cricketxy(1,i),Data(vid).cricketxy(2,i),'g*');
    hold on
    plot(Data(vid).mouse_xy(1,i + (-startLag:0)), Data(vid).mouse_xy(2,i + (-10:0)),'b', 'Linewidth',2)
    if approach(i)
        plot(Data(vid).mouse_xy(1,i),Data(vid).mouse_xy(2,i),'go','LineWidth',2)
    end
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
