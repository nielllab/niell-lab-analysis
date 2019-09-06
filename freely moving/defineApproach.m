close all
framerate =30;
for m = 1:length(useData);

    vid = useData(m);


deltaR = diff(dist{vid})*30;
figure
plot(deltaR);
vsmooth = conv(mouseSp{vid},ones(5,1)/5,'same');

figure
subplot(4,1,1);
plot(dist{vid}); title('range');
subplot(4,1,2);
plot(deltaR); title('deltaR');
subplot(4,1,3);
plot(az{vid}); title('az');
subplot(4,1,4);
plot(vsmooth); title('speed')

% %%% identify approach!!!
% 
% dRThresh=-10; %%%cm/sec
% vThresh=10;
% azThresh = pi/4;  %%% pi/4 = 45 deg
% approach = deltaR<dRThresh & vsmooth(1:end-1)>vThresh & abs(az{vid}(1:end-1))<azThresh;
% approach(1)=0; approach(end)=0; %%% boundary conditions
% 
% starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% find start;stop
% 
% for j= 1:length(ends)-1;  %%% stitch together approachs with gap of less than 5 frames
%     if (starts(j+1)-ends(j))<5 & (range{vid}(starts(j+1))- range{vid}(ends(j)))<3
%         approach(ends(j) : starts(j+1))=1;
%     end
% end
% starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% update start;stop
% 
% for j = 1:length(starts);   %%% remove short approaches (less than 10 frames)
%     if ends(j)-starts(j)<10
%         approach(starts(j):ends(j))=0;
%     end
% end
% starts = find(diff(approach)>0);  ends = find(diff(approach)<0);  %%% update start;stop

approach = appEpoch{m};
figure
plot(approach), ylim([0 1.5])



Rth = Rtheta{vid} - nanmedian(Rtheta{vid});
Lth = Ltheta{vid} - nanmedian(Ltheta{vid});
Rph = Rphi{vid} - nanmedian(Rphi{vid});
Lph = Lphi{vid} - nanmedian(Lphi{vid});

angles = 0:0.01:2*pi; R= 5; %%% substitute actual eye radii

startLag = 10; %%% offset from beginning (neededfor trails behind mouse)

figure
set(gcf,'Position',[100 0 560 840] )
for i = (1+startLag):length(Rth)-1
    
    %%% plot left eye position (with circle)
    subplot(5,2,3);
    plot(Lth,Lph); hold on
    plot(Lth(i)+R*cos(angles),Lph(i)+R*sin(angles),'Linewidth',2);
    axis equal; axis([-45 45 -45 45]); hold off; title('left eye'); xlabel('theta'); ylabel('phi')
    
    %%% plot right eye position (with circle)
    subplot(5,2,4);
    plot(Rth,Rph); hold on
    plot(Rth(i)+R*cos(angles),Rph(i)+R*sin(angles),'Linewidth',2);
    axis equal; axis([-45 45 -45 45]); hold off; title('right eye'); xlabel('theta'); ylabel('phi')
    
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
    plot(az{vid},'.'); hold on; plot(i,az{vid}(i),marker,'LineWidth',2); xlim([0 length(deltaR)]); ylim([-pi pi])
    plot(find(approach),zeros(sum(approach),1),'g.'); plot([1 length(deltaR)],[azThresh azThresh],'r:'); plot([1 length(deltaR)],[-azThresh -azThresh],'r:');
    title('azimuth');hold off
    
    %%% plot cricket and mouse tracks
    subplot(5,2,[1:2]);
    plot(cricketPos{vid}(1,i),cricketPos{vid}(2,i),'g*');
    hold on
    plot(mousePos{vid}(1,i + (-startLag:0)), mousePos{vid}(2,i + (-startLag:0)),'b', 'Linewidth',2)
    if approach(i)
        plot(mousePos{vid}(1,i),mousePos{vid}(2,i),'go','LineWidth',2)
    end
  
    %%% calculate gaze direction (head + eyes); assume each eye is centered
    %%% at 45deg (pi/4)
    
    eyeOffset = 20;
    rth = thetaHead{vid}(i) + eyeOffset + Rtheta{vid}(i);
    lth = thetaHead{vid}(i) - eyeOffset + Ltheta{vid}(i);
    %     rth = Data(vid).theta(i) + Data(vid).Rtheta(i)*pi/180;
    %     lth = Data(vid).theta(i) + Data(vid).Ltheta(i)*pi/180;
    
    plot(mousePos{vid}(1,i) + [0 400*cosd(rth) ] ,mousePos{vid}(2,i) + [0 400*sind(rth)],'c', 'Linewidth',2)
    plot(mousePos{vid}(1,i) + [0 400*cosd(lth) ] ,mousePos{vid}(2,i) + [0 400*sind(lth)],'m', 'Linewidth',2)
    
      %%% calculate head vector
    hx = 200*cosd(thetaHead{vid}(i));
    hy = 200*sind(thetaHead{vid}(i));
    plot(mousePos{vid}(1,i) + [0 hx ] ,mousePos{vid}(2,i) + [0 hy],'k', 'Linewidth',2)
    
    
    axis equal; axis([0 1600 0 1200]);  hold off
    
    drawnow
    mov(i-startLag) = getframe(gcf);
end


%% save video

f = sprintf('%s_%s_%s_%s.avi',ani{vid},date{vid},sess{vid},clipnum{vid})
movObj = VideoWriter(f);
open(movObj);
writeVideo(movObj,mov);
close(movObj);

close all

end
