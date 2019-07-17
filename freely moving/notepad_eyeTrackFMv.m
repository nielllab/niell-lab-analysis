%%

figure;
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';
vidfile = VideoWriter('EyeRadboth_approachcapture_040619.mp4','MPEG-4');
open(vidfile);
long = EllipseParams(:,3); short = EllipseParams(:,4);
rad = long+short./2*length(short)
%for  v=1:size(EyeVid,3)
%  subplot(4,2,[5 6 7 8])
%  plot(rad,'-k');hold on;
for v=5800:6800%size(EyeVid,3)
    subplot(2,4,[1 2])
    imagesc(EyeVidL(:,:,v)); colormap gray; hold on; axis equal off;
    scatter(PointsxL(v,:),PointsyL(v,:),100,'.b'); hold off; % Uncomment if to plot DLC points too
%     if ~any((Pointsx(v,:)<10 | Pointsy(v,:)<10)) % check points
%         e_t = fit_ellipse(Pointsx(v,:),Pointsy(v,:),ax);
%     end

    subplot(2,4,[3 4])
    imagesc(EyeVid(:,:,v)); colormap gray; hold on; axis equal off;
    scatter(Pointsx(v,:),Pointsy(v,:),100,'.r'); hold off; % Uncomment if to plot DLC points too
%     if ~any((Pointsx(v,:)<10 | Pointsy(v,:)<10)) % check points
%         e_t = fit_ellipse(Pointsx(v,:),Pointsy(v,:),ax);
%     end
    title(sprintf('Time(s) = %d',round(v/30)));
    drawnow limitrate;
    subplot(2,4,[5 6 7 8])
    plot(normalize(medfilt1(radL,10)),'b'); xlim([5500 7000]);
    hold on
    plot(normalize(medfilt1(rad,10)),'r'); 
    
    plot(v,.9,'g*','MarkerSize', 15); ylim([.35 1]);
    plot([5940 5940],[0 1],'-k','Linewidth',4)    
    plot([6450 6450],[0 1],'-k','Linewidth',4)    

    %set(gca,'xtick',([1:1800:size(EyeVid,3)]))% tick every 30 seconds
    set(gca,'xtick',([5800:450:6800]))% tick every 30 seconds
    set(gca,'xticklabels',({'15','30','45'}),'FontSize',10)

    %set(gca,'xticklabels',({'1','2','3','4','5','6','7','8','9','10','11','12','13'}),'FontSize',10)
    
    
    xlabel('Time (s)'); ylabel ('normalized pupil radius');
    legend('left eye', 'right eye');
    drawnow limitrate;
    hold off
    F(v) = getframe(gcf); 
    writeVideo(vidfile,F(v));
    fprintf('Time = %d \n',v);
end
close(vidfile)

%3:18:3:35 approach and capture = frames 5940:6450 (run from 5800:6800)

%%



figure;
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';
vidfile = VideoWriter('EyeMvboth_approachcapture_040619.mp4','MPEG-4');
open(vidfile);
%xR = medfilt1(EllipseParams(:,1),10); yR = medfilt1(EllipseParams(:,2),10);
%xL = medfilt1(EllipseParamsL(:,1),10); yL = medfilt1(EllipseParamsL(:,2),10);



xR = (EllipseParams(:,1)); yR = (EllipseParams(:,2));
xL = (EllipseParamsL(:,1)); yL = (EllipseParamsL(:,2));

%app=(5940:6450) 
% figure
% subplot(2,2,1)
% [rxy,lags] = xcorr(xR(app),yR(app));
% plot(lags,rxy,'b'); axis square;% hold on
% title('Right eye xy corr')
% clear lags
% subplot(2,2,2)
% [lxy,lags] = xcorr(xL(app),yL(app));
% plot(lags,lxy,'r');axis square;
% title('left eye xy corr')
% 
% subplot(2,2,3)
% [lrx,lags] = xcorr(xR(app),xL(app));
% plot(lags,rxy,'k');axis square; %hold on %x both eyes
% title('x position corr both eyes')
% clear lags
% subplot(2,2,4)
% [lry,lags] = xcorr(yR(app),yL(app));
% plot(lags,lxy,'g');axis square; %y both eyes
% title('y position corr both eyes')

% xyc = corrcoef(xR(app),xL(app))




%for  v=1:size(EyeVid,3)
%  subplot(4,2,[5 6 7 8])
%  plot(rad,'-k');hold on;
for v=5800:6800%size(EyeVid,3)
    subplot(2,4,[1 2])
    imagesc(EyeVidL(:,:,v)); colormap gray; hold on; axis equal off;
    scatter(PointsxL(v,:),PointsyL(v,:),100,'.b'); hold off; % Uncomment if to plot DLC points too
%     if ~any((Pointsx(v,:)<10 | Pointsy(v,:)<10)) % check points
%         e_t = fit_ellipse(Pointsx(v,:),Pointsy(v,:),ax);
%     end

    subplot(2,4,[3 4])
    imagesc(EyeVid(:,:,v)); colormap gray; hold on; axis equal off;
    scatter(Pointsx(v,:),Pointsy(v,:),100,'.r'); hold off; % Uncomment if to plot DLC points too
%     if ~any((Pointsx(v,:)<10 | Pointsy(v,:)<10)) % check points
%         e_t = fit_ellipse(Pointsx(v,:),Pointsy(v,:),ax);
%     end
    title(sprintf('Time(s) = %d',round(v/30)));
    drawnow limitrate;
    subplot(2,4,[5 6 7 8])
     plot(xR(v),yR(v),'r.'); axis square;ylim([125 300]); xlim([300 475]);
     hold on
     plot(xL(v),yL(v),'b.')
     if v==5940
         title('APPROACH')
     end

if v>=6450
    title('CAPTURE')
end

    xlabel('x centroid'); ylabel ('y centroid');
 %   legend({'right eye', 'left eye'});
    drawnow limitrate;
    
    F(v) = getframe(gcf); 
    writeVideo(vidfile,F(v));
    fprintf('Frame = %d \n',v);
end
close(vidfile)

%%

figure;
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';
vidfile = VideoWriter('full_approachcapture_040719c.mp4','MPEG-4');
open(vidfile);

clear xR yR xL yL
xR = medfilt1(EllipseParams(:,1),10); yR = medfilt1(EllipseParams(:,2),10);
xL = medfilt1(EllipseParamsL(:,1),10); yL = medfilt1(EllipseParamsL(:,2),10);


xR = xR-mean(xR); yR = yR-mean(yR);
xL = xL-mean(xL); yL = yL-mean(yL);

%10410:10680

for v=10300:10850%5800:6800%size(EyeVid,3)
    subplot(2,3,1)
    imagesc(EyeVidL(:,:,v)); colormap gray; hold on; axis equal off;
   s= scatter(PointsxL(v,:),PointsyL(v,:),20,'.b'); hold off; % Uncomment if to plot DLC points too
 

%     if ~any((Pointsx(v,:)<10 | Pointsy(v,:)<10)) % check points
%         e_t = fit_ellipse(Pointsx(v,:),Pointsy(v,:),ax);
%     end

    subplot(2,3,2)
    imagesc(EyeVid(:,:,v)); colormap gray; hold on; axis equal off;
   s= scatter(Pointsx(v,:),Pointsy(v,:),20,'.r'); hold off; % Uncomment if to plot DLC points too
%     if ~any((Pointsx(v,:)<10 | Pointsy(v,:)<10)) % check points
%         e_t = fit_ellipse(Pointsx(v,:),Pointsy(v,:),ax);
%     end
    title(sprintf('Time(s) = %d',round(v/30)));
    drawnow limitrate;
    subplot(2,3,3)
     plot(xR(v),yR(v),'r.'); axis square;%ylim([-60 60]); xlim([-100 100]);
     hold on
     plot(xL(v),yL(v),'b.')
     if v==10410%5940
         title('APPROACH')
     end

if v>=10680 %6450 
    title('CAPTURE')
end

    xlabel('x centroid'); ylabel ('y centroid');
 %   legend({'right eye', 'left eye'});
    drawnow limitrate;
    
    
    
     subplot(2,3,[4 5 6])
    plot(normalize(medfilt1(radL,10)),'b'); xlim([10000 11000]);
    hold on
    plot(normalize(medfilt1(rad,10)),'r'); 
    
    plot(v,.9,'g*','MarkerSize', 15); ylim([.35 1]);
%     plot([5940 5940],[0 1],'-k','Linewidth',4)    
%     plot([6450 6450],[0 1],'-k','Linewidth',4)    

   plot([10390 10390],[0 1],'-k','Linewidth',4)    
    plot([10720 10720],[0 1],'-k','Linewidth',4)  

    %set(gca,'xtick',([1:1800:size(EyeVid,3)]))% tick every 30 seconds
    %set(gca,'xtick',([5800:450:6800]))% tick every 30 seconds
    set(gca,'xtick',([10250:450:10800]))% tick every 30 seconds

    set(gca,'xticklabels',({'15','30','45'}),'FontSize',10)

    %set(gca,'xticklabels',({'1','2','3','4','5','6','7','8','9','10','11','12','13'}),'FontSize',10)
    
    
    xlabel('Time (s)'); ylabel ('normalized pupil radius');
 %   legend('left eye', 'right eye');
    drawnow limitrate;
    hold off
    
    
    F(v) = getframe(gcf); 
    writeVideo(vidfile,F(v));
    fprintf('Frame = %d \n',v);
end
close(vidfile)



%%
clear xR yR xL yL
xR = medfilt1(EllipseParams(:,1),10); yR = medfilt1(EllipseParams(:,2),10);
xL = medfilt1(EllipseParamsL(:,1),10); yL = medfilt1(EllipseParamsL(:,2),10);


xR = xR-mean(xR); yR = yR-mean(yR);
xL = xL-mean(xL); yL = yL-mean(yL);

app = 5940:6450
app2 =10410:10680 %5:47-5:56
app3= 18270:18720 %10:09-10:24

figure
subplot(2,2,1)
plot(xR,yR,'.'); hold on;
plot(xR(app),yR(app),'.r')
plot(xR(app2),yR(app2),'.g')
plot(xR(app3),yR(app3),'.y')
xlim([-200 200]); ylim([-100 100]);
title('xy position, R eyes')
xlabel('x position'); ylabel('yposition')


subplot(2,2,2)
plot(xL,yL,'.'); hold on
plot(xL(app),yL(app),'r.');
plot(xL(app2),yL(app2),'g.');
plot(xL(app3),yL(app3),'y.');
xlim([-200 200]); ylim([-150 150]);
xlabel('x position'); ylabel('yposition')
title('xy position, L eye')

subplot(2,2,3)
plot(xR,xL(1:end-1),'.');hold on ;%xlim([200 600]); ylim([200 600]);
plot(xR(app),xL(app),'r.');
plot(xR(app2),xL(app2),'g.');
plot(xR(app3),xL(app3),'y.');
xlim([-100 100]);% ylim([-100 100]);
xlabel('x position'); ylabel('yposition')
title('x position, both eyes')

subplot(2,2,4)
plot(yR,yL(1:end-1),'.'); hold on
plot(yR(app),yL(app),'r.')
plot(yR(app2),yL(app2),'g.')
plot(yR(app3),yL(app3),'y.')
xlim([-100 100]);% ylim([-100 100]);
xlabel('x position'); ylabel('yposition')

title('y position, both eyes')

%%

%Left Eye
[thetaL,phiL,EllipseParamsL,~] = EyeCameraCalc1(EyeVidL,PointsxL,PointsyL);

%Right Eye
[thetaR,phiR,EllipseParamsR,~] = EyeCameraCalc1(EyeVid,Pointsx,Pointsy);


%%






figure;plot(radR(1:20:length(radR)),radL(1:20:length(radL)),'ko','Markersize', 4);axis square; xlabel('Right Eye Normalized Pupil Radius');
ylabel('Left Eye Normalized Pupil Radius');


% xR = xR-mean(xR); yR = yR-mean(yR);
% xL = xL-mean(xL); yL = yL-mean(yL);
% 
app = 5940:6450
app2 =10410:10680 %5:47-5:56
app3= 18270:18720 %10:09-10:24

figure;plot(xR, yR,'.r');hold on; plot(xL,yL,'.b');
plot(xR(app), yR(app),'.g');hold on; plot(xL(app),yL(app),'.g');
plot(xR(app2), yR(app2),'.y');hold on; plot(xL(app2),yL(app2),'.y');
plot(xR(app3), yR(app3),'.c');hold on; plot(xL(app3),yL(app3),'.c');

velThetaL=diff(thetaL); velThetaR=diff(thetaR);
%%

figure;
subplot(1,2,1)
plot(velThetaR,velThetaL,'.'); hold on
plot(velThetaR(app), velThetaL(app),'.g'); axis equal; axis square
plot(velThetaR(app2), velThetaL(app2),'.y');
plot(velThetaR(app3), velThetaL(app3),'.m');
xlabel('Theta R Eye');ylabel('Theta L Eye');

velPhiL=diff(phiL); velPhiR=diff(phiR);

subplot(1,2,2)
plot(velPhiR,velPhiL(1:end-1),'.'); hold on
plot(velPhiR(app), velPhiL(app),'.g'); axis square
plot(velPhiR(app2), velPhiL(app2),'.y');
plot(velPhiR(app3), velPhiL(app3),'.m');
xlabel('Phi R Eye');ylabel('Phi L Eye');

%%
clear cricketX
cricketX=PointsxT(:,7);
cln =LHT(:,7)>.90


cricketY=PointsyT(:,7)
clny =LHT(:,7)>.90
figure;plot(cricketX(cln),cricketY(clny)); hold on;
plot(PointsxT(:,6), PointsyT(:,6));

%%
cricketX(cln==0) =NaN;

distX=medfilt1(cricketX-PointsxT(:,6),10); distY=medfilt1(cricketY-PointsyT(:,6),10);
figure;plot(normalize(distX)); hold on; plot(radR);

%%

ap1 = 30:40:270
ap2 = 2700:40:3690
ap3 = 6690:40:7230
ap4 = 10470:40:10550
ap5 = 13440:40:13560 %7:28 -7:32
ap6 = 14220:40:14430% 7:54 - 8:01
ap7 = 14820:40:14910  %8:14-8:17
ap8 = 15210:40:15600;      %8:27-8:40

% ap9 = 15240:15540 
mnRadR = (EllipseParamsR(:,3))./(mean(EllipseParamsR(:,3)));
mnRadL = (EllipseParamsL(:,3))./(mean(EllipseParamsL(:,3)));


radR = medfilt1(mnRadR,20); radL = medfilt1(mnRadL,20);
figure;
plot(radR,'r'); hold on; plot(radL,'b')
sec = length(radR)./30; min = sec./60
set(gca,'xtick',([0:1800:length(radR)]))% tick every 1 min
set(gca,'xticklabels',({'0','1','2','3','4','5','6','7','8','9','10'}),'FontSize',10)
xlabel('Time(min)'); ylabel('normalized pupil radius')

plot(ap1,radR(ap1),'g'); plot(ap2,radR(ap2),'g');plot(ap3,radR(ap3),'g');
plot(ap4,radR(ap4),'g'); plot(ap5,radR(ap5),'g'); plot(ap6,radR(ap6),'g');plot(ap7,radR(ap7),'g');
plot(ap8,radR(ap8),'g');

plot(ap1,radL(ap1),'g'); plot(ap2,radL(ap2),'g');plot(ap3,radL(ap3),'g');
plot(ap4,radL(ap4),'g'); plot(ap5,radL(ap5),'g'); plot(ap6,radL(ap6),'g');plot(ap7,radL(ap7),'g');
plot(ap8,radL(ap8),'g');
%%


figure;plot(radR(1:40:length(radR)),radL(1:40:length(radL)),'ko','Markersize', 4);axis square; xlabel('Right Eye Normalized Pupil Radius');
ylabel('Left Eye Normalized Pupil Radius');
hold on; plot(radR(ap1),radL(ap1),'og'); plot(radR(ap2),radL(ap2),'go');
plot(radR(ap3),radL(ap3),'og'); plot(radR(ap4),radL(ap4),'go');
plot(radR(ap5),radL(ap5),'og'); plot(radR(ap6),radL(ap6),'go');
plot(radR(ap7),radL(ap7),'og'); plot(radR(ap8),radL(ap8),'go');

%%
clear xR yR xL yL
xR = medfilt1(EllipseParamsR(:,1),10); yR = medfilt1(EllipseParamsR(:,2),10);
xL = medfilt1(EllipseParamsL(:,1),10); yL = medfilt1(EllipseParamsL(:,2),10);

figure;plot(xR, yR,'.r');hold on; plot(xL,yL,'.b');
plot(xR(ap1), yR(ap1),'.g');hold on; plot(xL(ap1),yL(ap1),'.g');
plot(xR(ap2), yR(ap2),'.g');hold on; plot(xL(ap2),yL(ap2),'.g');
plot(xR(ap3), yR(ap3),'.g');hold on; plot(xL(ap3),yL(ap3),'.g');

%%
ap1 = 30:30:270
ap2 = 2700:30:3690
ap3 = 6690:30:7230
ap4 = 10470:30:10550
ap5 = 13440:30:13560 %7:28 -7:32
ap6 = 14220:30:14430% 7:54 - 8:01
ap7 = 14820:30:14910  %8:14-8:17
ap8 = 15210:30:15600;      %8:27-8:40

%thetaR = medfilt1(thetaR,30);thetaL = medfilt1(thetaL,30)
%phiR = medfilt1(phiR,30);phiL = medfilt1(phiL,30)

thetaRfilt=thetaR(1:30:length(thetaR));thetaLfilt=thetaL(1:30:length(thetaL))
phiRfilt=phiR(1:30:length(phiR));phiLfilt=phiL(1:30:length(phiL))

figure;plot(thetaRfilt, thetaLfilt,'ok'); hold on; axis square; %xlim([-80 80]); ylim([-80 80]);
plot(thetaR(ap1),thetaL(ap1),'go'); plot(thetaR(ap2),thetaL(ap2),'go');
plot(thetaR(ap3),thetaL(ap3),'go'); plot(thetaR(ap4),thetaL(ap4),'go');
plot(thetaR(ap5),thetaL(ap5),'go'); plot(thetaR(ap6),thetaL(ap6),'go');
plot(thetaR(ap7),thetaL(ap7),'go'); plot(thetaR(ap8),thetaL(ap8),'go');
xlabel ('R theta'); ylabel('L theta');


figure;plot(phiRfilt, phiLfilt,'ok'); hold on; axis square; %xlim([-80 80]); ylim([-80 80]);
plot(phiR(ap1),phiL(ap1),'go'); plot(phiR(ap2),phiL(ap2),'go');
plot(phiR(ap3),phiL(ap3),'go'); plot(phiR(ap4),phiL(ap4),'go');
plot(phiR(ap5),phiL(ap5),'go'); plot(phiR(ap6),phiL(ap6),'go');
plot(phiR(ap7),phiL(ap7),'go'); plot(phiR(ap8),phiL(ap8),'go');
xlabel ('R phi'); ylabel('L phi'); xlim([-80 80]);ylim([-80 80]);


%%


figure;plot(thetaR); hold on; plot(thetaL)
plot(ap1,thetaR(ap1),'g');hold on; plot(ap1,thetaL(ap1),'g');
plot(ap2,thetaR(ap2),'g');hold on; plot(ap2,thetaL(ap2),'g');
plot(ap3,thetaR(ap3),'g');hold on; plot(ap3,thetaL(ap3),'g');
plot(ap4,thetaR(ap4),'g');hold on; plot(ap4,thetaL(ap4),'g');
plot(ap5,thetaR(ap5),'g');hold on; plot(ap5,thetaL(ap5),'g');
plot(ap6,thetaR(ap6),'g');hold on; plot(ap6,thetaL(ap6),'g');
plot(ap7,thetaR(ap7),'g');hold on; plot(ap7,thetaL(ap7),'g');
plot(ap8,thetaR(ap8),'g');hold on; plot(ap8,thetaL(ap8),'g');


%%

ap1 = 30:270
ap2 = 2700:3690
ap3 = 6690:7230
ap4 = 10470:10550
ap5 = 13440:13560 %7:28 -7:32
ap6 = 14220:14430% 7:54 - 8:01
ap7 = 14820:14910  %8:14-8:17
ap8 = 15210:15600;    
velThetaL=diff(thetaL); velThetaR=diff(thetaR);
velThetaL= velThetaL(1:30:length(velThetaL)); velThetaR = velThetaR(1:30:length(velThetaR))

figure;
subplot(1,2,1)
plot(velThetaR,velThetaL,'o'); hold on
plot(velThetaR(ap1), velThetaL(ap1),'og'); axis equal; axis square
plot(velThetaR(ap2), velThetaL(ap2),'og');
plot(velThetaR(ap3), velThetaL(ap3),'og');
plot(velThetaR(ap4), velThetaL(ap4),'og');
plot(velThetaR(ap5), velThetaL(ap5),'og');
plot(velThetaR(ap6), velThetaL(ap6),'og');
plot(velThetaR(ap7), velThetaL(ap7),'og');
plot(velThetaR(ap8), velThetaL(ap8),'og');
ylim([-9 9]); xlim([-15 15]);

xlabel('Theta R Eye');ylabel('Theta L Eye');

velPhiL=diff(phiL); velPhiR=diff(phiR);

subplot(1,2,2)
plot(velPhiR,velPhiL,'o'); hold on
plot(velPhiR(ap1), velPhiL(ap1),'og'); axis square
plot(velPhiR(ap2), velPhiL(ap2),'og');
plot(velPhiR(ap3), velPhiL(ap3),'og');
plot(velPhiR(ap4), velPhiL(ap4),'og');
plot(velPhiR(ap5), velPhiL(ap5),'og');
plot(velPhiR(ap6), velPhiL(ap6),'og');
plot(velPhiR(ap7), velPhiL(ap7),'og');
plot(velPhiR(ap8), velPhiL(ap8),'og');
ylim([-11 11]); xlim([-25 25]);


xlabel('Phi R Eye');ylabel('Phi L Eye');


%%
clr = jet(length(cricketX));
figure
plot(cricketX(cln),cricketY(clny));hold on
%scatter(cricketX(cln),cricketY(clny),10,clr);colormap jet


%%
ear(:,1,1) = PointsxT(:,2);  %%% left
ear(:,2,1) = PointsyT(:,2);
ear(:,3,1) = LHT(:,2);

ear(:,1,2) = PointsxT(:,4);  %%% left
ear(:,2,2) = PointsyT(:,4);
ear(:,3,2) = LHT(:,4);


%angle =atan2(left ear x - right ear x, left ear y-right ear Y) 

headangle= radtodeg(atan2(PointsxT(:,2)-PointsxT(:,4),PointsyT(:,2)-PointsyT(:,4)))
HA= medfilt1(diff(headangle),10)

figure;subplot(2,1,1);plot(headangle);hold on; ylabel('head angle (degrees)')
plot(ap1,(headangle(ap1)),'g')
plot(ap2,(headangle(ap2)),'g')
plot(ap3,(headangle(ap3)),'g')
plot(ap4,(headangle(ap4)),'g')
plot(ap5,(headangle(ap5)),'g')
plot(ap6,(headangle(ap6)),'g')
plot(ap7,(headangle(ap7)),'g')
plot(ap8,(headangle(ap8)),'g')


figure;plot(HA); title('velocity'); hold on
plot(ap1,(HA(ap1)),'g')
plot(ap2,(HA(ap2)),'g')
plot(ap3,(HA(ap3)),'g')
plot(ap4,(HA(ap4)),'g')
plot(ap5,(HA(ap5)),'g')
plot(ap6,(HA(ap6)),'g')
plot(ap7,(HA(ap7)),'g')
plot(ap8,(HA(ap8)),'g')

set(gca,'xtick',([0:1800:length(HA)]))% tick every 1 min
set(gca,'xticklabels',({'0','1','2','3','4','5','6','7','8','9','10'}),'FontSize',10)
xlabel('Time (min)'); ylabel('head angle velocity')

%%

figure;histogram(HA);hold on;
histogram(HA(ap1),'Facecolor','g')
histogram(HA(ap2),'Facecolor','g')
histogram(HA(ap3),'Facecolor','g')
histogram(HA(ap4),'Facecolor','g')
histogram(HA(ap5),'Facecolor','g')
histogram(HA(ap6),'Facecolor','g')
histogram(HA(ap7),'Facecolor','g')
histogram(HA(ap8),'Facecolor','g')

%%

nHA = normalize(headangle)
figure;plot(nHA); hold on;
plot(normalize(thetaR))


figure;plot(headanglefilt(1:end-1),thetaRfilt,'.');% hold on; plot(thetaR)

%%
clear xR yR xL yL
xR = medfilt1(EllipseParamsR(:,1),10); yR = medfilt1(EllipseParamsR(:,2),10);
xL = medfilt1(EllipseParamsL(:,1),10); yL = medfilt1(EllipseParamsL(:,2),10);


xR = xR-mean(xR); yR = yR-mean(yR);
xL = xL-mean(xL); yL = yL-mean(yL);

ap1 = 30:270
ap2 = 2700:3690
ap3 = 6690:7230
ap4 = 10470:10550
ap5 = 13440:13560 %7:28 -7:32
ap6 = 14220:14430% 7:54 - 8:01
ap7 = 14820:14910  %8:14-8:17
ap8 = 15210:15600

figure;plot(xR,xL,'o'); hold on; xlim([-50 50]);ylim([-50 50])
plot(xR(ap1),xL(ap1),'go')
plot(xR(ap2),xL(ap2),'go')
plot(xR(ap3),xL(ap3),'go')
plot(xR(ap4),xL(ap4),'go')
plot(xR(ap5),xL(ap5),'go')
plot(xR(ap6),xL(ap6),'go')
plot(xR(ap7),xL(ap7),'go')
plot(xR(ap8),xL(ap8),'go')
xlabel('x R eye');ylabel('x Left eye')

figure;plot(xR,yR,'o'); hold on; xlim([-50 50]);ylim([-50 50])
plot(xR(ap1),yR(ap1),'go')
plot(xR(ap2),yR(ap2),'go')
plot(xR(ap3),yR(ap3),'go')
plot(xR(ap4),yR(ap4),'go')
plot(xR(ap5),yR(ap5),'go')
plot(xR(ap6),yR(ap6),'go')
plot(xR(ap7),yR(ap7),'go')
plot(xR(ap8),yR(ap8),'go')
xlabel('x R eye');ylabel('Y right eye')

figure;plot(xL,yL,'o'); hold on; xlim([-50 50]);ylim([-50 50])
plot(xL(ap1),yL(ap1),'go')
plot(xL(ap2),yL(ap2),'go')
plot(xL(ap3),yL(ap3),'go')
plot(xL(ap4),yL(ap4),'go')
plot(xL(ap5),yL(ap5),'go')
plot(xL(ap6),yL(ap6),'go')
plot(xL(ap7),yL(ap7),'go')
plot(xL(ap8),yL(ap8),'go')
xlabel('x L eye');ylabel('Y L eye')


figure;plot(yR,yL,'o'); hold on; ylim([-50 50]);ylim([-50 50])
plot(yR(ap1),yL(ap1),'go')
plot(yR(ap2),yL(ap2),'go')
plot(yR(ap3),yL(ap3),'go')
plot(yR(ap4),yL(ap4),'go')
plot(yR(ap5),yL(ap5),'go')
plot(yR(ap6),yL(ap6),'go')
plot(yR(ap7),yL(ap7),'go')
plot(yR(ap8),yL(ap8),'go')
xlabel('y R eye');ylabel('Y L eye')


%%
ap1 = 30:270
ap2 = 2700:3690
ap3 = 6690:7230
ap4 = 10470:10550
ap5 = 13440:13560 %7:28 -7:32
ap6 = 14220:14430% 7:54 - 8:01
ap7 = 14820:14910  %8:14-8:17
ap8 = 15210:15600;      %8:27-

xT= PointsxT(:,6); yT= PointsyT(:,6)
xC= PointsxT(:,7); yC= PointsyT(:,7)


figure;plot(xT,yT);
hold on;
plot(xT(ap1),yT(ap1),'g')
plot(xT(ap2),yT(ap2),'g')
plot(xT(ap3),yT(ap3),'g')
plot(xT(ap4),yT(ap4),'g')
plot(xT(ap5),yT(ap5),'g')
plot(xT(ap6),yT(ap6),'g')
plot(xT(ap7),yT(ap7),'g')
plot(xT(ap8),yT(ap8),'g')

plot(xC(ap1(1)),yC(ap1(1)),'or')
plot(xC(ap2(1)),yC(ap2(1)),'or')
plot(xC(ap3(1)),yC(ap3(1)),'or')
plot(xC(ap4(1)),yC(ap4(1)),'or')
plot(xC(ap5(1)),yC(ap5(1)),'or')
plot(xC(ap6(1)),yC(ap6(1)),'or')
plot(xC(ap7(1)),yC(ap7(1)),'or')
plot(xC(ap8(1)),yC(ap8(1)),'or')



%%
figure;plot(thetaRfilt, thetaLfilt,'ok'); hold on; axis square; %xlim([-80 80]); ylim([-80 80]);
plot(thetaR(app1),thetaL(app1),'go'); plot(thetaR(app2),thetaL(app2),'go');
plot(thetaR(app3),thetaL(app3),'go'); plot(thetaR(app4),thetaL(app4),'go');
plot(thetaR(app5),thetaL(app5),'go'); plot(thetaR(app6),thetaL(app6),'go');
plot(thetaR(app7),thetaL(app7),'go'); plot(thetaR(app8),thetaL(app8),'go');
plot(thetaR(app9),thetaL(app9),'go');
xlabel ('R theta'); ylabel('L theta');


figure;plot(phiRfilt, phiLfilt,'ok'); hold on; axis square; %xlim([-80 80]); ylim([-80 80]);
plot(phiR(app1),phiL(app1),'go'); plot(phiR(app2),phiL(app2),'go');
plot(phiR(app3),phiL(app3),'go'); plot(phiR(app4),phiL(app4),'go');
plot(phiR(app5),phiL(app5),'go'); plot(phiR(app6),phiL(app6),'go');
plot(phiR(app7),phiL(app7),'go'); plot(phiR(app8),phiL(app8),'go');
plot(phiR(app9),phiL(app9),'go');
xlabel ('R phi'); ylabel('L phi'); xlim([-80 80]);ylim([-80 80]);


