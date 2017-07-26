function [h_d targ body r thetaR thetaSac vhead vbody abody vtarg atarg head sub] = analyzeCapture(fname,fps,scale,sub)
if isempty(fname)
[f p] = uigetfile('*.*','behav data');
fname = fullfile(p,f);
end
try
data = xlsread(fname);
catch
    data =dlmread(fname);
end

 missing = isnan(mean(data(:,7),2));
 data = data(~missing,:);
 data=data(4:end,:);%excludes boundry coordinates and column label

figure

% read in points
body = data(:,5:6);
pt1 = data(:,3:4);%left ear coordinates
pt2 = data(:,1:2);%right ear
targ = data(:,7:8);
head = 0.5*(pt1+pt2);

% final = min(find(isnan(targ(:,1))))-1;
% if ~isempty(final)
%     head = head(1:final,:);
% pt1 = pt1(1:final,:); pt2=pt2(1:final,:);
% body = body(1:final,:);
% targ=targ(1:final,:);
% end

% plot head and target trajectories
subplot(2,3,1)
plot(head(:,1),head(:,2));
hold on
plot(targ(:,1),targ(:,2),'g');
axis([ 0 2000 0 1200])
title (sub)

% distance between head and target
r = sqrt((head(:,1) - targ(:,1)).^2 + (head(:,2) - targ(:,2)).^2);
r=r./scale;
subplot(2,3,2)
plot((1:length(r))/fps,r); ylabel('distance to target (cm)');
xlim([1 length(r)]/fps); ylim([0 80]);hold on
% plot((1:length(isTouch))/fps,isTouch*10,'r')

%% mouse and target speed


%% mouse velocity
vbody = diff(body); vbody = sqrt(vbody(:,1).^2 + vbody(:,2).^2);%distance/60hz
abody= diff(vbody);
vbody = conv(vbody,ones(1,12),'same'); %filters/smoothes the velocity data
vscale=vbody/5 % to put on scale more comparable to acceleration
abody = conv(abody,ones(1,12),'same');

%  figure
%  plot(vbody); ylabel('body speed');hold on
%  plot(vscale,'g');hold on
%  plot(abody,'r');hold on

 %% prey speed
 
vtarg = diff(targ); vtarg = sqrt(vtarg(:,1).^2 + vtarg(:,2).^2);
atarg= diff(vtarg);
vtarg = conv(vtarg,ones(1,12),'same'); %filters the velcity data
vscaleT=vtarg/5
atarg= conv(atarg,ones(1,12),'same');
ascale=atarg/2

% figure
%  plot(ascale,'r');hold on
%  plot(abody,'g');hold on

%%acceleration correlation
% [axc lags] = xcorr(abody,atarg,'coeff');
% %vxc = conv(vxc,ones(1,100),'same')/100;
% axc = axc(lags>=-300 & lags<=300);
% lags = lags(lags>=-300 & lags<=300);
% subplot(2,3,4)
% plot(lags/fps,axc,'Linewidth',2);
% xlabel('secs')
% ylabel('acceleration body corr'); %legend('head vs targ angular vel')
% axis([ -2 2 -0.25 0.5 ]); hold on; plot([0 0],[-0.25 0.5],'r:')
 
%% get angles
headvec = pt2-pt1;
headtheta = getSmoothAngle(headvec)-pi/2;%head angle relative to coordinate system

targvec = targ-head;
targtheta = getSmoothAngle(targvec); %target angle relative to coordinate system

neckvec = head-body;
necktheta = getSmoothAngle(neckvec);

saccadetheta=necktheta-headtheta;
Rtheta=headtheta-targtheta;

%% head angle and target position
% subplot(2,3,3)
% plot((1:length(headtheta))/fps,mod(headtheta*180/pi,360)); hold on; plot((1:length(headtheta))/fps,mod(targtheta*180/pi,360),'g');
% legend('head','targ'); ylabel('theta'); xlabel('time');
% xlim([1 length(headtheta)]/fps)

thetaH=mod(headtheta*180/pi,360);
thetaT=mod(targtheta*180/pi,360);
thetaR=mod(Rtheta*180/pi,360);
thetaSac= mod(saccadetheta*180/pi,360);

%% head angle and change in distance
thetaR_sym=mod(thetaR-180,360);
thetaSac_sym=mod(thetaSac-180,360);

subplot(2,3,3);
%figure
binr=0:5:50;
bint=0:10:360;
hrt=hist2(r,thetaR_sym,binr,bint);
hrt= hrt/sum(hrt(:));
imagesc(flipud(hrt));
title 'ThetaTarget'

subplot(2,3,4);
binr=0:5:50;
bint=0:10:360;
hrt=hist2(r,thetaSac_sym,binr,bint);
hrt= hrt/sum(hrt(:));
imagesc(flipud(hrt));
title 'head saccade/thetaBody'

%% calculate correlation between head angle and target location
% opt = 'unbiased';
% [xc lags] = xcorr(headtheta(r>200),targtheta(r>200),opt); t0= find(lags==0);
% cc = corrcoef(headtheta(r>200),targtheta(r>200));
% subplot(2,3,4)
% hold on; 
% [headxc] = xcorr(headtheta(r>200),headtheta(r>200),opt);
% hold on; plot(lags,headxc/headxc(t0),'Linewidth',2);
% [targxc] = xcorr(targtheta(r>200),targtheta(r>200),opt);
% plot(lags,targxc/targxc(t0),'g','Linewidth',2)
% plot(lags,xc*cc(1,2)/xc(t0),'k','Linewidth',2); 
% plot([0 0], [0 1.2],'r');
% legend('head','targ','head vs targ')
% xlabel('lag');
% ylabel('correlation')
% headxc = headxc(lags>=-100 & lags<=100);
% targxc = targxc(lags>=-100 & lags<=100);
% xc = xc(lags>=-100 & lags<=100);

%% if I want to plot head angle correlations between prey and mouse
% lags = -300:300;
% subplot(2,3,4); hold off
% xc=lagCircRMS(headtheta,targtheta,lags);
% plot(lags/fps,xc,'k');
% hold on
% plot(lags/fps,lagCircRMS(headtheta,headtheta,lags),'b')
% plot(lags/fps,lagCircRMS(targtheta,targtheta,lags),'g')
% axis([lags(1)/fps lags(end)/fps -1 1]); plot([0 0], [-1 1],'r:');
% ylabel('angle corr')

%% correlation between target movement and head movement
dt = .1*fps;

vhead = headtheta((1+dt):end)- headtheta(1:end-dt);
vtarg = targtheta((1+dt):end)- targtheta(1:end-dt);

% subplot(2,3,5); hold off
% plot(vhead,'b');
% hold on
% plot(vtarg,'g')
% xlim([1 length(vtarg)])

% %plot(vhead); hold on; plot(vtarg,'g')
% [vxc lags] = xcorr(vhead,vtarg,'coeff');
% %vxc = conv(vxc,ones(1,100),'same')/100;
% vxc = vxc(lags>=-120 & lags<=120);
% lags = lags(lags>=-120 & lags<=120);
% subplot(2,3,5)
% plot(lags/fps,vxc,'Linewidth',2);
% xlabel('secs')
% ylabel('speed corr'); %legend('head vs targ angular vel')
% axis([ -2 2 -0.25 0.5 ]); hold on; plot([0 0],[-0.25 0.5],'r:')
% 

%%% how far off? head angle and distance
% figure
% plot(r, 'Linewidth',2); hold on; plot((headtheta-targtheta)*180/pi,':','Linewidth',2)

%% location of target relative to animal's view

if max(r)>=6
r_dis=r(r>=6 & r<50);
Rtheta_dis=Rtheta(r>=6 & r<50);
x_d = r_dis.*cos(Rtheta_dis); y_d = r_dis.*sin(Rtheta_dis);
subplot(2,3,6)
plot(y_d,x_d); hold on; plot(0,0,'r*'); axis([-50 50 -50 50]);

%figure 
bins = -80:5:80;
h_d = hist2(x_d,y_d,bins,bins);
h_d= h_d/sum(h_d(:));

else 
    h_d(:,:)=NaN;
    
% x = r.*cos(Rtheta); y = r.*sin(Rtheta);
% subplot(2,3,6)
% plot(y,x); hold on; plot(0,0,'r*'); axis([-1000 1000 -1000 1000])
% %figure 
% bins = -1025:50:1025;
% h = hist2(x_d,y_d,bins,bins);
% h= h_d/sum(h(:)); 
end
% imagesc(flipud(h))
% hold on
% plot(floor(length(bins)/2) ,floor(length(bins)/2) + 1,'ro','Linewidth',2)
%%
figure
hist(r,1:2:80);
close all



