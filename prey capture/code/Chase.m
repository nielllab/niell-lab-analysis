function [targ, head, r, thetaR_sym, vhead, vtarg ]= Chase(Left, Right, B, CK,scale,fps,sub )

%% what to do with Nan 
% missing = isnan(mean(data(:,[1,3,7]),2));
%  data = data(~missing,:);
%  data=data(4:end,:);%excludes boundry coordinates and column label

body = B;
pt1 = Left;%left ear coordinates
pt2 = Right;%right ear
targ = CK;
head = 0.5*(pt1+pt2);

% plot head and target trajectories
subplot(2,3,1)
plot(head(:,1),head(:,2));
hold on
plot(targ(:,1),targ(:,2),'b');
axis([ 0 1000 0 1200])
title (sub)

% distance between head and target
r = sqrt((head(:,1) - targ(:,1)).^2 + (head(:,2) - targ(:,2)).^2);
r=r./scale;
subplot(2,3,2)
plot((1:length(r))/fps,r); ylabel('distance to target (cm)');
xlim([1 length(r)]/fps); ylim([0 80]);hold on
% plot((1:length(isTouch))/fps,isTouch*10,'r')

% mouse velocity
vbody = diff(body); vbody = sqrt(vbody(:,1).^2 + vbody(:,2).^2);%distance/60hz
abody= diff(vbody);
vbody = conv(vbody,ones(1,12),'same'); %filters/smoothes the velocity data
vscale=vbody/5; % to put on scale more comparable to acceleration
abody = conv(abody,ones(1,12),'same');

 %% prey speed
 
vtarg = diff(targ); vtarg = sqrt(vtarg(:,1).^2 + vtarg(:,2).^2);
atarg= diff(vtarg);
vtarg = conv(vtarg,ones(1,12),'same'); %filters the velcity data
vscaleT=vtarg/5;
atarg= conv(atarg,ones(1,12),'same');
ascale=atarg/2;


%% get angles
headvec = pt2-pt1;
headtheta = getSmoothAngle(headvec)-pi/2;%head angle relative to coordinate system

targvec = targ-head;
targtheta = getSmoothAngle(targvec); %target angle relative to coordinate system

neckvec = head-body;
necktheta = getSmoothAngle(neckvec);

saccadetheta=necktheta-headtheta;
Rtheta=headtheta-targtheta;

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


dt = .1*fps;

vhead = headtheta((1+dt):end)- headtheta(1:end-dt);
vtarg = targtheta((1+dt):end)- targtheta(1:end-dt);


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
