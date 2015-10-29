function [xc vxc h r] = analyzeCapture(fname)
if isempty(fname)
[f p] = uigetfile('*.*','behav data');
fname = fullfile(p,f);
end
try
data = xlsread(fname);
catch
    data =dlmread(fname);
end

missing = isnan(mean(data,2));
data = data(~missing,:);

figure

%%% read in points
body = data(:,5:6);
%body(245:246,1)=350; body(245:246,2)=530; %%%% fill in Nans
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

%%% plot head and target trajectories
subplot(2,3,1)
plot(head(:,1),head(:,2));
hold on
plot(targ(:,1),targ(:,2),'g')
axis([ 0 2000 0 1200])


%%% distance between head and target
r = sqrt((head(:,1) - targ(:,1)).^2 + (head(:,2) - targ(:,2)).^2);
subplot(2,3,2)
plot(r); ylabel('distance to target')
xlim([1 length(r)]); ylim([0 2000])

%%%% mouse speed
vbody = diff(body); vbody = sqrt(vbody(:,1).^2 + vbody(:,2).^2);
vbody = conv(vbody,ones(1,10),'same');
% figure
% plot(vbody); ylabel('body speed')

%%% get angles
headvec = pt2-pt1;
headtheta = getSmoothAngle(headvec)-pi/2;

targvec = targ-head;
targtheta = getSmoothAngle(targvec);

neckvec = head-body;
necktheta = getSmoothAngle(neckvec);


% plot((necktheta-headtheta)*180/pi)
% ylabel('head angle wrt body'); ylim([-45 45])


%%% head angle and target position
subplot(2,3,3)
plot(mod(headtheta*180/pi,360)); hold on; plot(mod(targtheta*180/pi,360),'g');
legend('head','targ'); ylabel('theta'); xlabel('time');
xlim([1 length(headtheta)])

% %%% angle between head and target
% figure
%  plot((headtheta-targtheta)*180/pi,'k','Linewidth',2); hold on
%  plot([1 length(targtheta)],[0 0 ],'r');
% ylabel('head wrt target')

%%% calculate correlation between head angle and target location
opt = 'unbiased';
[xc lags] = xcorr(headtheta,targtheta,opt); t0= find(lags==0);
cc = corrcoef(headtheta,targtheta);
subplot(2,3,4)
hold on; 
[headxc] = xcorr(headtheta,headtheta,opt);
hold on; plot(lags,headxc/headxc(t0),'Linewidth',2);
[targxc] = xcorr(targtheta,targtheta,opt);
plot(lags,targxc/targxc(t0),'g','Linewidth',2)
plot(lags,xc*cc(1,2)/xc(t0),'k','Linewidth',2); 
plot([0 0], [0 1.2],'r');
legend('head','targ','head vs targ')
xlabel('lag');
ylabel('correlation')
headxc = headxc(lags>=-100 & lags<=100);
targxc = targxc(lags>=-100 & lags<=100);
xc = xc(lags>=-100 & lags<=100);

subplot(2,3,4); hold off
xc=lagCircRMS(headtheta,targtheta,-300:300);
plot(-300:300,xc,'k');
hold on
plot(-300:300,lagCircRMS(headtheta,headtheta,-300:300),'b')
plot(-300:300,lagCircRMS(targtheta,targtheta,-300:300),'g')
axis([-300 300 -1 1])


%%% correlation between target movement and head movement
dt = 10;

vhead = headtheta((1+dt):end)- headtheta(1:end-dt);
vtarg = targtheta((1+dt):end)- targtheta(1:end-dt);
% subplot(2,3,4); hold off
% plot(vhead,'b');
% hold on
% plot(vtarg,'g')
% xlim([1 length(vtarg)])

%plot(vhead); hold on; plot(vtarg,'g')
[vxc lags] = xcorr(vhead,vtarg,'coeff');
vxc = conv(vxc,ones(1,100),'same')/100;
vxc = vxc(lags>=-300 & lags<=300);
lags = lags(lags>=-300 & lags<=300);
subplot(2,3,5)
plot(lags/60,vxc,'Linewidth',2);
xlabel('secs')
ylabel('xcorr'); legend('head vs targ angular vel')
axis([ -5 5 -0.1 0.2 ])

%%% how far off? head angle and distance
% figure
% plot(r, 'Linewidth',2); hold on; plot((headtheta-targtheta)*180/pi,':','Linewidth',2)

%%% location of target relative to animal's view
x = r.*cos(headtheta-targtheta); y = r.*sin(headtheta-targtheta);
subplot(2,3,6)
plot(y,x); hold on; plot(0,0,'r*'); axis([-1600 1600 -1600 1600])

% figure 
bins = -1025:50:1025;
h = hist2(x,y,bins,bins);
h= h/sum(h(:));
% imagesc(flipud(h))
% hold on
% plot(floor(length(bins)/2) ,floor(length(bins)/2) + 1,'ro','Linewidth',2)


