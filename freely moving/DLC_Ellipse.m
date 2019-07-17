clear; close all;

[f,p] = uigetfile({'*.avi';},'choose eye EyeFilt');
% Read in Eye Tracking Video
if f ~=0
    TempVid = VideoReader(fullfile(p,f));
    frame=1; k=1;
    while hasFrame(TempVid)
        EyeVid(:,:,frame) = rgb2gray(readFrame(TempVid));
        if mod(frame,100)==0
            fprintf('frame = %d\n',frame)
        end
        frame=frame+1;
        
    end
end
fprintf('Done Reading in Video \n');
EyeVid2 = EyeVid;%movmean(EyeVid(:,:,1:5000),Wind,3);
% load(fullfile(p,'EllipsePoints2.mat'));
vidsize = 1:size(EyeVid2,3); %1950:2850; %
[f2,p2] = uigetfile({'*.csv';},'choose DLC Points File');
Pointsxy = csvread(fullfile(p2,f2),3,0);
Pointsx = Pointsxy(vidsize,[2:3:end]); Pointsy = Pointsxy(vidsize,[3:3:end]); 
LH =  Pointsxy(vidsize,[4:3:end]); 
j=1;
for i=1:size(LH,1)
    if isempty(find(LH(i,:)<.1))~=1
        bdfit(j,1) = i; j=j+1;
    end
end
EyeVid2(:,:,bdfit)=[]; Pointsx(bdfit,:)=[]; Pointsy(bdfit,:)=[];
% [f3,p3] = uigetfile({'*.csv';},'choose DLC Points File');
% EyeShift = csvread(fullfile(p3,f3),3,0);
% EyeShiftx = EyeShift(:,[2:3:end]); EyeShifty = EyeShift(:,[3:3:end]);
% EyeShiftx = POSITION_X; EyeShifty = POSITION_Y;
% Pointsx=Pointsx(1:5000,:);% movmean(Pointsx(1:5000,:),Wind,1);
% Pointsy=Pointsy(1:5000,:); %movmean(Pointsy(1:5000,:),Wind,1);
%%
i = randi([1 size(EyeVid,3)],1);
figure; 
imagesc(EyeVid2(:,:,i));
uiwait(msgbox('Draw a region with left button down.  Double click inside to finish it'));
h = imrect;
vertices = wait(h);
vertices = round(vertices);
CropEye = double(EyeVid2(vertices(2):vertices(2)+vertices(4),vertices(1):vertices(1)+vertices(3),:));
mx = (max(max(max(CropEye,3))));
CropEye = CropEye./mx;
close all; 
fprintf('y = %d : %d \nx = %d : %d\n',vertices(2),vertices(2)+vertices(4),vertices(1),vertices(1)+vertices(3))

%% bright spot
radius = 5; decomposition = 0;
se = strel('disk', radius, decomposition);
% figure;
centroid = zeros(size(EyeVid2,3),2);
rad = zeros(size(EyeVid2,3),1);
test = zeros(size(EyeVid2,3),1);

for i=1:size(EyeVid2,3)
    temp = CropEye(:,:,i);
%     B = imgaussfilt(temp,4);
    [BW,maskedImage] = segmentImage(temp);
%     BW1 = imopen(B, se);
%     BW = imbinarize(B, 'adaptive', 'Sensitivity', 0.420000, 'ForegroundPolarity', 'bright');
%     BW2 = edge(BW,'Canny',[.90 .990]);
%     BW2 = imclearborder(BW);
    try
    [center,radii,metric] = imfindcircles(BW,[3 12],'ObjectPolarity','bright','Sensitivity',0.8);
    catch
        test(i) = i;
    end
    [~,idx] = max(metric); % pick the circle with best score
    if isempty(idx)==0
        centroid(i,:) = center(idx,:);
        rad(i,:) = radii(idx);
    end
    imshowpair(temp,BW); hold on; scatter(centroid(i,1),centroid(i,2),100,'b.');hold off;
    drawnow limitrate;
end

Centroid2 = centroid; Rad2 = rad;
zrs = find(Centroid2(:,1)==0);
Centroid2(zrs,:)= []; Rad2(zrs)=[];
EyeVid2(:,:,zrs)=[]; Pointsx(zrs,:)=[]; Pointsy(zrs,:)=[];
Centroid2(:,1) = Centroid2(:,1)+ vertices(1); Centroid2(:,2)=Centroid2(:,2)+vertices(2); 
Centroid2 = movmean(Centroid2,5,1);
Cen_Var = max(std(Centroid2,[],1));
Cen_mean = mean(Centroid2,1);
fprintf('Done\n')
%% Shift Video to Stabalize
% shiftxy(1,:) = nanmean(EyeShiftxb,1); shiftxy(2,:) = nanmean(EyeShiftyb,1);
% shiftx=Cen_mean(1,1)-Centroid2(:,1);
shiftx=table2array(Stable2(:,1));
shifty=table2array(Stable2(:,2));
% shifty=Cen_mean(1,2)-Centroid2(:,2);
% shiftx(BdFrames2,:)=0; shifty(BdFrames2,:)=0;

% PupCentx = EllipseParams(:,1); PupCenty=EllipseParams(:,2);
parfor i =1: size(EyeVid2,3)-1
    EyeVid2(:,:,i) = imtranslate(EyeVid2(:,:,i),[shiftx(i,1) shifty(i,1)]);
%     PupCentx(i,1) = PupCentx(i,1) + shiftx(i,1);
%     PupCenty(i,1) = PupCenty(i,1) + shifty(i,1);
    Pointsx(i,:) = Pointsx(i,:) + shiftx(i,1);
    Pointsy(i,:) = Pointsy(i,:) + shifty(i,1);   
end
fprintf('done \n')
% EllipseParams(:,1) = PupCentx(:,1);
% EllipseParams(:,2) = PupCenty(:,1);
%% Calculate Ellipse Fits
figure;
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';
axis([1 size(EyeVid2,1) 1 size(EyeVid2,2)] )
EllipseParams = zeros(size(EyeVid2,3),7);
% imshow(EyeVid2(:,:,1)); axis equal off; hold on;

parfor v=1:size(EyeVid2,3)
    if ~any((Pointsx(v,:)<10 | Pointsy(v,:)<10)| isnan(Pointsx(v,:)))
        e_t = fit_ellipse(Pointsx(v,:),Pointsy(v,:),ax);
        if isempty(e_t.status)==1 &&  ((e_t.short_axis/2)/(e_t.long_axis/2) >.5)
            EllipseParams(v,:)=[e_t.X0_in, e_t.Y0_in, e_t.long_axis/2, e_t.short_axis/2,  e_t.angleToX*pi/180, e_t.angleFromX*pi/180, e_t.phi];
        end
    end
end
fprintf('done \n')
close all;

%%
% BdFrames=find(EllipseParams(:,1)<=10);
% 
% EllipseParams(BdFrames,:)=[];
% EyeVid2(:,:,BdFrames)=[];
% Pointsx(BdFrames,:) =[]; Pointsy(BdFrames,:)=[];
% EyeShiftx(BdFrames,:)=[]; EyeShifty(BdFrames,:)=[];

% shiftxy(1,:) = median(EyeShiftx,1); shiftxy(2,:) = median(EyeShifty,1);
% shiftx=shiftxy(1,:)-EyeShiftx;
% shifty=shiftxy(2,:)-EyeShifty;
% BdFrames2 = find(abs(shiftx)>30 | abs(shifty)>30);
% EyeShiftxb=EyeShiftx; EyeShiftyb=EyeShifty;
% EyeShiftxb(BdFrames2,:)=NaN; EyeShiftyb(BdFrames2,:)=NaN;
% EyeShiftx(BdFrames2,:)=0; EyeShifty(BdFrames2,:)=0;

% EllipseParams(BdFrames2,:)=[];
% EyeVid2(:,:,BdFrames2)=[];
% Pointsx(BdFrames2,:) =[]; Pointsy(BdFrames2,:)=[];
% EyeShiftx(BdFrames2,:)=[]; EyeShifty(BdFrames2,:)=[];

% Wind=5;

% EllipseParams(:,1) = movmean(EllipseParams(:,1),Wind);
% EllipseParams(:,2) = movmean(EllipseParams(:,2),Wind);
% EllipseParams(:,5) = medfilt1(EllipseParams(:,5),5);
% EllipseParams = medfilt1(EllipseParams,Wind,[],1);
%%



%% Plot Ellipse Fits on Image
figure;
axis tight manual
ax = gca;
ax.NextPlot = 'replaceChildren';
axis([1 size(EyeVid2,1) 1 size(EyeVid2,2)] )
% EllipseParams = zeros(size(EyeVid2,3),6);

for v=1000:size(EyeVid2,3)
%     v=10135;z
    imagesc(EyeVid2(:,:,v)); colormap gray; axis equal off; hold on;
    scatter(Pointsx(v,:),Pointsy(v,:),100,'.r')
%     if ~any((Pointsx(v,:)<10 | Pointsy(v,:)<10)| isnan(Pointsx(v,:)))
%         e_t = fit_ellipse(Pointsx(v,:),Pointsy(v,:),ax);
% %         scatter(EyeShiftx(v,:),EyeShifty(v,:),100,'.r');
% %         if isempty(e_t.status)==1 &&  ((e_t.short_axis/2)/(e_t.long_axis/2) >.6)
% %             EllipseParams(v,:)=[e_t.X0_in, e_t.Y0_in, e_t.long_axis/2, e_t.short_axis/2,  e_t.angleToX*pi/180,e_t.phi];
% %         end
%     end
    hold off;
    title(sprintf('Frame = %d',v));
    drawnow limitrate
%     pause
end
fprintf('done \n')
% close all;
%%  Calc Camera Center
R = linspace(0,2*pi,100);
list = find(EllipseParams(:,4)./EllipseParams(:,3)<.98); %randi([1 size(EllipseParams,1)],50);%  1:size(EllipseParams,1); %
A = [cos(EllipseParams(list,5)),sin(EllipseParams(list,5))]; 
b=diag(A*EllipseParams(:,1:2)');
CamCent=(A'*A)\A'*b;
% 
% scale = nansum(sqrt(1-(EllipseParams(:,4)./EllipseParams(:,3)).^2)'*vecnorm([EllipseParams(:,1)';EllipseParams(:,2)']-CamCent,2,1)')/...
%     nansum(1-(EllipseParams(:,4)./EllipseParams(:,3)).^2);
Ellipticity = EllipseParams(list,4)./EllipseParams(list,3);
scale = nansum(sqrt(1-(Ellipticity).^2).*vecnorm(EllipseParams(list,1:2)'-CamCent,2,1)')./nansum(1-(Ellipticity).^2);

theta = asin((EllipseParams(:,1)-CamCent(1))*1/scale);
phi = asin((EllipseParams(:,2)-CamCent(2))./cos(theta)*1/scale);

%%
i=1350;
w = EllipseParams(i,5); L=EllipseParams(i,3); l=EllipseParams(i,4); x_C=EllipseParams(i,1); y_C=EllipseParams(i,2);
Rotation1 = [cos(w),-sin(w);sin(w),cos(w)];
L1 = [L,0;0,l];
C1 = [x_C; y_C];
q = [cos(R);sin(R)];
q_star = Rotation1*L1*q+C1;
qcirc = [L/scale,0;0,L/scale]*q;
qcirc2 = [L,0;0,L]*q+CamCent;

theta2 = asin((EllipseParams(i,1)-CamCent(1))*1/scale);
phi2 = asin((EllipseParams(i,2)-CamCent(2))./cos(theta2)*1/scale);

newCent = scale*[sin(theta2); sin(phi2).*cos(theta2)]+CamCent;
PointsRot = newCent + scale*[cos(theta2), 0; -sin(theta2)*sin(phi2), cos(phi2)]*qcirc;

figure; imagesc(EyeVid2(:,:,i)); colormap gray; axis equal off; hold on;
plot(q_star(1,:),q_star(2,:),'g','LineWidth',2); scatter(EllipseParams(i,1),EllipseParams(i,2),100,'og');
scatter(Pointsx(i,:),Pointsy(i,:),100,'.m');
scatter(CamCent(1),CamCent(2),100,'or'); scatter(qcirc2(1,:),qcirc2(2,:),50,'.r');
scatter(newCent(1),newCent(2),100,'xb')
scatter(PointsRot(1,:),PointsRot(2,:),100,'.b');

% scatter(EyeShift(i,2),EyeShift(i,3),100,'.y')

%%
figure;
plot(vecnorm(EllipseParams(:,1:2)' - CamCent),scale*sqrt(1-(EllipseParams(:,4)./EllipseParams(:,3)).^2),'.')
axis equal;hold on; 
plot(linspace(0,250),linspace(0,250),'r')
title('Scale Factor Calibration')
axis([0 100 0 100])

x1=vecnorm(EllipseParams(:,1:2)' - CamCent)';
y1=scale*sqrt(1-(EllipseParams(:,4)./EllipseParams(:,3)).^2);
fitvars = polyfit(x1,y1, 1);
m = fitvars(1);
c = fitvars(2);

delta = (CamCent-EllipseParams(:,1:2)');
figure; 
plot(vecnorm(delta,2,1),((delta(1,:)'.*cos(EllipseParams(:,5)))+(delta(2,:)'.*sin(EllipseParams(:,5))))./vecnorm(delta,2,1)','.')
% axis equal;hold on; 
% plot(linspace(0,100),linspace(0,100),'r');
% plot(linspace(-100,0),linspace(100,0),'r');
title('Camera Center Calibration')
ylabel('abs([PC-EC]).[cosw;sinw]');
xlabel('abs(PC-EC)')

%% Theta and Phi
dtheta = diff(theta);
dtheta(dtheta>pi/2) = dtheta(dtheta>pi/2)-pi;
dtheta(dtheta<-pi/2) = dtheta(dtheta<-pi/2)+pi;
newtheta = theta(1) + cumsum(dtheta);
dphi = diff(phi);
dphi(dphi>pi/2) = dphi(dphi>pi/2)-pi;
dphi(dphi<-pi/2) = dphi(dphi<-pi/2)+pi;
newphi = phi(1) + cumsum(dphi);

figure; subplot(121);% yyaxis left; 
plot(real(newtheta)*180/pi); hold on; plot(diff(headang))
% yyaxis right; plot(Vx)
subplot(122); %yyaxis left; 
plot(real(newphi)*180/pi);   hold on; plot(diff(headang))
%  yyaxis right; plot(Vy)
%% angle from X
dtheta = diff(EllipseParams(:,5));
dtheta(dtheta>pi/2) = dtheta(dtheta>pi/2)-pi;
dtheta(dtheta<-pi/2) = dtheta(dtheta<-pi/2)+pi;
newtheta = theta(1) + cumsum(dtheta);
% figure; yyaxis left;
% plot(EllipseParams(:,5)*180/pi); 
% yyaxis right; plot(Vx)
%% Radius
figure; 
plot(EllipseParams(:,3))

%%
% n_m = 1./(sqrt(L^2.*sin(R).^2+l^2.*cos(R).^2)).*[L*sin(w).*sin(R)-l*cos(w).*cos(R);L*cos(w).*sin(R)-l*sin(w).*cos(R)];


% %%
% psi=0; phi=atan((CamCent(2)-EllipseParams(:,2))./(CamCent(1)-EllipseParams(:,1)));
% theta = acos(1./sqrt((CamCent(1)-EllipseParams(:,1)).^2+(CamCent(2)-EllipseParams(:,2)).^2+1));
% % psi=0; phi=atan((CamCent(2)-C1(2))/(CamCent(1)-C1(1))); theta = acos(1/sqrt((CamCent(1)-C1(1)).^2+(CamCent(2)-C1(2)).^2+1));
% R_phi = [1, 0, 0; 0, cos(phi), sin(phi); 0, sin(phi), cos(phi)];
% R_theta = [cos(theta),0,-sin(theta); 0,1,0; sin(theta),1,cos(theta)];
% R_psi = [cos(psi),sin(psi),0; sin(psi),cos(psi),0; 0,0,1];
% R_eye = R_phi*R_theta*R_psi;
% 
% gaze = [sin(theta); sin(phi)*cos(theta); -cos(phi)*cos(theta)];
%%
% PupCent = scale*[sin(theta);sin(phi)*cos(theta)]+CamCent;
% points = PupCent+scale*[cos(theta(1)),0; -sin(theta(1))*sin(phi(1)), cos(phi(1))]*[q_star(1,:);q_star(2,:)];
% 
% 
% gazeE = [points;-ones(1,size(points,2))];
% R3 = R_eye*gazeE;
% r = sqrt(points(1,:).^2+points(2,:));
% 
% points2 = r./sqrt(1-(EllipseParams(1,4)./EllipseParams(1,3)).^2) .* [-sin(phi); cos(phi)*sin(theta)];
%%
figure; sphere; colormap gray; hold on; axis equal;

scatter3(-cos(0)*cos(0),sin(0)*cos(0),sin(0),1000,'.r')
scatter3(-cos(real(phi)).*cos(real(theta)),sin(real(phi)).*cos(real(theta)),sin(real(theta)),100,'.b')

% 
% plot3(sin(theta).*cos(phi),sin(theta).*sin(phi),cos(theta),'.'); axis([-1 1 -1 1 -1 1]);
% scatter3(CamCent(1)/norm(CamCent),CamCent(2)/norm(CamCent),0,100,'.g');
% scatter3(R3(1,:)/norm(R3),R3(2,:)/norm(R3),R3(3,:)/norm(R3),100,'r');











