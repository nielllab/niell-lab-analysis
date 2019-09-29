%%%%% Eye Camera Calculations
function [newtheta,newphi,EllipseParams,ExtraParams usegood ngood calibrationR calibrationM scale] = EyeCameraCalc1(numFrames,Pointsx,Pointsy,Likelihood, psfilename, eyethresh)

% Inputs:
%   Vid1 - 3D grayscale array of video frames
%   Pointsx - The X component of Deeplabcut tracking
%   Pointsy - The Y component of Deeplabcut tracking
%    eyethresh - likelihood threshold for including eye points

% Outputs:
%   newtheta - the horizontal angle difference between Camera center and
%   pupil center.
%   newphi - the vertical angle between the camera center and the pupil
%   center.
%   EllipseParams - Parameters for the ellipse fit (# of frames by 7)
%       Column1 (X0_in) - the x coordinate of the ellipse
%       Column2 (Y0_in) - the y coordinate of the ellipse
%       Column3 (long_axis/2) - Radius of the long axis
%       Column4 (short_axis/2) - Radius of the short axis
%       Column5 (angleToX) - angle from long axis to horizontal plane
%       Column6 (angleFromX) - angle from hoizontal plane to long axis
%       Column7 (phi) - the tilt of the ellipse
%   ExtraParams - Extra parameters to explore
%  usegood - are all 8 points good?
%  ngood = number of good points
% calibrationR, calibrationM = correlation coeff and slope for calibration values (to be used as diagnostic)
% scale = scale factor to convert pix to deg based on ellipticity

if ~exist('eyethresh','var')
    eyethresh = 0.97;
end

if exist('psfilename','var')
    savePDF=1;
else
    savePDF=0;
end



%% Calculate Ellipse Fits
% figure;
% axis tight manual
% ax = gca;
% % ax.NextPlot = 'replaceChildren';
% axis([1 size(Vid1,1) 1 size(Vid1,2)])


good = Likelihood>=eyethresh; %likelihood of all 8 pts must be >=.95 
ngood = sum(good,2);
usegood = ngood>=8;


EllipseParams = zeros(numFrames,7);
ExtraParams = zeros(numFrames,6);


parfor v=1:(numFrames)
    if usegood(v)==1 && (~any((Pointsx(v,:)<10) | (Pointsy(v,:)<10)| (isnan(Pointsx(v,:)))))
        e_t = fit_ellipse2(Pointsx(v,:),Pointsy(v,:));
        if isempty(e_t.status)==1 &&  ((e_t.short_axis/2)/(e_t.long_axis/2) >.5)
            EllipseParams(v,:)=[e_t.X0_in, e_t.Y0_in, e_t.long_axis/2, e_t.short_axis/2,  e_t.angleToX*pi/180, e_t.angleFromX*pi/180, e_t.phi];
            ExtraParams(v,:) = [e_t.X0, e_t.Y0, e_t.a, e_t.b, e_t.cos_phi, e_t.sin_phi];
        end
    end
end
fprintf('done \n')
close all;


figure
plot(Likelihood); ylim([0 1])
ylabel('likelihood'); title('likelihood all eyepoints')
if savePDF
    set(gcf, 'PaperPositionMode', 'auto');  print('-dpsc',psfilename,'-append'); close(gcf)
end


%%% plot number of good points
figure
subplot(2,2,3)
plot(ngood); ylim([0 9]);
ylabel('# good eyepoints'); title(sprintf('%0.3f good thresh %0.2f',mean(usegood),eyethresh))
if savePDF
  set(gcf, 'PaperPositionMode', 'auto'); print('-dpsc',psfilename,'-append'); close(gcf)
end


%%% plot values for good points

subplot(2,2,1);
plot(Pointsx(usegood,:)); ylabel('x'); title('trace of  good timepoints')
subplot(2,2,2);
plot(Pointsy(usegood,:)); ylabel('y');

efitb = find(EllipseParams(:,1)==0);
EllipseParams(efitb,:)=NaN;
ExtraParams(efitb,:) = NaN;
% EllipseParams = fillmissing(EllipseParams,'linear',1);  %%% dangerous - not sure we want to do this
% ExtraParams = fillmissing(ExtraParams,'linear',1);

e_thresh=0.9; %%% ellipticity threshold

%%% plot ellipticity histogram
EllRange= (EllipseParams(usegood,4)./EllipseParams(usegood,3))
subplot(2,2,4)
hist(EllRange);
title(sprintf('ellipticity thresh = %.2f',e_thresh))


if savePDF
    set(gcf, 'PaperPositionMode', 'auto');  print('-dpsc',psfilename,'-append'); close(gcf)
end


%%  Calc Camera Center
R = linspace(0,2*pi,100);
list = find(EllipseParams(:,4)./EllipseParams(:,3)<e_thresh); %randi([1 size(EllipseParams,1)],50);%  1:size(EllipseParams,1); %
R = linspace(0,2*pi,100); 
list = find(usegood & EllipseParams(:,4)./EllipseParams(:,3)<e_thresh); %randi([1 size(EllipseParams,1)],50);%  1:size(EllipseParams,1); %
A = [cos(EllipseParams(list,5)),sin(EllipseParams(list,5))]; %%% cosw  + sinw
b=diag(A*EllipseParams(list,1:2)');
CamCent=(A'*A)\(A'*b)   %%% camcent*A = EllipseParams(list,1:2)*A

w= EllipseParams(:,5);

figure
subplot(2,2,1)
hold on
for i = 1:length(list);
    plot(EllipseParams(list(i),1) + [-5*cos(w(list(i))) 5*cos(w(list(i)))], EllipseParams(list(i),2) + [-5*sin(w(list(i))) 5*sin(w(list(i)))])
end
axis equal; plot(CamCent(1),CamCent(2),'r*','Markersize',8);
title('eye axes relative to center');

% scale = nansum(sqrt(1-(EllipseParams(:,4)./EllipseParams(:,3)).^2)'*vecnorm([EllipseParams(:,1)';EllipseParams(:,2)']-CamCent,2,1)')/...
%     nansum(1-(EllipseParams(:,4)./EllipseParams(:,3)).^2);
Ellipticity = EllipseParams(list,4)./EllipseParams(list,3);  %%% use all good points, since we want a range of ellipticity
scale = nansum(sqrt(1-(Ellipticity).^2).*vecnorm(EllipseParams(list,1:2)'-CamCent,2,1)')./nansum(1-(Ellipticity).^2);
theta = asin((EllipseParams(:,1)-CamCent(1))*1/scale); %in radians
thetad = asind((EllipseParams(:,1)-CamCent(1))*1/scale); %in deg
phi = asind((EllipseParams(:,2)-CamCent(2))./cos(theta)*1/scale); %in deg



%% calculate one example image
i=50;
w = EllipseParams(i,5); L=EllipseParams(i,3); l=EllipseParams(i,4); x_C=EllipseParams(i,1); y_C=EllipseParams(i,2);
Rotation1 = [cos(w),-sin(w);sin(w),cos(w)];
L1 = [L,0;0,l];
C1 = [x_C; y_C];
q = [cos(R);sin(R)];
q_star = Rotation1*L1*q+C1;
qcirc = [L/scale,0;0,L/scale]*q;
qcirc2 = [L,0;0,L]*q+CamCent;

theta2 = real(asin((EllipseParams(i,1)-CamCent(1))*1/scale));
phi2 = real(asin((EllipseParams(i,2)-CamCent(2))./cos(theta2)*1/scale));

newCent = scale*[sin(theta2); sin(phi2).*cos(theta2)]+CamCent;
PointsRot = newCent + scale*[cos(theta2), 0; -sin(theta2)*sin(phi2), cos(phi2)]*qcirc;

subplot(2,2,2)
%imagesc(Vid1(:,:,i)); colormap gray;
axis equal off;
plot(q_star(1,:),q_star(2,:),'g','LineWidth',2);hold on; scatter(EllipseParams(i,1),EllipseParams(i,2),100,'og');
scatter(Pointsx(i,:),Pointsy(i,:),100,'.m'); %DLC pts
scatter(CamCent(1),CamCent(2),100,'or'); % theoretical cam center
scatter(qcirc2(1,:),qcirc2(2,:),50,'.r'); %theoretical circle of camera axis (eye is straight ahead)
scatter(newCent(1),newCent(2),100,'xb'); % new center 
scatter(PointsRot(1,:),PointsRot(2,:),100,'.b'); %transformation from red to green, good if matches green
axis equal; %axis([250 400 250 400])
title(sprintf('omega = %.2f',EllipseParams(i,5)*180/pi))

%  scatter(EyeShift(i,2),EyeShift(i,3),100,'.y')




%% Check Calibration

subplot(2,2,3)
xvals = vecnorm(EllipseParams(usegood,1:2)' - CamCent); yvals = scale*sqrt(1-(EllipseParams(usegood,4)./EllipseParams(usegood,3)).^2);
plot(xvals,yvals,'.');
[calibrationR calibrationM b] = regression(xvals,yvals')
axis equal;hold on; xlabel('pupil camera dist'); ylabel('scale * ellipticity');
plot(linspace(0,250),linspace(0,250),'r')
title(sprintf('Scale=%0.1f r=%0.1f m=%0.1f',scale,calibrationR,calibrationM))
axis([0 50 0 50])

%%% calibration fit test
% x1=vecnorm(EllipseParams(:,1:2)' - CamCent)';
% y1=scale*sqrt(1-(EllipseParams(:,4)./EllipseParams(:,3)).^2);
% fitvars = polyfit(x1,y1, 1);
% m = fitvars(1);
% c = fitvars(2);

subplot(2,2,4)

delta = (CamCent-EllipseParams(:,1:2)');

plot(vecnorm(delta(:,usegood),2,1),((delta(1,usegood)'.*cos(EllipseParams(usegood,5)))+(delta(2,usegood)'.*sin(EllipseParams(usegood,5))))./vecnorm(delta(:,usegood),2,1)','.')
hold on
plot(vecnorm(delta(:,list),2,1),((delta(1,list)'.*cos(EllipseParams(list,5)))+(delta(2,list)'.*sin(EllipseParams(list,5))))./vecnorm(delta(:,list),2,1)','r.')
title('Camera Center Calibration')
ylabel('abs([PC-EC]).[cosw;sinw]');
xlabel('abs(PC-EC)');
legend('all points','list points')
if savePDF
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    close(gcf)
end



%% Theta and Phi

newtheta = (real(thetad));
newphi = (real(phi));

figure; subplot(121);
plot(diff(movmean(newtheta,10)),'LineWidth',2); %hold on; plot(diff(newlongang),'LineWidth',1)
title('Change in smoothed theta')
subplot(122);
plot(diff(movmean(newphi,10)),'LineWidth',2);  % hold on; plot(diff(newlongang),'LineWidth',1)
title('Change in smoothed phi')
% if savePDF
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfilename,'-append');
%     close(gcf)
% end






