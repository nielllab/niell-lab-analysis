%%%%% Eye Camera Calculations
function [newtheta,newphi,EllipseParams,ExtraParams usegood ngood calibrationR calibrationM scale] = EyeCameraCalc1(numFrames,Pointsx,Pointsy,Likelihood, psfilename, eyethresh)
% convert DeepLacCut (DLC) pupil circumference points into theta/phi rotation of the eye
%
% Overview
% 1) first, fits ellipse to pupil points
% 2) then compute camera center and scale factor from ellipses as delineated in Supplemental Material, p.26 
% 3) finally calculates theta/phi fo rech timepoint
% Steps 2,3 follow Wallace et al 2013, Supplemental Information p.26

% Inputs:
%   Vid1 - 3D grayscale array of video frames (used only for figures)
%   Pointsx - The X component of Deeplabcut tracking
%   Pointsy - The Y component of Deeplabcut tracking
%    eyethresh - likelihood threshold for including eye points
%  note : last 3 are result of DLC tracking of 8 evenly spaced point around the pupil

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

if ~exist('eyethresh','var')  %%% DLC likelihood threshold for pupil points
    eyethresh = 0.97;
end

if exist('psfilename','var')  %% save pdfs?
    savePDF=1;
else
    savePDF=0;
end



%% Calculate Ellipse Fits

%%% select timepoints that have all 8 good pupil points from DLC
good = Likelihood>=eyethresh; %likelihood of all 8 pts must be >=.95 
ngood = sum(good,2);
usegood = ngood>=8;


EllipseParams = zeros(numFrames,7);
ExtraParams = zeros(numFrames,6);

%%% fit ellipse for all good timepoints
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

% bad fits return zeros, set to NaNs
efitb = find(EllipseParams(:,1)==0);
EllipseParams(efitb,:)=NaN;
ExtraParams(efitb,:) = NaN;

e_thresh=0.9; %%% ellipticity threshold - pupils that are not slightly elliptical can gve spurious results


%%% some diagnostic plots

%%% likelihood of all points over time
figure
plot(Likelihood); ylim([0 1])
ylabel('likelihood'); title('likelihood all eyepoints')
if savePDF
    set(gcf, 'PaperPositionMode', 'auto');  print('-dpsc',psfilename,'-append'); close(gcf)
end


%%% plot number of good points over time
figure
subplot(2,2,3)
plot(ngood); ylim([0 9]);
ylabel('# good eyepoints'); title(sprintf('%0.3f good thresh %0.2f',mean(usegood),eyethresh))
if savePDF
  set(gcf, 'PaperPositionMode', 'auto'); print('-dpsc',psfilename,'-append'); close(gcf)
end


%%% plot x/y values for good points
subplot(2,2,1);
plot(Pointsx(usegood,:)); ylabel('x'); title('trace of  good timepoints')
subplot(2,2,2);
plot(Pointsy(usegood,:)); ylabel('y');


%%% plot ellipticity histogram
EllRange= (EllipseParams(usegood,4)./EllipseParams(usegood,3))
subplot(2,2,4)
hist(EllRange);
title(sprintf('ellipticity thresh = %.2f',e_thresh))



%% calibrate eye camera based on ellipse fits (see Wallace 2013 supplement)

%%%  Calc Camera Center 
%%% least-squares solution of the equation at the bottom of Wallace p. 27)

list = find(usegood & EllipseParams(:,4)./EllipseParams(:,3)<e_thresh); % only use good timepoints with ellipticity < e_threhs
A = [cos(EllipseParams(list,5)),sin(EllipseParams(list,5))]; %%% cosw  + sinw
b=diag(A*EllipseParams(list,1:2)');
CamCent=(A'*A)\(A'*b)   %%% camcent*A = EllipseParams(list,1:2)*A

w= EllipseParams(:,5);  %%% ellipse angle, omega in Wallace 2013

% a set of anlges for plotting
R = linspace(0,2*pi,100); 


%%% plot pupil ellipse axis for each timepoint, relative to camera center location
%%% ellipse axes should be orthogonal to camera center
figure
subplot(2,2,1)
hold on
for i = 1:length(list);
    plot(EllipseParams(list(i),1) + [-5*cos(w(list(i))) 5*cos(w(list(i)))], EllipseParams(list(i),2) + [-5*sin(w(list(i))) 5*sin(w(list(i)))])
end
axis equal; plot(CamCent(1),CamCent(2),'r*','Markersize',8);
title('eye axes relative to center');


%%% calculate scale factor f/z0
%%% Wallace, p28
Ellipticity = EllipseParams(list,4)./EllipseParams(list,3);  %%% use all good points, since we want a range of ellipticity
scale = nansum(sqrt(1-(Ellipticity).^2).*vecnorm(EllipseParams(list,1:2)'-CamCent,2,1)')./nansum(1-(Ellipticity).^2);

%%% calculate theta and phi
%%% invert equation at bottom of Wallace, p. 24
theta = asin((EllipseParams(:,1)-CamCent(1))*1/scale); %in radians
thetad = asind((EllipseParams(:,1)-CamCent(1))*1/scale); %in deg
phi = asind((EllipseParams(:,2)-CamCent(2))./cos(theta)*1/scale); %in deg


%%% calculate one example image
%%% transform a unit circle following equations on Wallace p. 23
i=50; % timepoint to use
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

%%% plot the transformed circle
subplot(2,2,2)
axis equal off;
plot(q_star(1,:),q_star(2,:),'g','LineWidth',2);hold on; scatter(EllipseParams(i,1),EllipseParams(i,2),100,'og');
scatter(Pointsx(i,:),Pointsy(i,:),100,'.m'); %DLC pts
scatter(CamCent(1),CamCent(2),100,'or'); % theoretical cam center
scatter(qcirc2(1,:),qcirc2(2,:),50,'.r'); %theoretical circle of camera axis (eye is straight ahead)
scatter(newCent(1),newCent(2),100,'xb'); % new center 
scatter(PointsRot(1,:),PointsRot(2,:),100,'.b'); %transformation from red to green, good if matches green
axis equal; %axis([250 400 250 400])
title(sprintf('omega = %.2f',EllipseParams(i,5)*180/pi))


%% Check Calibration
%%% we do two tests to make sure calibration works, verify the two least-squares fits

%%% first check least-squares fit for scale factor
%%% if it worked correctly, left and right side of equation in middle of p28 should be equal
%%% i.e. these points should lie on identity line
%%% we also keep values from a regression fit to this, to use as diagnostic
subplot(2,2,3)
xvals = vecnorm(EllipseParams(usegood,1:2)' - CamCent); yvals = scale*sqrt(1-(EllipseParams(usegood,4)./EllipseParams(usegood,3)).^2);
plot(xvals,yvals,'.');
[calibrationR calibrationM b] = regression(xvals,yvals')
axis equal;hold on; xlabel('pupil camera dist'); ylabel('scale * ellipticity');
plot(linspace(0,250),linspace(0,250),'r')
title(sprintf('Scale=%0.1f r=%0.1f m=%0.1f',scale,calibrationR,calibrationM))
axis([0 50 0 50])


%%% second diagnostic checks least squares fit for camera center
%%% if correct, vector from ellipse center to camera center should be
%%% orthogonal to ellipse axis vector (cos w, sin w) (bottom of p. 27)
%%% check this by taking dot-product of the thwo vectors, and plotting vs ellipse-camera dist
%%% it should be zero for sufficient ellipse-camera distances (fails for shor distances due to low ellipticity)
%%% note - we plot this separately for the points that were included in fit above, versus all points.
subplot(2,2,4)
delta = (CamCent-EllipseParams(:,1:2)');  %%% camera-ellipse vector

plot(vecnorm(delta(:,usegood),2,1),((delta(1,usegood)'.*cos(EllipseParams(usegood,5)))+(delta(2,usegood)'.*sin(EllipseParams(usegood,5))))./vecnorm(delta(:,usegood),2,1)','.')
hold on
plot(vecnorm(delta(:,list),2,1),((delta(1,list)'.*cos(EllipseParams(list,5)))+(delta(2,list)'.*sin(EllipseParams(list,5))))./vecnorm(delta(:,list),2,1)','r.')
title('Camera Center Calibration')
ylabel('[PC-EC].[cosw;sinw]');
xlabel('abs(PC-EC)');
legend('all points','list points')
if savePDF
    set(gcf, 'PaperPositionMode', 'auto');   print('-dpsc',psfilename,'-append'); close(gcf)
end



%% Theta and Phi Final values

%%% occasionally get complex values if calibration is way off
%%% (this step may not be a good idea, since it masks those, but keeps code from crashing)
newtheta = (real(thetad));
newphi = (real(phi));

%%% plot smoothed traces of theta/phi
figure; subplot(121);
plot(diff(movmean(newtheta,10)),'LineWidth',2); %hold on; plot(diff(newlongang),'LineWidth',1)
title('Change in smoothed theta')
subplot(122);
plot(diff(movmean(newphi,10)),'LineWidth',2);  % hold on; plot(diff(newlongang),'LineWidth',1)
title('Change in smoothed phi')







