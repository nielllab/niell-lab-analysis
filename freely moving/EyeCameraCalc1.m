%%%%% Eye Camera Calculations
function [newtheta,newphi,EllipseParams,ExtraParams] = EyeCameraCalc1(numFrames,Pointsx,Pointsy,psfilename)

% Inputs:
%   Vid1 - 3D grayscale array of video frames
%   Pointsx - The X component of Deeplabcut tracking
%   Pointsy - The Y component of Deeplabcut tracking

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

if exist('psfilename','var')
    savePDF=1;
end



%% Calculate Ellipse Fits
% figure;
% axis tight manual
% ax = gca;
% % ax.NextPlot = 'replaceChildren';
% axis([1 size(Vid1,1) 1 size(Vid1,2)])
EllipseParams = zeros(numFrames,7);
ExtraParams = zeros(numFrames,6);

parfor v=1:(numFrames)
    if ~any((Pointsx(v,:)<10 | Pointsy(v,:)<10)| isnan(Pointsx(v,:)))
        e_t = fit_ellipse2(Pointsx(v,:),Pointsy(v,:));
        if isempty(e_t.status)==1 &&  ((e_t.short_axis/2)/(e_t.long_axis/2) >.5)
            EllipseParams(v,:)=[e_t.X0_in, e_t.Y0_in, e_t.long_axis/2, e_t.short_axis/2,  e_t.angleToX*pi/180, e_t.angleFromX*pi/180, e_t.phi];
            ExtraParams(v,:) = [e_t.X0, e_t.Y0, e_t.a, e_t.b, e_t.cos_phi, e_t.sin_phi];
        end
    end
end
fprintf('done \n')
close all;

efitb = find(EllipseParams(:,1)==0);
EllipseParams(efitb,:)=NaN;
ExtraParams(efitb,:) = NaN;
EllipseParams = fillmissing(EllipseParams,'linear',1);
ExtraParams = fillmissing(ExtraParams,'linear',1);


%%  Calc Camera Center
R = linspace(0,2*pi,100);
list = find(EllipseParams(:,4)./EllipseParams(:,3)<.95); %randi([1 size(EllipseParams,1)],50);%  1:size(EllipseParams,1); %
A = [cos(EllipseParams(list,5)),sin(EllipseParams(list,5))];
b=diag(A*EllipseParams(:,1:2)');
CamCent=(A'*A)\A'*b;
%
% scale = nansum(sqrt(1-(EllipseParams(:,4)./EllipseParams(:,3)).^2)'*vecnorm([EllipseParams(:,1)';EllipseParams(:,2)']-CamCent,2,1)')/...
%     nansum(1-(EllipseParams(:,4)./EllipseParams(:,3)).^2);
Ellipticity = EllipseParams(list,4)./EllipseParams(list,3);
scale = nansum(sqrt(1-(Ellipticity).^2).*vecnorm(EllipseParams(list,1:2)'-CamCent,2,1)')./nansum(1-(Ellipticity).^2);

theta = asin((EllipseParams(:,1)-CamCent(1))*1/scale);
thetad =asind((EllipseParams(:,1)-CamCent(1))*1/scale);
phi = asind((EllipseParams(:,2)-CamCent(2))./cos(theta)*1/scale);

%%
i=50;
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

figure;
%imagesc(Vid1(:,:,i)); colormap gray;
axis equal off;
plot(q_star(1,:),q_star(2,:),'g','LineWidth',2); scatter(EllipseParams(i,1),EllipseParams(i,2),100,'og');hold on;
scatter(Pointsx(i,:),Pointsy(i,:),100,'.m');
scatter(CamCent(1),CamCent(2),100,'or'); scatter(qcirc2(1,:),qcirc2(2,:),50,'.r');
scatter(newCent(1),newCent(2),100,'xb')
scatter(PointsRot(1,:),PointsRot(2,:),100,'.b');

% scatter(EyeShift(i,2),EyeShift(i,3),100,'.y')

if savePDF
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end
close(gcf)


%% Check Calibration
figure;
plot(vecnorm(EllipseParams(:,1:2)' - CamCent),scale*sqrt(1-(EllipseParams(:,4)./EllipseParams(:,3)).^2),'.')
axis equal;hold on;
plot(linspace(0,250),linspace(0,250),'r')
title('Scale Factor Calibration')
axis([0 100 0 100])
if savePDF
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end
close(gcf)

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
if savePDF
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end
close(gcf)


%% Theta and Phi

newtheta = (real(thetad));
newphi = (real(phi));

figure; subplot(121);
plot(diff(movmean(newtheta,10)),'LineWidth',2); %hold on; plot(diff(newlongang),'LineWidth',1)
title('Change in smoothed theta')
subplot(122);
plot(diff(movmean(newphi,10)),'LineWidth',2);  % hold on; plot(diff(newlongang),'LineWidth',1)
title('Change in smoothed phi')
if savePDF
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
end
close(gcf)




