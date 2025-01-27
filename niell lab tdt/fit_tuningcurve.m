function [theta_pref cvOSI cvDSI A1 A2 w B orth_rate yfit ] = fit_tuningcurve(R, theta_ind);

% temp_baseline =min(R); 
% R=R-temp_baseline;  %%%need to keep things positive for fourier analysis

theta_ind
if abs(min(R))>max(R)  %%% check for units that are inhibited more than excited
    R = -R;
    inverted=1;
else
    inverted=0;
end

if min(R)<0
    R = R-min(R);
end


i=sqrt(-1);
mu = (sum(R.*exp(2*i*theta_ind)))/sum(abs(R));
dsi = (sum(R.*exp(i*theta_ind)))/sum(abs(R));



%osi = sum(R.*exp(sqrt(-1)*theta'))/sum(R)
%osi = sum(R.*exp(sqrt(-1)*theta))/sum(R);

if isnan(mu);
    mu=0;
end

if isnan(dsi);
    dsi=0;
end
cvOSI = abs(mu)
cvDSI = abs(dsi)
theta_pref = mod(angle(mu)/2,2*pi)
 
  


%R= R+temp_baseline;

if max(theta_ind)<pi;
    theta_ind(1,size(theta_ind,2)+1 : 2*size(theta_ind,2)) = theta_ind+pi;
    R(1,size(R,2)+1 : 2*size(R,2)) = R;
end

theta_ind(1,size(theta_ind,2)+1) = 2*pi;
R(1,size(R,2)+1)=R(1,1);

p0(1) = interp1(theta_ind,R,theta_pref,'spline') ;  %%% A1
p0(2) =  interp1(theta_ind,R,mod(theta_pref+pi,2*pi),'spline');  %%%A2
p0(3) =pi/8; %%%width
p0(4) = interp1(theta_ind,R,mod(theta_pref+pi/2,2*pi),'spline') ; %%%baseline
p0(1)= p0(1)-p0(4);
p0(2) = p0(2)-p0(4);

p0

dtheta = 2*pi/32;
clear x
x(1,:) = 0:dtheta:2*pi-dtheta % angle/stim at which you get the response
y = interp1(theta_ind,R,x(1,:),'spline'); %firing response magnitude
x(2,:) = theta_pref;

p = nlinfit(x,y,@wrapped_gaussian,p0);

xfit(1,:)=theta_ind;
xfit(2,:)=theta_pref;
yfit = wrapped_gaussian(p,xfit);

% figure
% plot(x(1,1:31)*180/pi,y(1:31))
% hold on
% plot(x(1,1:2:31)*180/pi,y(1:2:31),'o');
% plot(x(1,:)*180/pi,wrapped_gaussian(p,x),'g');


A1=max(p(1),p(2));
A2=min(p(1),p(2));
if abs(p(2))>abs(p(1))
    theta_pref = theta_pref+pi;
    A1=p(2);
    A2=p(1);
else
    A1=p(1);
    A2=p(2);
end

w=p(3);
B = p(4);

if inverted
    A1=-A1;
    A2=-A2;
    B=-B;
    yfit = -yfit;
end
if theta_pref<0
    theta_pref=theta_pref+2*pi;
elseif theta_pref>2*pi
    theta_pref = theta_pref-2*pi;
end

orth_rate = interp1(theta_ind,yfit,mod(theta_pref+pi/2,2*pi));
%null_rate = interp1(theta_ind,yfit,mod(theta_pref+pi,2*pi));


A1
A2 
B

