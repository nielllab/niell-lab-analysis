function [tsamp vsmooth groomstate] = getBlockVelocity_groom(tank,block);
%%% reads in mouse samples and converts to real-world movement
%%% written by cmn, based on TankGetMouseEvents by sr


%%% read data
[dTS, dData] = TankGetMouseEvents(tank,block);

C=18000;  %%% circumference in ticks
scale_factor = C/(pi*8*2.54);  %%% ball is 8 inches, so this is ticks/cm


fprintf('Processing...\n');

%%% parse data into x,y, and times for each mouse
Mouse1_Idx = find(dData(1,:)==0);
Mouse2_Idx = find(dData(1,:)==1);
M1_dX = dData(3,Mouse1_Idx);
M1_dY = dData(4,Mouse1_Idx);

M2_clicks= find((dData(1,:)==1)&dData(5,:)~=0);

clicktimes = dTS(M2_clicks)
clickevent = dData(5,M2_clicks)
click_on=4;
click_off=1;

M1_X = cumsum(M1_dX);
M1_Y = cumsum(M1_dY);

M1_T = dTS(Mouse1_Idx);
M2_T = dTS(Mouse2_Idx);

Mouse2_Idx = find(dData(1,:)==1);
M2_dX = dData(3,Mouse2_Idx);
M2_dY = dData(4,Mouse2_Idx);

M2_X = cumsum(M2_dX);
M2_Y = cumsum(M2_dY);

%%sw= input(' optical mice = horizontal (0) or vertical (1) : ');
display('using vertical mice')
sw = 1;

if sw ==0;
    m = M1_X;
    M1_X= M1_Y;
    M1_Y = m;
    
    m = M2_X;
    M2_X=M2_Y;
    M2_Y=m;
end

%%% resample data into even time intervals
dt = 0.1;
Tmax = max(max(M1_T),max(M2_T))
Tmin =max(min(M1_T),min(M2_T))

tsamp = Tmin:dt:Tmax;
if M2_T(1) == M2_T(2);
    M2_T = M2_T(2:end); %%% somehow first samp can get duplicated
    M2_X = M2_X(2:end);
    M2_Y = M2_Y(2:end);
end

if M1_T(1) == M1_T(2);
    M1_T = M1_T(2:end); %%% somehow first samp can get duplicated
    M1_X = M1_X(2:end);
    M1_Y = M1_Y(2:end);
end

M1_Xsamp = interp1(M1_T,M1_X,tsamp);
M2_Xsamp = interp1(M2_T,M2_X,tsamp);
M1_Ysamp = interp1(M1_T,M1_Y,tsamp);
M2_Ysamp = interp1(M2_T,M2_Y,tsamp);

M1_dXsamp = diff(M1_Xsamp);
M2_dXsamp = diff(M2_Xsamp);
M1_dYsamp = diff(M1_Ysamp);
M2_dYsamp = diff(M2_Ysamp);

%%% calculate real-world position by integrating head orientation
%%% and forward/orthogonal movement

theta =zeros(size(M1_dYsamp));
x=zeros(size(M1_dYsamp)); y=zeros(size(M1_dYsamp));

dTheta=-2*pi*0.5*(M1_dXsamp+M2_dXsamp)/C;
M2_dYsamp =-1*M2_dYsamp;  %%% to keep axes following righthand rule

theta(1)=0+dTheta(1);
x(1)=0+M1_dYsamp(1);
y(1) = 0+M2_dYsamp(1);
for t=2:length(M1_dYsamp);
    x(t) = x(t-1) + M1_dYsamp(t-1)*cos(theta(t-1)) - M2_dYsamp(t-1)*sin(theta(t-1));
    y(t) = y(t-1) + M1_dYsamp(t-1)*sin(theta(t-1)) + M2_dYsamp(t-1)*cos(theta(t-1));
    theta(t) = theta(t-1)+dTheta(t);
end
  

%%% plot trajectory
figure
plot(x/(100*scale_factor*dt),y/(100*scale_factor*dt),'.')
axis equal

%%% plot velocity as a function of time
v = sqrt(diff(x).^2 + diff(y).^2)/(scale_factor*dt);
vsmooth = conv(v,ones(1,11))/10;
vsmooth=vsmooth(6:length(vsmooth)-4);
tsamp = tsamp(1:length(tsamp)-1)+dt/2;
figure
% plot(tsamp,v)
% hold on
plot(tsamp,vsmooth,'g')
groomstate = zeros(1,length(tsamp));

for i = 1:length(clicktimes)
    if clickevent(i)==click_on;
        groomstate(tsamp>clicktimes(i))=1;
    else
      groomstate(tsamp>clicktimes(i))=0;
    end
end
figure
plot(tsamp,vsmooth,'g')
hold on
plot(tsamp,groomstate*50);
%axis([0  max(tsamp) -0.2 1.2])
  %%% because velocity is at center of time samples

%%% plot ball rotational velocities
% figure
% plot(M1_dXsamp);
% hold on
% plot(M1_dYsamp,'g');
% 
% figure
% plot(M2_dXsamp);
% hold on
% plot(M2_dYsamp,'g');

