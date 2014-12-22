function [tsamp vsmooth] = getBlockVelocity_general(tank,block,version);
%%% reads in mouse samples and converts to real-world movement
%%% written by cmn, based on TankGetMouseEvents by sr
%%% updated to use both versions of mouse sender software 
%%% ie. original for TDT (0), and update for contstim (1) which packs data

if version==0
C=18000;  %%% circumference in ticks
else
    C=9800;  %%% circumference in ticks
    %%% C=7000;
end

scale_factor = C/(pi*8*2.54);  %%% ball is 8 inches, so this is ticks/cm

%%% read data
if version ==0 %%% old version, asynchronous differential position
    [dTS, dData] = TankGetMouseEvents(tank,block);

    fprintf('Processing...\n');

    %%% parse data into x,y, and times for each mouse
    Mouse1_Idx = find(dData(1,:)==0);
    Mouse2_Idx = find(dData(1,:)==1);
    M1_dX = dData(3,Mouse1_Idx);
    M1_dY = dData(4,Mouse1_Idx);

    M1_T = dTS(Mouse1_Idx);
    M2_T = dTS(Mouse2_Idx);

    Mouse2_Idx = find(dData(1,:)==1);
    M2_dX = dData(3,Mouse2_Idx);
    M2_dY = dData(4,Mouse2_Idx);
    
    M1_X = cumsum(M1_dX);
    M1_Y = cumsum(M1_dY);

    M2_X = cumsum(M2_dX);
    M2_Y = cumsum(M2_dY);
else %%% new 30hz cumulative position
    [dTS, dData] = TankGetUDP_SyncEvents(tank,block);

    fprintf('Processing...\n');
    dData
    
    Mouse1_Idx = find(dData(1,:)==2); % all events on channel 2
    dXY = dData(2,Mouse1_Idx);
    [M1_X  M1_Y] = GetMouseTraj(dXY);
    M1_T = dTS(Mouse1_Idx);
    
     Mouse2_Idx = find(dData(1,:)==3); % all events on channel 3
    dXY = dData(2,Mouse2_Idx);
    [M2_X  M2_Y] = GetMouseTraj(dXY);
    M2_T = dTS(Mouse2_Idx);
end
    

figure
plot(M1_T)


figure
plot(M2_T)



%%% resample data into even time intervals
dt = 0.1;
Tmax = max(max(M1_T),max(M2_T))
Tmin =max(min(M1_T),min(M2_T))

tsamp = Tmin:dt:Tmax;
M1_Xsamp = interp1(M1_T,M1_X,tsamp);
M2_Xsamp = interp1(M2_T,M2_X,tsamp);
M1_Ysamp = interp1(M1_T,M1_Y,tsamp);
M2_Ysamp = interp1(M2_T,M2_Y,tsamp);

max(M1_T)
max(M2_T)

M1_dXsamp = diff(M1_Xsamp);
M2_dXsamp = diff(M2_Xsamp);
M1_dYsamp = diff(M1_Ysamp);
M2_dYsamp = diff(M2_Ysamp);

%%% calculate real-world position by integrating head orientation
%%% and forward/orthogonal movement

theta =zeros(size(M1_dYsamp));
x=zeros(size(M1_dYsamp)); y=zeros(size(M1_dYsamp));

figure
plot(M1_dXsamp,M2_dXsamp,'o')

length(M2_dXsamp)
find(isnan(M2_dXsamp));

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
plot(x/scale_factor,y/scale_factor,'.')
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
total_distance = sum(v)/(100*scale_factor*dt);
sprintf('ran %f meters in %f mins',total_distance,max(tsamp)/60)


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

