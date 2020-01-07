fn = 15

load('J462a_originalACC.mat')
d1 = Data(fn);
load('J462a_deInterACC.mat')
d2 = Data(fn);

n = min(length(d1.Rtheta),length(d2.Rtheta));

figure
plot(d1.Rtheta(1:n),d2.Rtheta(1:n),'.')

figure
plot(diff(d1.Rtheta(1:n)),diff(d2.Rtheta(1:n)),'.')

figure
plot(d1.Rtheta(1:n),'b'); hold on; plot(d2.Rtheta(1:n),'r')


figure
plot(nanxcorr(d1.dxRTheta,diff(d1.theta),30,'coeff'))
hold on
plot(nanxcorr(d2.dxRTheta,diff(d2.theta),30,'coeff'))
title('dEyeTheta vs dHeadTheta');

figure
plot(d1.dxRTheta,diff(d1.theta)*180/pi,'b.');
axis equal
hold on
plot(d2.dxRTheta(1:n-1),diff(d1.theta)*180/pi,'r.');
axis equal
xlabel('dEye theta'); ylabel('dHead theta');

figure
plot(nanxcorr(d1.dxRTheta,d1.dxLTheta,30,'coeff'));
hold on
plot(nanxcorr(d2.dxRTheta,d2.dxLTheta,30,'coeff'))
title('dRtheta vs dLtheta');


