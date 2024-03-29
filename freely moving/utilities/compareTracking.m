fn =1

load('J462c_originalACC_new.mat')
d1 = Data(fn);
load('J462c_deInterACC_NN.mat')
d2 = Data(fn);

n = min(length(d1.Rtheta),length(d2.Rtheta));

figure
plot(d1.Rtheta(1:n),d2.Rtheta(1:n),'.')

figure
plot(diff(d1.Rtheta(1:n)),diff(d2.Rtheta(1:n)),'.')

figure
plot(d1.Rtheta(1:n),'b'); hold on; plot(d2.Rtheta(1:n),'r')


figure
plot(nanxcorr(d1.dxRTheta,d1.dtheta(1:end-1),30,'coeff'))
hold on
plot(nanxcorr(d2.dxRTheta,d2.dtheta(1:end-1),30,'coeff'))
title('dEyeTheta vs dHeadTheta');


figure
plot(nanxcorr(d1.dxRTheta,d1.accShift(1:end-1,6),30,'coeff'))
hold on
plot(nanxcorr(d2.dxRTheta,d2.accShift(1:end-1,6),30,'coeff'))
title('dEyeTheta vs head acc');

g3 = 0.5*(d1.accShift(1:end-1,6) + d1.accShift(2:end,6));
figure
plot(nanxcorr(d1.dxRTheta,g3,30,'coeff'))

figure
plot(d1.dxRTheta,d1.dtheta(1:end-1)*180/pi,'b.');
axis equal
hold on
plot(d2.dxRTheta(1:n-1),d2.dtheta(1:end-1)*180/pi,'r.');
axis equal
xlabel('dEye theta'); ylabel('dHead theta');

figure
plot(d1.dxRTheta,d1.accShift(1:end-1,6),'b.');
axis equal
hold on
plot(d2.dxRTheta(1:n-1),d2.accShift(1:end-1,6),'r.');
axis equal
xlabel('dEye theta'); ylabel('dHead theta');


figure
plot(nanxcorr(d1.dxRTheta,d1.dxLTheta,30,'coeff'));
hold on
plot(nanxcorr(d2.dxRTheta,d2.dxLTheta,30,'coeff'))
title('dRtheta vs dLtheta');

sprintf('%d nans raw vs %d nans deInter',sum(isnan(d1.Rtheta)),sum(isnan(d2.Rtheta)))

figure
plot(nanxcorr(d1.accShift(:,6),d1.dtheta,30,'coeff'))

figure
plot(d1.accShift(:,6));
hold on
plot(d1.dtheta*180/pi)
