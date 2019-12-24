load('J462a_cmnTest.mat')
d1 = Data(3);
load('J462a_cmnTest_deInter.mat')
d2 = Data(3);

figure
plot(d1.Rtheta,d2.Rtheta,'.')

figure
plot(d1.Rtheta); hold on; plot(d2.Rtheta)

figure
plot(nanxcorr(d1.dxRTheta,diff(d1.theta),30,'coeff'))
hold on
plot(nanxcorr(d2.dxRTheta,diff(d2.theta),30,'coeff'))

figure
plot(d1.dxRTheta,diff(d1.theta)*180/pi,'.');
axis equal

figure
plot(d2.dxRTheta,diff(d1.theta)*180/pi,'.');
axis equal

