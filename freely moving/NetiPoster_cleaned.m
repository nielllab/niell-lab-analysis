[thetaL,phiL,EllipseParamsL,~] = EyeCameraCalc1(EyeVidL,PointsxL,PointsyL);

%Right Eye
[thetaR,phiR,EllipseParamsR,~] = EyeCameraCalc1(EyeVid,Pointsx,Pointsy);


%%

app1=30:60
app2 = (6*30):(30*30)
app3 = 30*(120+42):30*(120+49)
app4=30*(120+51):30*(120+54)
app5= 30*(180+58):30*(180+1)
app6= 30*(180+5):30*(180+10)
app7=30*(300+33):30*(300+39)
app8 = 30*(420+30):30*(420+46)
app9= 30*(540+34):30*(540+36)

apps=[app1, app2, app3, app4, app5, app5, app6, app7, app8, app9]

7:30 7:46


mnRadR = (EllipseParamsR(:,3))./(mean(EllipseParamsR(:,3)));
mnRadL = (EllipseParamsL(:,3))./(mean(EllipseParamsL(:,3)));

radR = medfilt1(mnRadR,20); radL = medfilt1(mnRadL,20);
radR=radR(1:end-1)
figure;
plot(radR,'r'); hold on; plot(radL,'b')
sec = length(radR)./30; min = sec./60
set(gca,'xtick',([0:1800:length(radR)]))% tick every 1 min
set(gca,'xticklabels',({'0','1','2','3','4','5','6','7','8','9','10'}),'FontSize',10)
xlabel('Time(min)'); ylabel('normalized pupil radius')

figure;plot(radR(1:40:length(radR)),radL(1:40:length(radL)),'ko','Markersize', 4);axis square; xlabel('Right Eye Normalized Pupil Radius');
ylabel('Left Eye Normalized Pupil Radius'); xlim([.7 1.5]); ylim([.7 1.5])
hold on; plot(radR(app1),radL(app1),'og'); plot(radR(app2),radL(app2),'go');
plot(radR(app3),radL(app3),'og'); plot(radR(app4),radL(app4),'go');
plot(radR(app5),radL(app5),'og'); plot(radR(app6),radL(app6),'go');
plot(radR(app7),radL(app7),'og'); plot(radR(app8),radL(app8),'go');
plot(radR(app9),radL(app9),'go');


%%
figure;
plot(radR,'r'); hold on; plot(radL,'b')
sec = length(radR)./30; min = sec./60
set(gca,'xtick',([0:1800:length(radR)]))% tick every 1 min
set(gca,'xticklabels',({'0','1','2','3','4','5','6','7','8','9','10'}),'FontSize',10)
xlabel('Time(min)'); ylabel('normalized pupil radius')

plot(app1,radR(app1),'g'); plot(app2,radR(app2),'g');plot(app3,radR(app3),'g');
plot(app4,radR(app4),'g'); plot(app5,radR(app5),'g'); plot(app6,radR(app6),'g');plot(app7,radR(app7),'g');
plot(app8,radR(app8),'g');plot(app9,radR(app9),'g');

plot(app1,radL(app1),'g'); plot(app2,radL(app2),'g');plot(app3,radL(app3),'g');
plot(app4,radL(app4),'g'); plot(app5,radL(app5),'g'); plot(app6,radL(app6),'g');plot(app7,radL(app7),'g');
plot(app8,radL(app8),'g');plot(app9,radL(app9),'g');


%%

app1=30:10:60
app2 = (6*30):10:(30*30)
app3 = 30*(120+42):10:30*(120+49)
app4=30*(120+51):10:30*(120+54)
app5= 30*(180+58):10:30*(180+1)
app6= 30*(180+5):10:30*(180+10)
app7=30*(300+33):10:30*(300+39)
app8 = 30*(420+30):10:30*(420+46)
app9= 30*(540+34):10:30*(540+36)

thetaRfilt=thetaR(1:10:length(thetaR));thetaLfilt=thetaL(1:10:length(thetaL))
phiRfilt=phiR(1:10:length(phiR));phiLfilt=phiL(1:10:length(phiL))
 %thetaR=thetaR(1:end-1)
%phiR=phiR(1:end-1)

figure;plot(thetaRfilt, thetaLfilt,'ok'); hold on; axis square; %xlim([-80 80]); ylim([-80 80]);
plot(thetaR(app1),thetaL(app1),'go'); plot(thetaR(app2),thetaL(app2),'go');
plot(thetaR(app3),thetaL(app3),'go'); plot(thetaR(app4),thetaL(app4),'go');
plot(thetaR(app5),thetaL(app5),'go'); plot(thetaR(app6),thetaL(app6),'go');
plot(thetaR(app7),thetaL(app7),'go'); plot(thetaR(app8),thetaL(app8),'go');
plot(thetaR(app9),thetaL(app9),'go');
xlabel ('R theta'); ylabel('L theta');


figure;plot(phiRfilt, phiLfilt,'ok'); hold on; axis square; %xlim([-80 80]); ylim([-80 80]);
plot(phiR(app1),phiL(app1),'go'); plot(phiR(app2),phiL(app2),'go');
plot(phiR(app3),phiL(app3),'go'); plot(phiR(app4),phiL(app4),'go');
plot(phiR(app5),phiL(app5),'go'); plot(phiR(app6),phiL(app6),'go');
plot(phiR(app7),phiL(app7),'go'); plot(phiR(app8),phiL(app8),'go');
plot(phiR(app9),phiL(app9),'go');
xlabel ('R phi'); ylabel('L phi'); xlim([-80 80]);ylim([-80 80]);

%%



app1=30:60
app2 = (6*30):(30*30)
app3 = 30*(120+42):30*(120+49)
app4=30*(120+51):30*(120+54)
app5= 30*(180+58):30*(180+1)
app6= 30*(180+5):30*(180+10)
app7=30*(300+33):30*(300+39)
app8 = 30*(420+30):30*(420+46)
app9= 30*(540+34):30*(540+36)

apps=[app1, app2, app3, app4, app5, app5, app6, app7, app8, app9]

du = find(thetaR<0 & thetaL>0)
uu = find(thetaR>0 & thetaL>0)
dd = find(thetaR<0 & thetaL<0)
ud = find(thetaR>0 & thetaL<0)


binRange = -80:5:80
hcx = histcounts(thetaR,[binRange Inf]);
hcy = histcounts(thetaL,[binRange Inf]);
hcxa = histcounts(thetaR(apps),[binRange Inf]);
hcya = histcounts(thetaL(apps),[binRange Inf]);



figure
b1=bar([hcx./length(thetaR);hcy./length(thetaL)]'); hold on


b2=bar([hcxa./length(apps);hcya./length(apps)]','r')


figure;plot(hcx./length(thetaR)+hcy./length(thetaL),'b'); hold on;
plot(hcxa./length(apps)+ hcya./length(apps),'r')




%%
velThetaL=diff(thetaL); velThetaR=diff(thetaR);
%velThetaL= velThetaL(1:30:length(velThetaL)); velThetaR = velThetaR(1:30:length(velThetaR))


app1=30:20:60
app2 = (6*30):20:(30*30)
app3 = 30*(120+42):20:30*(120+49)
app4=30*(120+51):20:30*(120+54)
app5= 30*(180+58):20:30*(180+1)
app6= 30*(180+5):20:30*(180+10)
app7=30*(300+33):20:30*(300+39)
app8 = 30*(420+30):20:30*(420+46)
app9= 30*(540+34):20:30*(540+36)

figure;
subplot(1,2,1)
plot(velThetaR,velThetaL,'o'); hold on
plot(velThetaR(app1), velThetaL(app1),'og'); axis equal; axis square
plot(velThetaR(app2), velThetaL(app2),'og');
plot(velThetaR(app3), velThetaL(app3),'og');
plot(velThetaR(app4), velThetaL(app4),'og');
plot(velThetaR(app5), velThetaL(app5),'og');
plot(velThetaR(app6), velThetaL(app6),'og');
plot(velThetaR(app7), velThetaL(app7),'og');
plot(velThetaR(app8), velThetaL(app8),'og');
plot(velThetaR(app9), velThetaL(app9),'og');

%ylim([-9 9]); xlim([-15 15]);

xlabel('Theta R Eye');ylabel('Theta L Eye');

velPhiL=diff(phiL); velPhiR=diff(phiR);

subplot(1,2,2)
plot(velPhiR,velPhiL,'o'); hold on
plot(velPhiR(app1), velPhiL(app1),'og'); axis square
plot(velPhiR(app2), velPhiL(app2),'og');
plot(velPhiR(app3), velPhiL(app3),'og');
plot(velPhiR(app4), velPhiL(app4),'og');
plot(velPhiR(app5), velPhiL(app5),'og');
plot(velPhiR(app6), velPhiL(app6),'og');
plot(velPhiR(app7), velPhiL(app7),'og');
plot(velPhiR(app8), velPhiL(app8),'og');
plot(velPhiR(app9), velPhiL(app9),'og');

%ylim([-11 11]); xlim([-25 25]);


xlabel('Phi R Eye');ylabel('Phi L Eye');

%%

app1=30:60
app2 = (6*30):(30*30)
app3 = 30*(120+42):30*(120+49)
app4=30*(120+51):30*(120+54)
app5= 30*(180+58):30*(180+1)
app6= 30*(180+5):30*(180+10)
app7=30*(300+33):30*(300+39)
app8 = 30*(420+30):30*(420+46)
app9= 30*(540+34):30*(540+36)

headangle= radtodeg(atan2(PointsxT(:,2)-PointsxT(:,4),PointsyT(:,2)-PointsyT(:,4)))
HA= medfilt1(diff(headangle),10)

figure;subplot(2,1,1);plot(headangle);hold on; ylabel('head angle (degrees)')
plot(app1,(headangle(app1)),'g')
plot(app2,(headangle(app2)),'g')
plot(app3,(headangle(app3)),'g')
plot(app4,(headangle(app4)),'g')
plot(app5,(headangle(app5)),'g')
plot(app6,(headangle(app6)),'g')
plot(app7,(headangle(app7)),'g')
plot(app8,(headangle(app8)),'g')


figure;plot(HA); title('head angle velocity'); hold on
plot(app1,(HA(app1)),'g')
plot(app2,(HA(app2)),'g')
plot(app3,(HA(app3)),'g')
plot(app4,(HA(app4)),'g')
plot(app5,(HA(app5)),'g')
plot(app6,(HA(app6)),'g')
plot(app7,(HA(app7)),'g')
plot(app8,(HA(app8)),'g')
plot(app9,(HA(app9)),'g')

set(gca,'xtick',([0:1800:length(HA)]))% tick every 1 min
set(gca,'xticklabels',({'0','1','2','3','4','5','6','7','8','9','10'}),'FontSize',10)
xlabel('Time (min)'); ylabel('head angle velocity')


%%
length(apps)/length(radR)

figure; plot(diff(PointsxT(:,6))./mean(diff(PointsxT(:,6))))
hold on; plot(radR)

%%
figure
plot(PointsxT(app1,6),PointsyT(app1,6),'k'); hold on
plot(PointsxT(app1,7),PointsyT(app1,7),'g')


%%
HAmn = HA./mean(HA)

[corrL,lagsL] = xcorr(HA, thetaL);

figure;plot(lagsL, corrL); xlim([-30 30]);

figure;plot(xcorr(HA, thetaR))

figure
plot(HA(1:18000),thetaR(11:18010),'.');

figure
plot(HA(1:18000),thetaR(11:18010),'.');

%%HAmn = HA./mean(HA)

[corrR,lagsR] = xcorr(HA(1:length(velThetaR)), velThetaR,'normalize'); 


figure;plot(lagsR, corrR); xlim([-30 30]);
figure;plot(lagsR, corrR); 

%figure;plot(xcorr(HA, velThetaR))

% figure
% plot(HA(1:18000),velThetaR(11:18010),'.');
% 
% figure
% plot(HA(1:18000),velThetaR(11:18010),'.');


HAmn = HA./mean(HA)

[corrL,lagsL] = xcorr(HA, velThetaL);

figure;plot(lagsL, corrL); xlim([-30 30]);
figure;plot(lagsL, corrL); xlim([-30 30]);

%figure;plot(xcorr(HA, velThetaL))
% 
% figure
% plot(HA(1:18000),velThetaL(11:18010),'.');
% 
% figure
% plot(HA(1:18000),velThetaL(11:18010),'.');

%%HAmn = HA./mean(HA)
HA = HA(1:length(velThetaR))


ap1 = 30:15:270
ap2 = 2700:15:3690
ap3 = 6690:15:7230
ap4 = 10470:15:10550
ap5 = 13440:15:13560 %7:28 -7:32
ap6 = 14220:15:14430% 7:54 - 8:01
ap7 = 14820:15:14910  %8:14-8:17
ap8 = 15210:15:15600;

apps= [ap1, ap2, ap3, ap4, ap5, ap6, ap7, ap8]

[corrR,lagsR] = xcorr(HA(1:length(velThetaR)), velThetaR,'normalize'); 
corrcoef(HA(1:length(velThetaR)), velThetaR)

figure; subplot(1,2,1);plot(lagsR, corrR,'r');xlim([-30 30]);title('right eye');
set(gca,'xtick',([-30:15:30]))% tick every 1 min
set(gca,'xticklabels',({'- 1','-.5','0','.5','1'}),'FontSize',10)
subplot(1,2,2);plot(lagsL, corrL); xlim([-30 30]); title('left eye')
set(gca,'xtick',([-30:15:30]))% tick every 500 ms
set(gca,'xticklabels',({'- 1','-.5','0','.5','1'}),'FontSize',10)


%%

figure;
subplot(1,2,1)
scatter(HA, velThetaR); hold on;scatter(HA(apps), velThetaR(apps),'g'); title('right eye'); axis square
subplot(1,2,2)
scatter(HA, velThetaL);  hold on;scatter(HA(apps), velThetaL(apps),'g'); title('left eye'); axis square
 
figure;
subplot(1,2,1)
scatter(velThetaR,HA); hold on;scatter(velThetaR(apps),HA(apps),'g'); title('right eye'); axis square
subplot(1,2,2)
scatter(velThetaL,HA);  hold on;scatter(velThetaL(apps),HA(apps),'g'); title('left eye'); axis square
%%
% subplot(2,2,3);plot(lagsR, corrR,'r'); title('right eye');
% subplot(2,2,4);plot(lagsL, corrL); title('left eye'); 

% set(gca,'xtick',([-30:15:30]))% tick every 1 min
% set(gca,'xticklabels',({'0','1','2','3','4','5','6','7'}),'FontSize',10)


%figure;plot(xcorr(HA, velThetaR))

% figure
% plot(HA(1:18000),velThetaR(11:18010),'.');
% 
% figure
% plot(HA(1:18000),velThetaR(11:18010),'.');


HAmn = HA./mean(HA)

[corrL,lagsL] = xcorr(HA(1:length(velThetaL)), velThetaL,'normalize'); 
corrcoef(HA(1:length(velThetaL)), velThetaL)


figure;plot(lagsL, corrL); xlim([-30 30]);
figure;plot(lagsL, corrL);


%%

[corrR,lagsR] = xcorr(HA(1:length(velPhiR)), velPhiR,'normalize'); 
[corrL,lagsL] = xcorr(HA(1:length(velPhiL)), velPhiL,'normalize'); 

figure; subplot(1,2,1);plot(lagsR, corrR,'r');xlim([-30 30]);title('right eye');
set(gca,'xtick',([-30:15:30]))% tick every 1 min
set(gca,'xticklabels',({'- 1','-.5','0','.5','1'}),'FontSize',10)
subplot(1,2,2);plot(lagsL, corrL); xlim([-30 30]); title('left eye')
set(gca,'xtick',([-30:15:30]))% tick every 500 ms
set(gca,'xticklabels',({'- 1','-.5','0','.5','1'}),'FontSize',10)
%%


%figure;plot(xcorr(HA, velThetaL))

% figure
% plot(HA(1:18000),velThetaL(11:18010),'.');
% 
% figure
% plot(HA(1:18000),velThetaL(11:18010),'.');

%% thetaRfilt=thetaR(1:30:length(thetaR));thetaLfilt=thetaL(1:30:length(thetaL))
phiRfilt=phiR(1:30:length(phiR));phiLfilt=phiL(1:30:length(phiL))

figure; subplot(1,2,1)
plot(thetaRfilt, phiRfilt,'ok'); hold on; axis square; %xlim([-80 80]); ylim([-80 80]);
plot(thetaR(ap1),phiR(ap1),'go'); plot(thetaR(ap2),phiR(ap2),'go');
plot(thetaR(ap3),phiR(ap3),'go'); plot(thetaR(ap4),phiR(ap4),'go');
plot(thetaR(ap5),phiR(ap5),'go'); plot(thetaR(ap6),phiR(ap6),'go');
plot(thetaR(ap7),phiR(ap7),'go'); plot(thetaR(ap8),phiR(ap8),'go');
xlabel ('R theta'); ylabel('R phi');


subplot(1,2,2); plot(thetaLfilt, phiLfilt,'ok'); hold on; axis square; %xlim([-80 80]); ylim([-80 80]);
plot(thetaL(ap1),phiL(ap1),'go'); plot(thetaL(ap2),phiL(ap2),'go');
plot(thetaL(ap3),phiL(ap3),'go'); plot(thetaL(ap4),phiL(ap4),'go');
plot(thetaL(ap5),phiL(ap5),'go'); plot(thetaL(ap6),phiL(ap6),'go');
plot(thetaL(ap7),phiL(ap7),'go'); plot(thetaL(ap8),phiL(ap8),'go');
xlabel ('L theta'); ylabel('L phi');

