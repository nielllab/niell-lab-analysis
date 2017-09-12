%% non shaded version of range vs theta plots overlaid for the four main conditions
figure; hold on
plot(dbin,nanmean(avg_theta_L,3),'b','LineWidth',2);hold on
plot(dbin,(nanmean(avg_theta_L,3)-stdE_L),'b','LineWidth',1);
plot(dbin,(nanmean(avg_theta_L,3)+ stdE_L),'b','LineWidth',1);
set(gca,'xdir','reverse','XLim',[0 40],'YLim',[0 180])

hold on
plot(dbin,nanmean(avg_theta_D,3),'k','LineWidth',2);hold on
plot(dbin,(nanmean(avg_theta_D,3)-stdE_D),'k','LineWidth',1);
plot(dbin,(nanmean(avg_theta_D,3)+ stdE_D),'k','LineWidth',1);
set(gca,'xdir','reverse','XLim',[0 40],'YLim',[0 180])

plot(dbin,nanmean(avg_theta_EPL,3),'g','LineWidth',2);hold on
plot(dbin,(nanmean(avg_theta_EPL,3)-stdE_EPL),'g','LineWidth',1);
plot(dbin,(nanmean(avg_theta_EPL,3)+ stdE_EPL),'g','LineWidth',1);
set(gca,'xdir','reverse','XLim',[0 40],'YLim',[0 180])

plot(dbin,nanmean(avg_theta_EPD,3),'r','LineWidth',2);hold on
plot(dbin,(nanmean(avg_theta_EPD,3)-stdE_EPD),'r','LineWidth',1);
plot(dbin,(nanmean(avg_theta_EPD,3)+ stdE_EPD),'r','LineWidth',1);
set(gca,'xdir','reverse','XLim',[0 40],'YLim',[0 180])

%% shaded error bar of range vs theta
figure; hold on
shadedErrorBar(dbin,nanmean(avg_theta_app,3),stdE_L,'b');hold on
set(gca,'xdir','reverse','XLim',[0 40],'YLim',[0 180])

hold on
shadedErrorBar(dbin,nanmean(avg_theta_app,3),stdE_L,'k');hold on
shadedErrorBar(dbin,nanmean(avg_theta_app,3),stdE_L,'g');hold on
shadedErrorBar(dbin,nanmean(avg_theta_app,3),stdE_L,'r');hold on


%% theta distributions over all session overalid for eye suture controls conditions
figure
plot(thetabins,nanmean(PreyE_mid_L,2),'b','LineWidth',3); hold on%eye suture light
plot(thetabins,nanmean(PreyE_mid_D,2),'k','LineWidth',3);%eye suture light
plot(thetabins,nanmean(PreyE_mid_ESL,2),'c','LineWidth',3); %eye suture light
plot(thetabins,nanmean(PreyE_mid_ESD,2),'m','LineWidth',3); %eye suture 


