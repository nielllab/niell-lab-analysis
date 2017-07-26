%close all
% clear all

 %captureBatch

for f = 1:length(files)
    lighting(f) = files(f).lighting;
    contrast(f) = files(f).contrast;
    group(f)=files(f).group
    sex(f)=files(f).sex
    
    [ location(:,:,f) targ{f} body{f} r{f} thetaR{f} thetaSac{f} vhead{f} vbody{f} abody{f} vtarg{f} atarg{f}] = analyzeCapture([pathname files(f).trackpts],files(f).fps,files(f).scale) ;
    latency(f) = length(r{f});
    dist(f) = mean(r{f});
    title(['animal ',num2str(files(f).subj)])
end

%time to capture

d = [mean(latency(lighting==1& group==1)) mean(latency(lighting==0&group==1)) mean(latency(lighting==1&group==2)) mean(latency(lighting==0&group==2)) mean(latency(group==7 & lighting==1)) mean(latency(group==7 & lighting==0)) mean(latency(group==4 & lighting==1)) mean(latency(group==4 & lighting==0))]/60;
err = [std(latency(lighting==1& group==1))/sqrt(sum(lighting==1& group==1)) std(latency(lighting==0& group==1))/sqrt(sum(lighting==0& group==1)) std(latency(lighting==1& group==2))/sqrt(sum(lighting==1& group==2)) std(latency(lighting==0& group==2))/sqrt(sum(lighting==0& group==2)) std(latency(lighting==1& group==7))/sqrt(sum(lighting==1& group==7)) std(latency(lighting==0& group==7))/sqrt(sum(lighting==0& group==7)) std(latency(group==4& lighting==1))/sqrt(sum(group==4&lighting==1)) std(latency(group==4 & lighting==0))/sqrt(sum(group==4 & lighting==0)) ]/60;
figure
barweb(d,err)
ylabel('time to caputure (s)')
legend('light','dark', 'ES light','ES dark')


%%extract measures as a function of test condition, light vs dark etc.

range_L = vertcat(r{lighting==1& group==1});
range_D = vertcat(r{lighting==0 & group==2});
distTarg_L=abs(vertcat(vtarg{lighting==1& group==1}));
distTarg_D=abs(vertcat(vtarg{lighting==0& group==2}));

% distTarg_L=distTarg_L./35; % put in terms of centimeters 35 pixels per centimeter and data come out in terms of pixels and speed
% distTarg_D=distTarg_D./35;
% 
% speedT_L=distTarg_L.*60; %cm/sec
% speedT_D=distTarg_D.*60;
% 
% speedT_L=speedT_L*10;%mm/sec more standard reporting insect locomotion
% speedT_D=speedT_D*10;

% figure
% nhist(speedT_D,'maxx',15,'smooth','pdf');hold on
% nhist(speedT_L,'maxx',15,'smooth','color','qualitative','pdf')


%range_ES=vertcat(r{group==3});
range_ES_L=vertcat(r{group==3 & lighting==1});
%range_ES_D=vertcat(r{group==2 & lighting==0});

% range_ESR_L=vertcat(r{group==7 & lighting==1});
% range_ESR_D=vertcat(r{group==7 & lighting==0});
% 
% range_bgd=vertcat(r{group==3 & lighting==1});
% 
% range_EP_L=vertcat(r{group==4 & lighting==1});
% range_EP_D=vertcat(r{group==4 & lighting==0});


%index those trial that are in the light versus those in the dark
idx_L=find(lighting==1 & group==1);
idx_D=find(lighting==0 & group==2);

idx_ESL=find(group==3 & lighting==1);


figure
nhist(range_L(range_L<35),'p'); ylim([0 0.40])
title 'range light'
xlabel ('cm')

figure
nhist(range_D(range_D<35),'p'); ylim([0 0.40])
title 'range dark'
xlabel ('cm')

%eye suture experiments
figure
nhist(range_ES_L(range_ES_L<35),'p'); ylim([0 0.30])
title 'range ES light'
xlabel ('cm')

% figure
% nhist(range_ES_D(range_ES_D<35),'p'); ylim([0 0.30])
% title 'range ES dark'
% xlabel ('cm')
% 
% figure
% nhist(range_ESR_L(range_ESR_L<35),'p'); ylim([0 0.30])
% title 'range ESR light'
% xlabel ('cm')

% figure
% nhist(range_ESR_D(range_ESR_D<35),'p'); ylim([0 0.30])
% title 'range ESR dark'
% xlabel ('cm')

%% Ear plug experiments
% figure
% nhist(range_EP_L(range_EP_L<35),'p'); ylim([0 0.30])
% title 'range EP light'
% xlabel ('cm')
% 
% figure
% nhist(range_EP_D(range_EP_D<35),'p'); ylim([0 0.30])
% title 'range EP dark'
% xlabel ('cm')

%%extract distribution of theta's as a function of test condition, light vs dark etc.

thetaR_L = vertcat(thetaR{lighting==1 & group==1});
thetaR_D = vertcat(thetaR{lighting==0 & group==2});
% thetaSac_L=vertcat(thetaSac{lighting==1 & group==1});
% thetaSac_D=vertcat(thetaSac{lighting==0 & group==1});


%thetaR_ES = vertcat(thetaR{group==2});
thetaR_ESL = vertcat(thetaR{group==3 & lighting==1});
%thetaR_ESD = vertcat(thetaR{group==2 & lighting==0});

% thetaR_EPL = vertcat(thetaR{group==4 & lighting==1});
% thetaR_EPD=vertcat(thetaR{group==4 & lighting==0});

%%put thetas in symetrical space and works out range differences between
%%light and dark
thetaR_L=mod(thetaR_L-180,360);
thetaR_D=mod(thetaR_D-180,360);
% thetaSac_L=mod(thetaSac_L-180,360);
% thetaSac_D=mod(thetaSac_D-180,360);
% 
% thetaR_ES=mod(thetaR_ES-180,360);
thetaR_ESL=mod(thetaR_ESL-180,360);
% thetaR_ESD=mod(thetaR_ESD-180,360);
% 
% thetaR_EPL=mod(thetaR_EPL-180,360);
% thetaR_EPD=mod(thetaR_EPD-180,360);

figure
nhist(thetaR_L(range_L>5 & range_L<25),'p');ylim([0 0.3]);
figure
nhist(thetaR_D(range_D>5 & range_D<25),'p');ylim([0 0.3]);


%all eye sutured animals
% figure
% nhist(thetaR_ES(range_ES>5 & range_ES<15),'p');ylim([0 0.2]);
% 
figure
nhist(thetaR_ESL(range_ES_L>5 & range_ES_L<25),'p');ylim([0 0.3])
% % figure
% % nhist(thetaR_ESD(range_ES_D>5 & range_ES_D<15),'smooth','pdf');
% 
% %ear plugged animals
% figure
% nhist(thetaR_EPL(range_EP_L>5 & range_EP_L<15),'p');ylim([0 0.2]);
% 
% figure
% nhist(thetaR_EPD(range_EP_D>5 & range_EP_D<15),'p');ylim([0 0.2]);

figure
nhist(thetaR_D(range_D<5 ),'p')
title 'thetaR dark <10 cm'

figure
nhist(thetaR_L(range_L<5),'p')
title 'thetaR light <10cm'

% figure
% nhist(thetaR_ESD(range_ES_D<5 ),'p');ylim([0 0.10])
% hold on
% nhist(thetaR_ESL(range_ES_L<5),'p');ylim([0 0.10] )
% title 'thetaR light <10cm'
% 
% figure
% nhist(thetaR_EPD(range_EP_D<5 ),'p');ylim([0 0.10])
% hold on
% nhist(thetaR_EPL(range_EP_L<5),'p');ylim([0 0.10] )
% title 'thetaR light <10cm'

%% get rid of points where it crosses 360deg boundary (for plotting tracks)

% for blk6 mice in the dark
d=diff(thetaR_D); dd= diff(range_D);
thetaR_D_clean = thetaR_D;
thetaR_D_clean(abs(d)>180 | abs(dd)>200)=NaN;

figure
plot(range_D,thetaR_D_clean) ; ylim([-10 370])

figure
plot(range_D,thetaR_D_clean) ; ylim([-10 370])


%for blk6 mice in the light
d=diff(thetaR_L);  dd= diff(range_L);
thetaR_L_clean = thetaR_L;
thetaR_L_clean(abs(d)>180 | abs(dd)>200)=NaN;

figure
plot(range_L,thetaR_L_clean) ; ylim([-10 370])

%for all animals eye sutured 
% d=diff(thetaR_ES); dd= diff(range_ES);
% thetaR_ES_clean = thetaR_ES;
% thetaR_ES_clean(abs(d)>180 | abs(dd)>200)=NaN;
% 
% figure
% plot(range_ES,thetaR_ES_clean) ; ylim([-10 370]);
% 
% %for eye suture in light
d=diff(thetaR_ESL); dd= diff(range_ES_L);
thetaR_ESL_clean = thetaR_ESL;
thetaR_ESL_clean(abs(d)>180 | abs(dd)>200)=NaN;
% 
% figure
% plot(range_ES_L,thetaR_ESL_clean) ; ylim([-10 370]);
% 
% %for eye sutured mice in the dark
% d=diff(thetaR_ESD); dd= diff(range_ES_D);
% thetaR_ESD_clean = thetaR_ESD;
% thetaR_ESD_clean(abs(d)>180 | abs(dd)>200)=NaN;
% 
% figure
% plot(range_ES_D,thetaR_ESD_clean) ; ylim([-10 370])
% 
% %for Ear plugged in light
% d=diff(thetaR_EPL); dd= diff(range_EP_L);
% thetaR_EPL_clean = thetaR_EPL;
% thetaR_EPL_clean(abs(d)>180 | abs(dd)>200)=NaN;
% 
% figure
% plot(range_EP_L,thetaR_EPL_clean) ; ylim([-10 370]);
% 
% %for Ear plugged mice in the dark
% d=diff(thetaR_EPD); dd= diff(range_EP_D);
% thetaR_EPD_clean = thetaR_EPD;
% thetaR_EPD_clean(abs(d)>180 | abs(dd)>200)=NaN;
% 
% figure
% plot(range_EP_D,thetaR_EPD_clean) ; ylim([-10 370])


%% index: approaches, departures, trial edges, ranges around contact events (contact triggered average)
 
range_end=[];
n=0;

  for i=1:length(idx_L);
      clear range
      range=size(r{idx_L(i)}(:,1));
      range_end(i,:)=n+range(:,1);
      n=range_end(i);
  end

 range=range_L;
 TH=2.6;
 fps=60;
 
[range_smooth, contact, move_away, hits, approachStart peri_approach_start, escape, time_recovery, contactDwell_time, numEscape_trial]=rangeAnalysis(range,TH,range_end,fps);

contact_L=contact;
move_away_L=move_away;
range_smooth_L=range_smooth;
hits_L=hits;
peri_approach_start_L=peri_approach_start;
escape_L=escape;
time_recovery_L=time_recovery;
contactDwell_time_L=contactDwell_time;
numEscape_trial_L=numEscape_trial;

%% Plot approach triggered averages (ATA)
  ATA_range_L={}; ATA_theta_L={};pre_hit=peri_approach_start(~isnan(peri_approach_start));
  thetaR_L_clean_smooth=conv(thetaR_L_clean,ones(1,5),'same')/5;
 
figure
for i=1:length(pre_hit);
    
      ATA_range_L{i}=range_smooth(pre_hit(i):hits(i));
      ATA_theta_L{i}= thetaR_L_clean_smooth(pre_hit(i):hits(i));
      
      plot(ATA_range_L{i}(:,:),ATA_theta_L{i}(:,:)) ; ylim([-10 370]); hold on
end

figure
for i=1:length(ATA_range_L);
     
     x = -ATA_range_L{i}(:,:).*cos(ATA_theta_L{i}(:,:)*pi/180); 
     y = ATA_range_L{i}(:,:).*sin(ATA_theta_L{i}(:,:)*pi/180);
     plot(y,x); hold on; plot(0,0,'r*');
     axis ([-40 40 -40 40])
end


%% Dark data

range_end=[];
n=0;

  for i=1:length(idx_D);
      clear range
      range=size(r{idx_D(i)}(:,1));
      range_end(i,:)=n+range(:,1);
      n=range_end(i);
  end

 range=range_D;
 TH=2.6;
 fps=60;

[range_smooth, contact, move_away, hits, peri_approach_start, escape, time_recovery, contactDwell_time, numEscape_trial]=rangeAnalysis(range,TH,range_end,fps);

 
contact_D=contact;
move_away_D=move_away;
range_smooth_D=range_smooth;
hits_D=hits;
peri_approach_start_D=peri_approach_start;
escape_D=escape;
time_recovery_D=time_recovery;
contactDwell_time_D=contactDwell_time;
numEscape_trial_D=numEscape_trial;

clear app pre_hit 
  ATA_range_D={}; ATA_theta_D={};pre_hit=peri_approach_start(~isnan(peri_approach_start));
  thetaR_D_clean_smooth=conv(thetaR_D_clean,ones(1,5),'same')/5;

figure
for i=1:length(pre_hit);
      ATA_range_D{i}=range_smooth_D(pre_hit(i):hits(i));
      ATA_theta_D{i}= thetaR_D_clean(pre_hit(i):hits(i));
      
      plot(ATA_range_D{i}(:,:),ATA_theta_D{i}(:,:)) ; ylim([-10 370]); hold on
end

figure
for i=1:length(ATA_range_D);
  
     x = -ATA_range_D{i}(:,:).*cos(ATA_theta_D{i}(:,:)*pi/180); 
     y = ATA_range_D{i}(:,:).*sin(ATA_theta_D{i}(:,:)*pi/180);
     plot(y,x); hold on; plot(0,0,'r*');
     axis ([-40 40 -40 40])
      
end


%%
figure
nhist(time_recovery_D,'p',1:30);ylim([0 1]);xlim([-5 45]);
hold on
figure
nhist(time_recovery_L,'p',1:30);ylim([0 1]);xlim([-5 45]);
title 'distribution of recovery times after an escape'


%% Eye suture light
range_end=[];
n=0;

  for i=1:length(idx_ESL);
      clear range
      range=size(r{idx_ESL(i)}(:,1));
      range_end(i,:)=n+range(:,1);
      n=range_end(i);
  end

 range=range_ES_L;
 TH=2.6;
 fps=60;

[range_smooth, contact, move_away, hits, peri_approach_start, escape, time_recovery, contactDwell_time, numEscape_trial]=rangeAnalysis(range,TH,range_end,fps);

 
contact_ESL=contact;
move_away_ESL=move_away;
range_smooth_ESL=range_smooth;
hits_ESL=hits;
peri_approach_start_ESL=peri_approach_start;
escape_ESL=escape;
time_recovery_ESL=time_recovery;
contactDwell_time_ESL=contactDwell_time;
numEscape_trial_ESL=numEscape_trial;

clear app pre_hit 
  ATA_range_ESL={}; ATA_theta_ESL={};pre_hit=peri_approach_start_ESL(~isnan(peri_approach_start_ESL));
  thetaR_ESL_clean_smooth=conv(thetaR_ESL_clean,ones(1,5),'same')/5;

figure
for i=1:length(pre_hit);
      ATA_range_ESL{i}=range_smooth_ESL(pre_hit(i):hits(i));
      ATA_theta_ESL{i}= thetaR_ESL_clean(pre_hit(i):hits(i));
      
      plot(ATA_range_ESL{i}(:,:),ATA_theta_ESL{i}(:,:)) ; ylim([-10 370]); hold on
end

figure
for i=1:length(ATA_range_ESL);
  
     x = -ATA_range_ESL{i}(:,:).*cos(ATA_theta_ESL{i}(:,:)*pi/180); 
     y = ATA_range_ESL{i}(:,:).*sin(ATA_theta_ESL{i}(:,:)*pi/180);
     plot(y,x); hold on; plot(0,0,'r*');
     axis ([-40 40 -40 40])
end


%% Eye suture dark data
range_end=[];
n=0;

  for i=1:length(idx_ESD);
      clear range
      range=size(r{idx_ESD(i)}(:,1));
      range_end(i,:)=n+range(:,1);
      n=range_end(i);
  end

 range=range_ES_D;
 TH=2.6;
 fps=60;
 
[range_smooth, contact, move_away, hits, peri_approach_start, escape, time_recovery, contactDwell_time, numEscape_trial]=rangeAnalysis(range,TH,range_end,fps);

contact_ESD=contact;
move_away_ESD=move_away;
range_smooth_ESD=range_smooth;
hits_ESD=hits;
peri_approach_start_ESD=peri_approach_start;
escape_ESD=escape;
time_recovery_ESD=time_recovery;
contactDwell_time_ESD=contactDwell_time;
numEscape_trial_ESD=numEscape_trial;
%% Plot approach triggered averages (ATA)
  ATA_range_ESD={}; ATA_theta_ESD={};pre_hit=peri_approach_start(~isnan(peri_approach_start));
  thetaR_ESD_clean_smooth=conv(thetaR_ESD_clean,ones(1,5),'same')/5;
 
figure
for i=1:length(pre_hit);
    
      ATA_range_ESD{i}=range_smooth_ESD(pre_hit(i):hits(i));
      ATA_theta_ESD{i}= thetaR_ESD_clean_smooth(pre_hit(i):hits(i));
      
      plot(ATA_range_ESD{i}(:,:),ATA_theta_ESD{i}(:,:)) ; ylim([-10 370]); hold on
end

figure
for i=1:length(ATA_range_ESD);
     
     x = -ATA_range_ESD{i}(:,:).*cos(ATA_theta_ESD{i}(:,:)*pi/180); 
     y = ATA_range_ESD{i}(:,:).*sin(ATA_theta_ESD{i}(:,:)*pi/180);
     plot(y,x); hold on; plot(0,0,'r*');
     axis ([-40 40 -40 40])
end


  %% eye suture all animals
  range_end=[];
n=0;

  for i=1:length(idx_ES_all);
      clear range
      range=size(r{idx_ES_all(i)}(:,1));
      range_end(i,:)=n+range(:,1);
      n=range_end(i);
  end
  
 range=range_ES;
 TH=2.6;
 fps=60;
 
[range_smooth, contact, move_away, hits, peri_approach_start, escape, time_recovery, contactDwell_time, numEscape_trial]=rangeAnalysis(range,TH,range_end,fps);

contact_ES=contact;
move_away_ES=move_away;
range_smooth_ES=range_smooth;
hits_ES=hits;
peri_approach_start_ES=peri_approach_start;
escape_ES=escape;
time_recovery_ES=time_recovery;
contactDwell_time_ES=contactDwell_time;
numEscape_trial_ES=numEscape_trial;
%% Plot approach triggered averages (ATA)
  ATA_range_ES={}; ATA_theta_ES={};pre_hit=peri_approach_start(~isnan(peri_approach_start));
  thetaR_ES_clean_smooth=conv(thetaR_ES_clean,ones(1,5),'same')/5;
 
figure
for i=1:length(pre_hit);
    
      ATA_range_ES{i}=range_smooth_ES(pre_hit(i):hits(i));
      ATA_theta_ES{i}= thetaR_ES_clean_smooth(pre_hit(i):hits(i));
      
      plot(ATA_range_ES{i}(:,:),ATA_theta_ES{i}(:,:)) ; ylim([-10 370]); hold on
end

figure
for i=1:length(ATA_range_ES);
     
     x = -ATA_range_ES{i}(:,:).*cos(ATA_theta_ES{i}(:,:)*pi/180); 
     y = ATA_range_ES{i}(:,:).*sin(ATA_theta_ES{i}(:,:)*pi/180);
     plot(y,x); hold on; plot(0,0,'r*');
     axis ([-40 40 -40 40])
end


%% Ear plugged data light

range_end=[];
n=0;

  for i=1:length(idx_EPL);
      clear range
      range=size(r{idx_EPL(i)}(:,1));
      range_end(i,:)=n+range(:,1);
      n=range_end(i);
  end

 range=range_EP_L;
 TH=2.6;
 fps=60;
 
[range_smooth, contact, move_away, hits, peri_approach_start, escape, time_recovery, contactDwell_time, numEscape_trial]=rangeAnalysis(range,TH,range_end,fps);

contact_EPL=contact;
move_away_EPL=move_away;
range_smooth_EPL=range_smooth;
hits_EPL=hits;
peri_approach_start_EPL=peri_approach_start;
escape_EPL=escape;
time_recovery_EPL=time_recovery;
contactDwell_time_EPL=contactDwell_time;
numEscape_trial_EPL=numEscape_trial;
%% Plot approach triggered averages (ATA)
  ATA_range_EPL={}; ATA_theta_EPL={};pre_hit=peri_approach_start(~isnan(peri_approach_start));
  thetaR_EPL_clean_smooth=conv(thetaR_EPL_clean,ones(1,5),'same')/5;
 
figure
for i=1:length(pre_hit);
    
      ATA_range_EPL{i}=range_smooth_EPL(pre_hit(i):hits(i));
      ATA_theta_EPL{i}= thetaR_EPL_clean_smooth(pre_hit(i):hits(i));
      
      plot(ATA_range_EPL{i}(:,:),ATA_theta_EPL{i}(:,:)) ; ylim([-10 370]); hold on
end

figure
for i=1:length(ATA_range_EPL);
     
     x = -ATA_range_EPL{i}(:,:).*cos(ATA_theta_EPL{i}(:,:)*pi/180); 
     y = ATA_range_EPL{i}(:,:).*sin(ATA_theta_EPL{i}(:,:)*pi/180);
     plot(y,x); hold on; plot(0,0,'r*');
     axis ([-40 40 -40 40])
end

%% Ear plugged dark 

range_end=[];
n=0;

  for i=1:length(idx_EPD);
      clear range
      range=size(r{idx_EPD(i)}(:,1));
      range_end(i,:)=n+range(:,1);
      n=range_end(i);
  end

 range=range_EP_D;
 TH=2.6;
 fps=60;
 
[range_smooth, contact, move_away, hits, peri_approach_start, escape, time_recovery, contactDwell_time, numEscape_trial]=rangeAnalysis(range,TH,range_end,fps);

contact_EPD=contact;
move_away_EPD=move_away;
range_smooth_EPD=range_smooth;
hits_EPD=hits;
peri_approach_start_EPD=peri_approach_start;
escape_EPD=escape;
time_recovery_EPD=time_recovery;
contactDwell_time_EPD=contactDwell_time;
numEscape_trial_EPD=numEscape_trial;
%% Plot approach triggered averages (ATA)
  ATA_range_EPD={}; ATA_theta_EPD={};pre_hit=peri_approach_start(~isnan(peri_approach_start));
  thetaR_EPD_clean_smooth=conv(thetaR_EPD_clean,ones(1,5),'same')/5;
 
figure
for i=1:length(pre_hit);
    
      ATA_range_EPD{i}=range_smooth_EPD(pre_hit(i):hits(i));
      ATA_theta_EPD{i}= thetaR_EPD_clean_smooth(pre_hit(i):hits(i));
      
      plot(ATA_range_EPD{i}(:,:),ATA_theta_EPD{i}(:,:)) ; ylim([-10 370]); hold on
end

figure
for i=1:length(ATA_range_EPD);
     
     x = -ATA_range_EPD{i}(:,:).*cos(ATA_theta_EPD{i}(:,:)*pi/180); 
     y = ATA_range_EPD{i}(:,:).*sin(ATA_theta_EPD{i}(:,:)*pi/180);
     plot(y,x); hold on; plot(0,0,'r*');
     axis ([-40 40 -40 40])
end

%% angle of neck (head relative to body)
neck_angle_all = vertcat(head_saccade{:})*180/pi;
neck_angle_all = mod(neck_angle_all+180,360)-180;
figure
nhist(neck_angle_all)

neck_L = vertcat(head_saccade{lighting==1})*180/pi; neck_L = mod(neck_L+180,360)-180;
neck_D = vertcat(head_saccade{lighting==0})*180/pi; neck_D = mod(neck_D+180,360)-180;




%%group data for head to target angle relative to distance from prey
figure
binr=4:4:32
bint=0:10:360
hrt=hist2(range_smooth_L(range_smooth_L>4),thetaR_L_clean_smooth(range_smooth_L>4),binr,bint);
hrt= hrt/sum(hrt(:));
imagesc(flipud(hrt));
xlabel 'distance'
ylabel 'thetaT'
title 'light'

figure
binr=4:4:32
bint=0:10:360
%hrt=hist2(r_samp, theta_samp,binr,bint)
hrt=hist2(range_smooth_D(range_smooth_D>4),thetaR_D_clean_smooth(range_smooth_D>4),binr,bint);
hrt= hrt/sum(hrt(:));
imagesc(flipud(hrt));
xlabel 'distance'
ylabel 'thetaT'
title 'dark'

figure
binr=4:4:32
bint=0:10:360
%hrt=hist2(r_samp, theta_samp,binr,bint)
hrt=hist2(range_smooth_ESL(range_smooth_ESL>4),thetaR_ESL_clean_smooth(range_smooth_ESL>4),binr,bint);
hrt= hrt/sum(hrt(:));
imagesc(flipud(hrt));
xlabel 'distance'
ylabel 'thetaT'
title 'ES light'

figure
binr=4:4:32
bint=0:10:360
%hrt=hist2(r_samp, theta_samp,binr,bint)
hrt=hist2(range_smooth_ESD(range_smooth_ESD>4),thetaR_ESD_clean_smooth(range_smooth_ESD>4),binr,bint);
hrt= hrt/sum(hrt(:));
imagesc(flipud(hrt));
xlabel 'distance'
ylabel 'thetaT'
title 'ES Dark'


figure
binr=4:4:32
bint=0:10:360
%hrt=hist2(r_samp, theta_samp,binr,bint)
hrt=hist2(range_smooth_EPL(range_smooth_EPL>4),thetaR_EPL_clean_smooth(range_smooth_EPL>4),binr,bint);
hrt= hrt/sum(hrt(:));
imagesc(flipud(hrt));
xlabel 'distance'
ylabel 'thetaT'
title 'EP light'

figure
binr=4:4:32
bint=0:10:360
%hrt=hist2(r_samp, theta_samp,binr,bint)
hrt=hist2(range_smooth_EPD(range_smooth_EPD>4),thetaR_EPD_clean_smooth(range_smooth_EPD>4),binr,bint);
hrt= hrt/sum(hrt(:));
imagesc(flipud(hrt));
xlabel 'distance'
ylabel 'thetaT'
title 'EP dark'
%%%target angles relative to head centered space as a function of light

figure
subplot(1,3,1);
location(isnan(location)) = 0 ;
%%needs to index those files with only zeros and subtract out figure out
%%nanmean for this plot
imagesc(flipud(mean(location(:,:,lighting==0 ),3)),[0 0.04]); hold on; plot(17,17,'ro','Markersize',4)
title('dark')
axis square
subplot(1,3,2);
imagesc(flipud(mean(location(:,:,lighting==1),3)),[0 0.04]);hold on; plot(17,17,'ro','Markersize',4)
title('light')
axis square
subplot(1,3,3);
imagesc(flipud(mean(location(:,:,group==2),3)),[0 0.04]);hold on; plot(17,17,'ro','Markersize',4)
title('blind')
axis square

%this cancels out the NANs that weight the relative intensity, there is a
%better function that imagesc that ignores NaNs, check Martin's code for
%this function

k=max(location);
idx_dist=find(max(k>0));
lighting_dist=lighting(idx_dist);
loc_dist=location(:,:,idx_dist);

figure
subplot(1,2,1);
%%needs to index those files with only zeros and subtract out
imagesc(flipud(mean(loc_dist(:,:,lighting_dist==0 ),3)),[0 0.05]);
title('dark')
axis equal

subplot(1,2,2);
imagesc(flipud(mean(loc_dist(:,:,lighting_dist==1),3)),[0 0.05]);
title('light')
axis equal


