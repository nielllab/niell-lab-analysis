%function analyze_excel_rfdata
clear all
fname = 'allunits_new0108';

t_width = xlsread(fname,1,'F2:F242');
t2peak = xlsread(fname,1,'G2:G242');
peaktroughamp = xlsread(fname,1,'H2:H242');
snip_endheight = xlsread(fname,1,'I2:I242');
ISIpeak = xlsread(fname,1,'J2:J242');

barOSI =  xlsread(fname,1,'R2:R242');
bartheta = xlsread(fname,1,'S2:S242');
barpeak = xlsread(fname,1,'V2:V242');
barspont = xlsread(fname,1,'W2:W242');
barDSI = xlsread(fname,1,'U2:U242');
bartheta_width = xlsread(fname,1,'T2:T242');
bartheta_width=abs(bartheta_width);
bar_rfwidth = xlsread(fname,1,'X2:X242');
barpeak(barpeak<0)=0;

driftOSI = xlsread(fname,1,'Y2:Y242');
drifttheta = xlsread(fname,1,'Z2:Z242');
driftpeak = xlsread(fname,1,'AF2:AF242');
driftspont = xlsread(fname,1,'AG2:AG242');
driftSF = xlsread(fname,1,'AC2:AC242');
driftF1F0 = xlsread(fname,1,'AE2:AE242')/2;
driftDSI = xlsread(fname,1,'AB2:AB242');
drifttheta_width =xlsread(fname,1,'AA2:AA242');
drifttheta_width=abs(drifttheta_width);

driftpeak(driftpeak<0)=0;

driftSF_width = xlsread(fname,1,'AD2:AD242');

flash_on = xlsread(fname,1,'BE2:BE242');

cp_tf_tuning  = xlsread(fname,1,'AZ2:BC242');
sb_spont = xlsread(fname,1,'AY2:AY242');
wx = xlsread(fname,1,'AV2:AV242');
wy = xlsread(fname,1,'AW2:AW242');
sb_x = xlsread(fname,1,'AT2:AT242');
sb_y = xlsread(fname,1,'AU2:AU242');

contrast_mod = xlsread(fname,1,'AL2:AL242');
contrast_phase = xlsread(fname,1,'AM2:AM242');
movieN = xlsread(fname,1,'AJ2:AJ242');
sta_freq = xlsread(fname,1,'AI2:AI242');
sta_found = xlsread(fname,1,'AR2:AR242');
sta_theta = xlsread(fname,1,'AH2:AH242');
sta_halfcontrast = xlsread(fname,1,'AO2:AO242');
sta_duration = xlsread(fname,1,'AN2:AN242');
sta_slope = xlsread(fname,1,'AP2:AP242');
sta_adaptation = xlsread(fname,1,'AQ2:AQ242');
sta_subunits =  xlsread(fname,1,'AS2:AS242');

meanspont = 0.333*(driftspont+barspont+sb_spont);
monitor_offset = xlsread(fname,1,'BT2:BT242');
monitor_offset(isnan(monitor_offset))=0;

simple_est =  xlsread(fname,1,'O2:O242');
ds_est = xlsread(fname,1,'N2:N242');
orient_est =xlsread(fname,1,'M2:M242');
responsive_est = xlsread(fname,1,'Q2:Q242');
inh_est = xlsread(fname,1,'P2:P242');
exc = (inh_est==0);
inh = (inh_est==1);

bar_osicirc = xlsread(fname,1,'BP2:BP242');
drift_osicirc = xlsread(fname,1,'BR2:BR242');
bar_dsicirc =xlsread(fname,1,'BQ2:BQ242');
drift_dsicirc = xlsread(fname,1,'BS2:BS242');

layer = xlsread(fname,1,'B2:B242');

layer(layer==2.5 ) =2;
layer(layer==3.5)=3;
layer(layer==4.5)=5;
layer(layer==5.5) = 5;
layer(layer>6)=NaN;

layertwo = layer<3;
layerthree= layer==3;
layerfour = layer==4;
layerfive = layer==5;
layersix=layer==6;
layertwothree = layertwo | layerthree;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% reading in data is finished here
figure
plot3(peaktroughamp(exc),snip_endheight(exc),t2peak(exc),'go');
hold on
plot3(peaktroughamp(inh),snip_endheight(inh),t2peak(inh),'ro');

figure
plot(peaktroughamp(exc),snip_endheight(exc),'go','MarkerSize',8,'LineWidth',1);
hold on
plot(peaktroughamp(inh),snip_endheight(inh),'bo','MarkerSize',8,'LineWidth',1);
axis square
axis([0.5 3.5 0 1])

figure
plot(t_width(exc)/25,t2peak(exc)/25,'go');
hold on
plot(t_width(inh)/25,t2peak(inh)/25,'ro');

expect =sum(inh==1)/sum(inh==0 | inh==1)

for l = 2:6
    l
    sum(inh&layer==l)
    sum(inh&layer==l)/(sum(inh&layer==l) + sum(exc&layer==l))
    (sum(inh&layer==l) - expect*sum(layer==l))^2 / (expect*sum(layer==l)) + ...
        (sum(exc&layer==l) - (1-expect)*sum(layer==l))^2 / ((1-expect)*sum(layer==l))
end

clear data; 
data(:,1) = peaktroughamp;
data(:,1) = (data(:,1)-min(data(:,1)))/(max(data(:,1))-min(data(:,1)));
data(:,2) = snip_endheight;
data(:,2) = (data(:,2)-min(data(:,2)))/(max(data(:,2))-min(data(:,2)));
data(:,3) = t2peak;
data(:,3) = (data(:,3)-min(data(:,3)))/(max(data(:,3))-min(data(:,3)));
used = find(~isnan(sum(data,2)))

[cl  c]= kmeans(zscore(data(used,1:3)),2);
clear clust
clust(used)=cl;
clust = clust';
% data=data(used,:);
% clust(used) = clusterdata(data,4)
figure
plot3(peaktroughamp(clust==1&exc),snip_endheight(clust==1&exc),t2peak(clust==1&exc),'go');
hold on
plot3(peaktroughamp(clust==1&inh),snip_endheight(clust==1&inh),t2peak(clust==1&inh),'ro');
hold on
plot3(peaktroughamp(clust==2&inh),snip_endheight(clust==2&inh),t2peak(clust==2&inh),'r*');
hold on
plot3(peaktroughamp(clust==2&exc),snip_endheight(clust==2&exc),t2peak(clust==2&exc),'g*');

used =~isnan(sum(data,2))

wv = xlsread(fname,1,'BU2:CM241');
figure
hold on
for i=1:size(wv,1)
        if inh(i)==1;
        colorstr='r';
    elseif inh(i)==0;
        colorstr='g';
    else colorstr='.';
        end
    if ~used(i)
        colorstr='b';
    end
    [y j] = min(wv(i,:));
    j
    if j ==5
        wvshift(i,:) = wv(i,1:17);
        plot(wv(i,1:17),colorstr);
    elseif j==6
           wvshift(i,:) = wv(i,2:18);
           plot(wv(i,2:18),colorstr);
    elseif j==7
        wvshift(i,:) = wv(i,3:19);
        plot(wv(i,3:19),colorstr);
    else
        wvshift(i,:) = wv(i,1:17);
    end

end
figure
imagesc(wvshift)

[s coeff u] = fastica(wvshift(find(used),1:17)','numOfIC',2,'lastEig',2,'stabilization','on','g','tanh','approach','symm');
figure
plot(coeff);

figure
plot(s(1,find(inh(used)==1)),s(2,find(inh(used)==1)),'ro')
hold on
plot(s(1,find(inh(used)==0)),s(2,find(inh(used)==0)),'go')

figure
plot3(s(1,find(inh==1)),s(2,find(inh==1)),s(3,find(inh==1)),'ro')
hold on
plot3(s(1,find(inh==0)),s(2,find(inh==0)),s(3,find(inh==0)),'go')

newendheight = wvshift(:,17)./max(wvshift,[],2);
figure
plot(newendheight,snip_endheight,'o')
figure
plot(newendheight(inh==1),max(wvshift(inh==1),[],2),'ro')
hold on
plot(newendheight(inh==0),max(wvshift(inh==0),[],2),'go')

used=1;

used =~isnan(sum(data,2))
newratio = (wvshift(:,11)+wv(:,12));
newratio(newratio>3)=3;
figure
hist(newratio,[0:0.05:1])
figure
plot3(peaktroughamp(inh==1&used),newratio(inh==1&used),t2peak(inh==1&used),'ro');
hold on
plot3(peaktroughamp(inh==0&used),newratio(inh==0&used),t2peak(inh==0&used),'go');
hold on
plot3(peaktroughamp(~used),newratio(~used),t2peak(~used),'bo');

new_endheight = wvshift(:,17)-wvshift(:,15);
new_midpoint = wvshift(:,11)+wvshift(:,12);
new_max = max(wvshift,[],2);

figure
plot(new_endheight(inh==1),new_max(inh==1),'ro');
hold on
plot(new_endheight(inh==0&used),new_max(inh==0&used),'go');

used =1;
figure
plot3(new_endheight(inh==1&used),new_midpoint(inh==1&used),new_max(inh==1&used),'ro');
hold on
plot3(new_endheight(inh==0&used),new_midpoint(inh==0&used),new_max(inh==0&used),'go');

figure
plot(newratio(inh==1&used),snip_endheight(inh==1&used),'ro');
hold on
plot(newratio(inh==0&used),snip_endheight(inh==0&used),'go');


figure
plot(peaktroughamp(clust==1),snip_endheight(clust==1),'go','MarkerSize',8,'LineWidth',1);
hold on
plot(peaktroughamp(clust==2),snip_endheight(clust==2),'bo','MarkerSize',8,'LineWidth',1);
axis square
axis([0.5 3.5 0 1])
hold on
plot(c(1:2,1:2)','o')




%%%%%%%%%%%%%%%%%%%%%% orientation
figure
% barpeak(barpeak<0.5)=0.5;
% driftpeak(driftpeak<1)=0.5;
plot(log10(barpeak),log10(driftpeak),'o')
figure
plot(barpeak,driftpeak,'o')

drift_thresh = 2;
bar_thresh=4;
 dp = driftpeak;
 bp = barpeak;
% 
dp = dp-drift_thresh/2;
 bp = bp-bar_thresh/2;

dp(dp<0)=0;
bp(bp<0)=0;
meanOSI = (2*dp.*driftOSI + bp.*barOSI)./(2*dp+bp);
used = driftpeak>drift_thresh | barpeak>bar_thresh ;

%meanOSI = (2*meanOSI)./(1+meanOSI);

meanOSIcirc = (2*dp.*drift_osicirc + bp.*bar_osicirc)./(2*dp+bp);

figure
plot(meanOSI(used),meanOSIcirc(used),'o');


figure
plot(meanOSI(exc&used),driftF1F0(exc&used),'go');
hold on
plot(meanOSI(inh&used),driftF1F0(inh&used),'ro');
title('osi vs f1f0')

used = (driftpeak>drift_thresh | barpeak>bar_thresh) & ~isnan(meanOSI) &~isnan(meanOSIcirc) &responsive_est ==1;

layer_scatter(meanOSI,used,layer,exc,inh);
sum(used)
xlabel('mean osi','FontSize',16)

layer_scatter(meanOSIcirc,used,layer,exc,inh);
sum(used)
xlabel('mean osi (cv)','FontSize',16)

layer_bar_inh(meanOSI,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('mean OSI','FontSize',16);

[p h] = ranksum(meanOSI(used&layertwo&exc),meanOSI(used&layerthree&exc))

layer_bar(meanOSIcirc,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('mean OSI (cv)','FontSize',16);



meanOSI(meanOSI>1.05)=1.05;
meanOSI(meanOSI<0)=0;
figure
hist(meanOSI(used),[0.05:0.1:1.05]);

figure
hist(meanOSIcirc(used & (layertwothree | layerfour) & exc & responsive_est),[0.05:0.1:1.05]);


meanOSI(meanOSI>1.05)=1.05;
meanOSI(meanOSI<0)=0;
figure
hist(meanOSIcirc(used),[0.025:0.05:1.05]);

used = (driftpeak>drift_thresh | barpeak>bar_thresh) & ~isnan(meanOSI) &responsive_est ==1 & ~isnan(meanOSIcirc);

figure
plot(orient_est(used&exc==1),meanOSIcirc(used&exc==1),'o');
hold on
plot(orient_est(used&exc==0),meanOSIcirc(used&exc==0),'ro');


figure
plot(orient_est(used&exc==1),meanOSI(used&exc==1),'o');
hold on
plot(orient_est(used&exc==0),meanOSI(used&exc==0),'ro');

% for lay=2:5
%     out(lay) = mean(meanOSI(used&inh&(layer==lay)));
%     plot([out(lay) out(lay)],[lay-0.25 lay+0.25],'r','LineWidth',2)
% end
% %plot(out(2:5),2:5,'r','Linewidth',2);
% for lay=2:5
%     out(lay) = mean(meanOSI(used&exc&(layer==lay)));
%      plot([out(lay) out(lay)],[lay-0.25 lay+0.25],'g','LineWidth',2)
% end
% %plot(out(2:5),2:5,'g','Linewidth',2);


used = (driftpeak>drift_thresh*2 & barpeak>bar_thresh*2) & ~isnan(meanOSI) &responsive_est ==1;
figure
plot(driftOSI(used),barOSI(used),'o')
hold on
plot([0 1], [0 1])
axis equal

figure
plot(orient_est(used),meanOSI(used),'o');

bmock = rand(size(bartheta))*360;
dmock = rand(size(drifttheta))*360;

barthetamod = mod(bartheta,180);
driftthetamod = mod(drifttheta,180);
dt = barthetamod-driftthetamod;

barthetamod(dt>90) = barthetamod(dt>90)-180;
barthetamod(dt<-90) = barthetamod(dt<-90)+180;

used = driftpeak>2 & barpeak>4 & responsive_est==1 & orient_est==1 &driftOSI>0.4 & barOSI>0.4
figure
plot(driftthetamod(used),barthetamod(used),'*')
hold on
plot([0 180],[0 180])
axis equal
axis([-20 200 -20 200])
sum(used)
[r p] = corrcoef( driftthetamod(used), barthetamod(used))
r.^2

meantheta = (2*dp.*driftthetamod + bp.*barthetamod)./(2*dp+bp);

meantheta_width = (2*dp.*drifttheta_width + bp.*bartheta_width)./(2*dp+bp);
meantheta_width = 1.17*meantheta_width;
used = (meanOSI>0.4) &(driftpeak>3 | barpeak>6) &~isnan(meantheta_width);
layer_scatter(meantheta_width,used,layer,exc,inh);
sum(used)
xlabel('tuning width (deg)','FontSize',16)

%meantheta = (2*dp.*driftthetamod + bp.*barthetamod)./(2*dp+bp);



used = (meanOSI>0.4) &(driftpeak>2 | barpeak>4) &~isnan(meantheta_width);
layer_bar_inh(meantheta_width,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('mean OSI','FontSize',16);

layer_bar_median(meantheta_width,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('mean tuning width','FontSize',16);

figure
barweb([mean(meantheta_width(used&exc&layertwo)), mean(meantheta_width(used&exc&layerthree)), mean(meantheta_width(used&exc&layerfour)) ,mean(meantheta_width(used&exc&layerfive)),mean(meantheta_width(used&exc&layersix)); 0 0 0 0 0]', ...
    [std(meantheta_width(used&exc&layertwo))/sqrt(sum(used&exc&layertwo )), std(meantheta_width(used&exc&layerthree))/sqrt(sum(used&exc&layerthree)) ,  ...
    std(meantheta_width(used&exc&layerfour))/sqrt(sum(used&exc&layerfour )), std(meantheta_width(used&exc&layerfive))/sqrt(sum(used&exc&layerfive )) , std(meantheta_width(used&exc&layersix))/sqrt(sum(used&exc&layersix)); 0 0 0 0 0]', ...
    2,{'2' '3' '4' '5' '6'},[],[],[],[0 0 .75]);
set(gca,'xticklabel', {'2', '3', '4', '5' '6'})
ylabel('mean tuning width');
[h p] = ttest2(meantheta_width(used&exc&layertwo),meantheta_width(used&exc&(layerthree)))


%meantheta_width(meantheta_width>65)=65;
used = (meanOSI>0.5) &(driftpeak>3 | barpeak>6) &~isnan(meantheta_width);
figure
hist(meantheta_width(used),2.5:5:70)
median(meantheta_width(used))
upper_quartile(meantheta_width(used))
lower_quartile(meantheta_width(used))

used = (meanOSI>0.5) &(driftpeak>3 | barpeak>6) &~isnan(meantheta_width) ;
figure
plot(meanOSIcirc(used),meantheta_width(used),'o')


figure
used = meantheta_width<20  &(driftpeak>3 | barpeak>6) &~isnan(meantheta_width) ;
plot(meanOSIcirc(used),meanOSI(used),'go');
hold on
used = meantheta_width>20 & meantheta_width<40 &(driftpeak>3 | barpeak>6) &~isnan(meantheta_width) ;
plot(meanOSIcirc(used),meanOSI(used),'bo');
used = meantheta_width>40 &(driftpeak>3 | barpeak>6) &~isnan(meantheta_width) ;
plot(meanOSIcirc(used),meanOSI(used),'ro');

figure
plot(meantheta_width(used),meanOSI(used),'o')
%%%%%%%%%%%%%%%%%%direction selectivity


drift_thresh = 3;
bar_thresh=6;
 dp = driftpeak;
 bp = barpeak;
% 
dp = dp-drift_thresh/2;
 bp = bp-bar_thresh/2;

dp(dp<0)=0;
bp(bp<0)=0;
meanDSI = (2*dp.*driftDSI.*driftOSI + bp.*barDSI)./(2*dp+bp);   %%% multiply by drift OSI fixes error in compile analysis
meanDSIcirc = (2*dp.*drift_osicirc + bp.*bar_osicirc)./(2*dp+bp);

used = driftpeak>drift_thresh | barpeak>bar_thresh & orient_est==1;

figure
plot(meanDSI(used),meanDSIcirc(used),'o')

figure
plot(ds_est(used),meanDSIcirc(used),'o')

figure
plot(meanDSI(exc&used),driftF1F0(exc&used),'go');
hold on
plot(meanDSI(inh&used),driftF1F0(inh&used),'ro');
title('dsi vs f1f0')

used = (driftpeak>drift_thresh | barpeak>bar_thresh) & ~isnan(meanDSI) &responsive_est ==1;

layer_scatter(meanDSI,used,layer,exc,inh);
sum(used)
xlabel('mean dsi','FontSize',16)

layer_bar(meanDSI,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('mean dsi','FontSize',16);
sum(used)

sum(meanDSI(used)>=0.5)
sum(meanDSI(used)<0.5)

meanDSI(meanDSI>1.05)=1.05;
meanDSI(meanDSI<0)=0;
figure
hist(meanDSI(used),[0.05:0.1:1.05]);
axis([0 1.1 0 50])
figure
plot(meanDSI(used),ds_est(used),'*')



%%%%%%%%%%%% linearity %%%%%%%%
used = driftpeak>2;

layer_scatter(driftF1F0,used,layer,exc,inh);
sum(used)
xlabel('linearity (F1/F0)','FontSize',16)

layer_bar(driftF1F0,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('linearity (F1/F0)','FontSize',16);
sum(used)

layer_bar_inh_median(driftF1F0,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('linearity (F1/F0)','FontSize',16);
sum(used)

figure
hist(driftF1F0(used),0.05:.1:1.05)

figure
plot(simple_est(used),driftF1F0(used),'*')

[p h] =ranksum(driftF1F0(used&(layertwothree|layerfour)&exc),driftF1F0(used&layersix&exc))
[p h] =ranksum(driftF1F0(used&(layertwothree)&exc),driftF1F0(used&layerfour&exc))
%%%% spatial frequency

used = driftpeak>1.5  ;
layer_scatter(logSF,used,layer,exc,inh);
sum(used)
xlabel('pref spatial freq','FontSize',16)

logSF = driftSF;
logSF(logSF<.005)=.005;
logSF = log2(logSF/.01);
figure
barweb([mean(logSF(used&exc&layertwo)), mean(logSF(used&exc&layerthree)), mean(logSF(used&exc&layerfour)) ,mean(logSF(used&exc&layerfive)),mean(logSF(used&exc&layersix)); 0 0 0 0 0]', ...
    [std(logSF(used&exc&layertwo))/sqrt(sum(used&exc&layertwo )), std(logSF(used&exc&layerthree))/sqrt(sum(used&exc&layerthree)) ,  ...
    std(logSF(used&exc&layerfour))/sqrt(sum(used&exc&layerfour )), std(logSF(used&exc&layerfive))/sqrt(sum(used&exc&layerfive )) , std(logSF(used&exc&layersix))/sqrt(sum(used&exc&layersix)); 0 0 0 0 0]', ...
    2,{'upper 2/3' 'lower 2/3' '4' '5' '6'},[],[],[],[0 0 .75]);
set(gca,'xticklabel', {'2', '3', '4', '5' '6'})
ylabel('pref SF');
[h p] = ranksum(logSF(used&exc&~layersix),logSF(used&exc&(layersix)))
[h p] = ranksum(logSF(used&inh),logSF(used&exc&~layersix))

used = driftpeak>2;
layer_bar_inh(logSF,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('pref spatial freq','FontSize',16);
sum(used)


figure
hist(log2(driftSF(used)/.01),[-1:.5:6])
title('pref SF')

% figure
% used = driftpeak>2 & ~layerfive & exc;
% plot(logSF,driftF1F0,'o')



driftSF_width = driftSF_width*1.17;
driftSF_width = abs(driftSF_width);
driftSF_width(driftSF_width<0.35)=0.35;
driftSF_width(driftSF_width>3)=3;

logSF=driftSF;
logSF(logSF<.005)=.005;
logSF= log2(logSF/.01);
used = sb_x>0 & driftpeak>2;
figure
plot(sqrt((sb_x(used)-monitor_offset(used)).^2 + (sb_y(used)-36).^2),logSF(used),'o');
xlabel('azimuth');
ylabel('SF');
corrcoef(sqrt((sb_x(used)-monitor_offset(used)).^2 + (sb_y(used)-36).^2),logSF(used))

[h p] =corrcoef(sqrt((sb_x(used)-monitor_offset(used)).^2 + (sb_y(used)-36).^2),logSF(used))

%%%% sf width
used = driftpeak>1.5;
figure
layer_scatter(driftSF_width,used,layer,exc,inh);
sum(used)
xlabel('SF width','FontSize',16)

layer_bar_inh_median(driftSF_width,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('Sf width','FontSize',16);
sum(used)

figure
hist((driftSF_width(used)),0:0.25:3)

median(driftSF_width(used))

%%%%%%%%spont

plot(log10(meanspont(inh)),layer(inh),'ro');
hold on
plot(log10(meanspont(exc)),layer(exc),'go');
title('spont')
axis ij


meanspont_floor = meanspont;
meanspont_floor(meanspont_floor==0) =.01;
figure
hist(log10(meanspont_floor))

used = ~isnan(meanspont);

logspont = log10(meanspont_floor);
[h p] =ttest2(logspont(used&exc&layertwothree),logspont(used&exc&layerfour))
[p h ] =ranksum(meanspont(used&exc&layertwothree),meanspont(used&exc&layerfour))

layer_scatter(logspont,used,layer,exc,inh);
sum(used)
xlabel('spont rate (spikes/sec)','FontSize',16)

used = ~isnan(meanspont);
layer_bar_inh_median(meanspont,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('spont rate (spikes/sec)','FontSize',16);
sum(used)

used = ~isnan(meanspont);
layer_bar_inh(meanspont,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('spont rate (spikes/sec)','FontSize',16);
sum(used)

figure
used = inh&driftpeak>1;
plot(logspont(used),meanOSI(used),'o');

%%% rf size


used = barpeak>6 & ~isnan(bar_rfwidth);

layer_scatter(bar_rfwidth,used,layer,exc,inh);
sum(used)
xlabel('rf size (bars)','FontSize',16)

layer_bar(bar_rfwidth,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('rf size (bars)','FontSize',16);
sum(used)

layer_bar_inh(bar_rfwidth,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('rf size (bars)','FontSize',16);
sum(used)

used = barpeak>6 & driftpeak>2 & ~layerfive;
figure
plot(bar_rfwidth(used),driftF1F0(used),'o')


used = wx>0;
sbw= (wx+wy)/2;

layer_scatter(sbw*1.17,used,layer,exc,inh);
sum(used)
xlabel('rf size (deg)','FontSize',16)

layer_bar(sbw*1.17,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('rf size (deg)','FontSize',16);
sum(used)

layer_bar_inh(sbw*1.17,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('rf size (deg)','FontSize',16);
sum(used)

% used = sbw>0 & driftpeak>2 & ~layerfive &exc;
% figure
% plot(sbw(used),driftF1F0(used),'o')
% 
% used = (sbw>0) & (barpeak>6) &~isnan(bar_rfwidth);
% figure
% plot(sbw(used),bar_rfwidth(used),'o')
% corrcoef(sbw(used),bar_rfwidth(used))

% 
% used = wx>0;
% figure
% plot(wx(used&~layerfive),wy(used&~layerfive),'o');
% hold on
% plot(wx(used&layerfive),wy(used&layerfive),'go');

% used = (flash_on~=0 & flash_on<1.5);
% figure
% plot(flash_on(used),layer(used),'o');
% title('flash on')


% 
% used = barpeak>3 ;
% figure
% plot(meanOSI(used&inh),bar_rfwidth(used&inh),'ro');
% hold on
% plot(meanOSI(used&exc),bar_rfwidth(used&exc),'go');

%%% temporal frequency
used = abs(cp_tf_tuning(:,1))>0 & max(cp_tf_tuning,[],2)>0;
[r tfmax] = max(cp_tf_tuning,[],2);
2^mean(tfmax(used&layerfour==1))/2
2^mean(tfmax(used&layertwothree==1))/2
2^mean(tfmax(used&layerfive))/2

for i = 1:size(cp_tf_tuning,1);
    if tfmax(i)>1 & tfmax(i)<4
         tfmaxfit(i) = tfmax(i) - 0.5*((cp_tf_tuning(i,tfmax(i)+1)-cp_tf_tuning(i,tfmax(i)-1))/(cp_tf_tuning(i,tfmax(i)+1)+cp_tf_tuning(i,tfmax(i)-1)-2*cp_tf_tuning(i,tfmax(i))))
        i
    else
        tfmaxfit(i) = tfmax(i);
    end
end

figure
hist(tfmaxfit(used),[1:.25:4])

tfmax = 2.^(tfmaxfit-1);
layer_scatter(tfmax,used,layer,exc,inh);
sum(used)
xlabel('peak temp freq (Hz)','FontSize',16)

layer_bar(tfmax,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('peak temp freq (Hz)','FontSize',16);
sum(used)

layer_bar_inh_median(tfmax,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('peak temp freq (Hz)','FontSize',16);
sum(used)


[h p] =ranksum(tfmax(used&layerfour&exc), tfmax(used&layertwothree&exc))

tfmax = 2.^(tfmax-1);

used = abs(cp_tf_tuning(:,1))>0 & max(cp_tf_tuning,[],2)>0;
cp_tf_norm = cp_tf_tuning;
cp_tf_norm(cp_tf_norm<0)=0;
for i=1:size(cp_tf_norm,1)
    if used(i)==1
%      if min(cp_tf_norm(i,:))<0
%          cp_tf_norm(i,:) = cp_tf_norm(i,:) - min(cp_tf_norm(i,:));
%      end
     cp_tf_norm(i,:) = cp_tf_norm(i,:)/(cp_tf_norm(i,1));
    end
end
cp_tf_norm = cp_tf_tuning;
figure
%errorbar(mean(cp_tf_norm(find(layerfive&used),:)),std(cp_tf_norm(find(layerfive&used),:))/sqrt(sum(layerfive&used)),'r');
hold on
errorbar(mean(cp_tf_norm(find(layerfour&used),:)),std(cp_tf_norm(find(layerfour&used),:))/sqrt(sum(layerfour&used)),'g');
errorbar(mean(cp_tf_norm(find(layertwothree&used),:)),std(cp_tf_norm(find(layertwothree&used),:))/sqrt(sum(layertwothree&used)),'b');
legend('layer 4', 'layer 2/3')

[h p ] =ttest2(cp_tf_norm(find(layerfour&used),:),cp_tf_norm(find(layertwothree&used),:))



%%% movies
c_mod = contrast_mod; 
figure
used = movieN>200 ;
polar(contrast_phase(used&exc&~layerfive),(c_mod(used&exc&~layerfive)),'go')
hold on
polar(contrast_phase(used&exc&layerfive),(c_mod(used&exc&layerfive)),'ko')
polar(contrast_phase(used&inh),(c_mod(used&inh)),'ro')
polar(contrast_phase(used&exc&layersix),(c_mod(used&exc&layersix)),'mo')

figure
polar(-contrast_phase(used),c_mod(used),'o')

used = movieN>200 ;
layer_bar_median(contrast_mod,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('movie modulation','FontSize',16);
sum(used)


layer_bar_inh(contrast_mod,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('movie modulation','FontSize',16);
sum(used)

[h p] =ranksum(contrast_mod(used &exc &layerfive),contrast_mod(used &exc &layerthree))

figure
plot(responsive_est(used),c_mod(used).*cos(contrast_phase(used)),'o')

figure
hist(sta_halfcontrast(used),[0:.05:1])
used = movieN>200;
figure
hist(sta_adaptation(used),[-.25:.05:1])


layer_bar(sta_adaptation,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('movie adaptation','FontSize',16);
sum(used)

layer_bar(sta_halfcontrast,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('movie contrast','FontSize',16);
sum(used)


figure
hist(movieN(layerfour | layertwothree),[0:100:2000])
mean(meanspont((layerfour | layertwothree) & exc&~isnan(meanspont)))
mean(movieN((layerfour | layertwothree) & exc&~isnan(movieN))/600)


figure
hist(c_mod(used))

used =movieN>1;
c_mod = contrast_mod.*movieN/600
figure
plot(responsive_est(used),(c_mod(used)),'o')

for l = 2:6;
    used = movieN>200 &layer==l & responsive_est==1&exc;
    sum(contrast_mod(used)<0.3)/(sum(used))
end
    

used =movieN>200;
c_mod = contrast_mod.*movieN/600
figure
plot(responsive_est(used),(contrast_mod(used)),'o')

sta_freq(sta_freq<.02)=.02;
driftSF(driftSF<.02) = .02;
used = (sta_found==1 & movieN>300 & driftpeak>3 &driftSF<.12 );
figure
plot(log2(sta_freq(used)*(70/60)/.01),log2(driftSF(used)/.01),'o') %%% includes correction because display is not 70deg, as calculated in compile_analysis
hold on
plot([0.5 4],[0.5 4])

[r p]=corrcoef(log2(sta_freq(used)*(70/60)/.01),log2(driftSF(used)/.01))
r.^2

figure
plot(sta_theta(used),meantheta(used),'o')

dt = sta_theta-meantheta;

sta_theta(dt>90) = sta_theta(dt>90)-180;
sta_theta(dt<-90) = sta_theta(dt<-90)+180;

meantheta(sta_theta<0)=meantheta(sta_theta<0)+180;
sta_theta(sta_theta<0)=sta_theta(sta_theta<0)+180;

used = (sta_found==1 & movieN>200 & driftpeak>2 );
figure
plot(sta_theta(used),meantheta(used),'o')
axis([-20 200 -20 200]);
hold on
plot([-20 200],[-20 200]);
axis square

r = corrcoef(sta_theta(used),meantheta(used))
r.^2

used =  sta_found==1 &driftpeak>3;
figure
plot(driftF1F0(used),log10(movieN(used)),'o');
hold on
used =  sta_found<1 &driftpeak>3;
plot(driftF1F0(used),log10(movieN(used)),'ro');

figure
used = driftF1F0>0.5;
figure
hist(movieN(used&sta_found==1));
hold on
hist(movieN(used&sta_found==0));

used = sta_found==1;
layer_bar(sta_subunits,used,exc,inh,layertwo,layerthree,layerfour,layerfive,layersix);
ylabel('# subunits','FontSize',16);
sum(used)


used = sta_found==1 & ~isnan(meantheta_width);
figure
plot(sta_subunits(used),meantheta_width(used),'o')
title('theta width vs # subunits');

used = sta_found==1 & driftpeak>2;
figure
plot(sta_subunits(used),abs(driftSF_width(used)),'o')
title('sf width vs # subunits');


used = sta_found==1 & ~isnan(meantheta_width);
figure
barweb([mean(meantheta_width(used(1:length(sta_subunits)) & sta_subunits==1)) mean(meantheta_width(used(1:length(sta_subunits))&sta_subunits==2)) mean(meantheta_width(used(1:length(sta_subunits))&sta_subunits==3)); 0 0 0]', ...
       [std(meantheta_width(used(1:length(sta_subunits)) & sta_subunits==1))/3.4 std(meantheta_width(used(1:length(sta_subunits))&sta_subunits==2))/6.2 std(meantheta_width(used(1:length(sta_subunits))&sta_subunits==3))/4.2; 0 0 0]', 2);

used = sta_found==1 & driftpeak>2;
figure
barweb([mean(driftSF_width(used(1:length(sta_subunits)) & sta_subunits==1)) mean(driftSF_width(used(1:length(sta_subunits))&sta_subunits==2)) mean(driftSF_width(used(1:length(sta_subunits))&sta_subunits==3)); 0 0 0]', ...
       [std(driftSF_width(used(1:length(sta_subunits)) & sta_subunits==1))/3.4 std(driftSF_width(used(1:length(sta_subunits))&sta_subunits==2))/6.2 std(driftSF_width(used(1:length(sta_subunits))&sta_subunits==3))/4.2; 0 0 0]', 2);

figure
bar([mean(meanOSI(used(1:length(sta_subunits)) & sta_subunits==1)) mean(meanOSI(used(1:length(sta_subunits))&sta_subunits==2)) mean(meanOSI(used(1:length(sta_subunits))&sta_subunits==3))]);

figure
substwo = [sum(layertwo(1:length(sta_subunits)) & sta_subunits==1) sum(layertwo(1:length(sta_subunits)) & sta_subunits==2) sum(layertwo(1:length(sta_subunits)) & sta_subunits==3) ]/sum(layertwo(1:length(sta_subunits)) & sta_subunits>0);
substhree= [sum(layerthree(1:length(sta_subunits)) & sta_subunits==1) sum(layerthree(1:length(sta_subunits)) & sta_subunits==2) sum(layerthree(1:length(sta_subunits)) & sta_subunits==3) ]/sum(layerthree(1:length(sta_subunits)) & sta_subunits>0);
subsfour= [sum(layerfour(1:length(sta_subunits)) & sta_subunits==1) sum(layerfour(1:length(sta_subunits)) & sta_subunits==2) sum(layerfour(1:length(sta_subunits)) & sta_subunits==3) ]/sum(layerfour(1:length(sta_subunits)) & sta_subunits>0);
bar([substwo; substhree; subsfour]')

used= driftF1F0>0.6;

b1 = hist((movieN(used & sta_found>0)), [100:200:4100]);

b2 = hist((movieN(used & sta_found==0)), [100:200:4100]);
figure
bar([b1; b2]',1.5);


% 
% 
% 
% meanDSI = (barpeak.*barsDSI + driftpeak.*driftDSI)./(barpeak + driftpeak);
% meanDSI(meanDSI>1)=1;
% 
% 
% 

% layer = xlsread(fname,1,'B2:B242');
% 
% layer(layer==2.5 ) =2;
% layer(layer==3.5)=3;
% layer(layer==4.5)=5;
% layer(layer==5.5) = 5;
% layer(layer>6)=NaN;
% 
% layertwo = layer<3;
% layerthree= layer==3;
% layerfour = layer==4;
% layerfive = layer==5;
% layersix=layer==6;
% layertwothree = layertwo | layerthree;

simple_est = round(simple_est+rand(size(simple_est))*.2 -0.1);
orient_est = round(orient_est+rand(size(simple_est))*.2 -0.1);
responsive_est(responsive_est<1)=0;
%for i = 1:4
% 
% orient_est = double(meanOSIcirc>0.2) ;
% orient_est(isnan(meanOSI)==1)=NaN;

for i=5:5  
    figure
 switch i
        case 1
            used = (exc==1) & layertwo
            title('layer 2')
     case 2  
            used = exc==1 & layertwothree
            title('layer 3')
     case 3
            used = exc==1 & layerfour;
            title('layer 4')
        case 4
            used = exc==1 &layerfive;
            title('layer 5');
      case 5
            used = exc==1 &layersix;
            title('layer 6');       
        case 6
            used = inh==1;
            title('inhibitory')
    end
    
    num_est = sum(used)

    nonresp = sum((responsive_est==0 ) & used)/num_est;
    nonclass = sum(used&responsive_est==1 & isnan(simple_est + orient_est))/num_est;
    simp_orient = sum(simple_est==1 & orient_est==1& used &responsive_est==1)/num_est;
    comp_orient = sum(simple_est==0 & orient_est==1& used&responsive_est==1)/num_est;
    simp_nonorient = sum(simple_est==1 & orient_est==0& used&responsive_est==1)/num_est;
    comp_nonorient = sum(simple_est==0 & orient_est==0& used&responsive_est==1)/num_est;


%    pie([nonresp nonclass simp_orient comp_orient simp_nonorient comp_nonorient ]+10^-8);
%     legend( {'non-responsive' 'non-classified' 'linear oriented' 'nonlinear oriented' 'linear non-oriented' 'nonlinear non-oriented'});

  pie([nonresp+nonclass simp_orient comp_orient simp_nonorient comp_nonorient ]+10^-8);
    legend( {'non-responsive' 'linear oriented' 'nonlinear oriented' 'linear non-oriented' 'nonlinear non-oriented'});

%     figure
%     pie([nonresp simp_orient comp_orient simp_nonorient comp_nonorient],{'nonresponsive' 'linear oriented' 'nonlinear oriented' 'linear non-oriented' 'nonlinear non-oriented'} );

end


used = layer<=6;
   num_est = sum(used)

    nonresp = sum((responsive_est==0 ) & used)/num_est;
    nonclass = sum(used&responsive_est==1 & isnan(simple_est + orient_est))/num_est;
    simp_orient = sum(simple_est==1 & orient_est==1& used &responsive_est==1)/num_est;
    comp_orient = sum(simple_est==0 & orient_est==1& used&responsive_est==1)/num_est;
    simp_nonorient = sum(simple_est==1 & orient_est==0& used&responsive_est==1)/num_est;
    comp_nonorient = sum(simple_est==0 & orient_est==0& used&responsive_est==1)/num_est;