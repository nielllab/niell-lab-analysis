close all

%%% define inh
inh = alldata(:,5)<8.0 & alldata(:,6)<1.6;
mid = alldata(:,5)>8.5 & alldata(:,5)<10 ;
exc= alldata(:,5)>10 & alldata(:,6)>0 & alldata(:,6)<4 ;

figure
plot(alldata(inh,5),alldata(inh,6),'ro'); hold on

%%% show pinp psth by layer
figure
col = 'bgrcm';
for l=2:6
    plot(histbins,mean(psth_pinp(:,:),1),col(l-1));
    hold on
end
xlim([-30 45]);
plot([1 1], [ 0 8],'color',[0.3 0.3 0.3])
xlabel('msec')
legend('lyr 2','3','4','5','6')
title('avg pinp hist by layer')

figure
plot(histbins,mean(psth_pinp(:,:),1))



%%% calculate baseline, evoked, and artifact by choosing windows
baseline = mean(psth_pinp(:,5:45),2);
baseStd = std(psth_pinp(:,5:45),[],2);
artifact = psth_pinp(:,51);
ev = mean(psth_pinp(:,52:56),2);
evoked = ev-baseline;

zscore =evoked./baseStd;
zscore(zscore>10)=10;

%%% plot zscore vs evoked
figure
plot(zscore,evoked,'o')
hold on
plot(zscore(lyr==4),evoked(lyr==4),'go')
plot(zscore(lyr==5),evoked(lyr==5),'ro');
plot(zscore(lyr==6),evoked(lyr==6),'yo');
legend('lyr 2/3','4','5','6')
xlabel('zscore'); ylabel('evoked')

% artscore = (artifact - baseline)./baseStd;
% artscore(artscore>10)=10;
% 
% figure
% hist(artscore',0:10)

%%% choose pinped
pinped = (zscore>2& evoked>5); 

figure
plot(baseline,ev,'o')
axis equal
hold on
plot(baseline(pinped),ev(pinped),'go');
hold on
plot(baseline(pinped&inh),ev(pinped&inh),'gs');

legend('non','pinp')
plot([0 40],[0 40])
xlabel('baseline'); ylabel('laser')
sprintf('%d pinped neurons total',sum(pinped))

%%% define responsive
resp = driftpeak(:,2)>2 & driftpeak(:,1)>2;

%%% average psth for pinp vs non
figure
plot(histbins-1,median(psth_pinp(~pinped,:),1));
hold on
plot(histbins-1,median(psth_pinp(pinped,:),1),'g');
legend('non','pinp')
plot([0 0], [0 25],'Color',[0.5 0.5 0.5]);
xlim([-20 40]); xlabel('msec'); ylabel('sp/sec')
plot([-50 50],[median(baseline(pinped)) median(baseline(pinped))],'g:')
plot([-50 50],[median(baseline(~pinped)) median(baseline(~pinped))],'b:')

%%% # spikes evoked
nbins=8;
ev_spikes = mean(psth_pinp(pinped,54 + (1:nbins)),2)-baseline(pinped);
ev_spikes = nbins*ev_spikes/1000;
figure
hist(ev_spikes)
xlabel('# evoked spikes - pinped');

nbins=8;
ev_spikes = mean(psth_pinp(pinped,54 + (1:nbins)),2)-baseline(pinped);
ev_spikes = nbins*ev_spikes/1000;
h1 = hist(ev_spikes,-0.05:0.02:0.25)/length(ev_spikes);
ev_spikes = mean(~psth_pinp(~pinped,54 + (1:nbins)),2)-baseline(~pinped);
ev_spikes = nbins*ev_spikes/1000;
h2 = hist(ev_spikes,-0.05:0.02:0.25)/length(ev_spikes);
figure
bar((-0.05:0.02:0.25)+0.01,[h1; h2]')

%%% define wt
wt = GT'==3;

%%% fraction pinp by layer
figure
for l=2:6
    fraction(l) = sum(pinped & lyr==l & wt)/sum(lyr==l & wt);
end
bar(fraction)
xlabel('layer')
ylabel('fraction pinped')


%%% fraction responsive
frac_resp(1) = sum(resp& lyr==4 & GT'==2 & ~pinped)/sum( lyr==4 & GT'==2 & ~pinped);
frac_resp(2) = sum(resp& lyr==4 & GT'==1  & pinped)/sum( lyr==4 & GT'==1  & pinped);
figure
bar(frac_resp); ylim([0 1]); ylabel('fraction resp >2'); title('lyr 4')
set(gca,'xticklabel',{'non','pinp'})

peak_run_p=nanmedian(driftpeak(GT'==3 & lyr==4 & pinped & exc,2))
peak_stat_p=nanmedian(driftpeak(GT'==3 & lyr==4 & pinped & exc,1))
gain_ind_p=(peak_run_p-peak_stat_p)/(peak_run_p+peak_stat_p)

peak_run=nanmedian(driftpeak(GT'==3 & lyr==4 & ~pinped & exc,2))
peak_stat=nanmedian(driftpeak(GT'==3 & lyr==4 & ~pinped & exc,1))
gain_ind=(peak_run-peak_stat)/(peak_run+peak_stat)
%%% plot grating respons data

plotPinpData(driftspont,wt,lyr==0,pinped,inh,1)
plot([0 20],[0 20])
title('spont')


plotPinpData(driftF1F0,wt,lyr==0,pinped,inh,resp)
plot([0 2],[0 2])
title('F1F0')

plotPinpData(cvOSI,wt,lyr==0,pinped,inh,resp)
plot([0 1],[0 1])
title('cv OSI')


plotPinpData(cvDSI,wt,lyr<=4,pinped,inh,resp)
plot([0 1],[0 1])
title('cv DSI')

plotPinpData(driftpeak,GT'==1,lyr==4,pinped,inh,resp)
hold on
plot([0 50],[0 50])
title('peak')

%%% show psth of all pinped cells
%%% lots of figures but useful

% pinplist = find(pinped);
% for i=1:length(pinplist)
% c = pinplist(i);
% 
%     figure
% plot(psth_pinp(c,:)')
% hold on
% plot([0 100],[baseline(c) baseline(c)],'g')
% plot([0 100],[baseline(c)+baseStd(c)*3 baseline(c)+baseStd(c)*3],'r')
% title(sprintf('evoked = %f  zscore = %f layer = %d inh = %d gt = %d',evoked(c),zscore(c),lyr(c),inh(c),GT(c)))
% end


%%% waveform of pinped vs non
figure
plot(wvform(exc,:)','b');
hold on
plot(wvform(inh,:)','r'); hold on
plot(wvform(pinped,:)','g')
figure
plot(wvform(pinped,:)','g')
title('pinped')

I=find(inh)
k=find(pinped)
% for j=length(I)
% figure
% imagesc(STA_peak{1,12});
% end

h=ismember(k,I)

%generate STA for pinped cells
for j=1:length(STA_peak)

    if ismember(j,k)
figure
imagesc(STA_peak{1,j},[-64 64]);
    end
end

%generate STA for inhibitory cells
for j=1:length(STA_peak)

    if ismember(j,I)
figure
imagesc(STA_peak{1,j});
    end
end

figure
for j=1:length(STA_peak)

    if ismember(j,I)
plot(alldata(j,5),alldata(j,6),'ro'); hold on
    end
end

figure
subplot(1,2,1);
pie([sum(wt&pinped&~inh) sum(wt&pinped&inh)],{'broad','narrow'})
title('wt pinped')
subplot(1,2,2);
pie([sum(~wt&pinped&~inh) sum(~wt&pinped&inh)],{'broad','narrow'})
title('2A/2B ko pinped')

figure
plot(wvform(pinped,:)')

