%function compile_analysis_movement
clear all
cells=0;
afiles = {'C:\data\matlab data\data 2009\021609_awake\allstim1_new\analysis.mat' ...
    'C:\data\matlab data\data 2009\021609_awake\allstim2_new\analysis2.mat' ...
    'C:\data\matlab data\data 2009\030609_awake_tet_analysis\allstim1\analysis.mat', ...
    'C:\data\matlab data\data 2009\022209_awake_linear_analysis\allstim3\analysis.mat', ....
    'C:\data\matlab data\data 2009\022309_awake_tet_analysis\allstim1\analysis.mat' ...
       'C:\data\matlab data\data 2009\022309_awake_tet_analysis\allstim2\analysis.mat'};
N =0;
for i = 1:length(afiles)
    
    i
    load(afiles{i});
    n_units = length(L_ratio);
    cellrange = N+1:N+n_units;
    N=N+n_units;

  
    alldata( cellrange,1:2) = cells;
    alldata( cellrange,3) = L_ratio;

    %%% waveform
    alldata( cellrange,4) = trough_width;
    alldata( cellrange,5) = trough2peak;
    alldata( cellrange,6) = -trough_depth./peak_height;

    %     for c = 1:n_units;
    %         ch = alldata(c,1);
    %         cl = alldata(c,2);
    %         min_t = min(mean_wvform (:,ch : ch+3,cl),[],1);
    %         [trough_allc trig_chan] = min(min_t);
    %         trig_chan = ch+trig_chan-1;
    %         alldata(c,7) = mean_wvform(size(mean_wvform,1),trig_chan,cl)/peak_height(c);
    %
    %         t1= squeeze(event_times_all(ch,find(idx_all(ch,:) == cl)));
    %         dt = diff(t1);
    %         dt =dt(dt<.02);
    %         n=hist(dt,.001:0.002:.02);
    %         [y alldata(c,8)] = max(n);
    %         n=hist(dt,.0005:0.001:.02);
    %         alldata(c,9) = max(n(3:8))./mean(n(15:20))
    %     end;


%         A1(cellrange,:)=bars_A1;
%         A2(cellrange,:)=bars_A2;
%         w(cellrange,:)=bars_w;
%         theta(cellrange,:)=bars_theta;
%         bspont(cellrange,:)=barspont;
%        
%         size(rf_width)
%      
%         rfw(cellrange,:) = rf_width*30;


size(driftorient_A1)
driftA1(cellrange,:) =driftorient_A1;
driftA2(cellrange,:) = driftorient_A2;
driftB(cellrange,:) = driftorient_B;
drift_theta_w(cellrange,:)=driftorient_thetawidth;
driftspont(cellrange,:) = drift_spont;

driftwpref(cellrange,:) = wpref_dog;
driftwbw(cellrange,:) = wbw_dog;


%size(wvform)
wvform(cellrange,:) = wv';

end

figure
bar(driftspont(1:29,:));

clear OSI width peak



for i = 1:N
    for j=1:2
        driftA1(i,j);
        driftA2(i,j);
        driftB(i,j);
        drift_theta_w(i,j);
    [OSI(i,j) width(i,j) peak(i,j)] = calculate_tuning(driftA1(i,j),driftA2(i,j),driftB(i,j),drift_theta_w(i,j))
    end
end

order = [18 28 34 10 16 17 1   19 30 33 44 2 3 4 5 6 7 8 9 11 12 13 14 15 20 21 22 23 24 25 26 27 ...
        29 31 32 35 36 37 38 39 40 41 42 43 45 46 47 48 49]
figure
bar(peak(order,:),1)

figure
bar(driftspont(order,:),1)


% figure
% plot(wvform(1:29,:)','g')
% hold on
% plot(wvform([2 9 15],:)','r');
% plot(wvform([11 14 25],:)','b');

inh = zeros(29,1)
inh([ 2 9 11 14 15 25])=1;
inh = zeros(29,1)
inh([ 2 9  15 ])=1;
figure
plot(wvform(inh==1,:)','r')
hold on
plot(wvform(inh==0,:)','g')

%peak(peak<0)=0;
order = [2 9 15 11 14 25 1 3 4 5 6 7 8 19 12 13 16 17 18 19 20 21 22 23 24 26 27 28 29];
figure
bar(peak(order,:))

%sp = (driftspont+bspont)/2;
sp = driftspont;
figure
bar(sp(order,:))

figure
plot(sp(order,2)-sp(order,1),peak(order,2)-peak(order,1),'o')

figure
barweb(median(driftspont(order(12:49),:),1),std(driftspont(order(12:49),:),[],1)/sqrt(37));
legend('stationary','running')

figure
barweb(mean(peak(order(12:49),:),1), std(peak(order(12:49),:),[],1)/sqrt(37));
legend('stationary','running')

data = [median(sp(order(12:49),:),1); median(peak(order(12:49),:),1)]
err = [std(sp(order(12:49),:),[],1)/sqrt(23) ; std(peak(order(12:49),:),[],1)/sqrt(37)]
figure
barweb(data,err)
legend('stationary','running')

figure
bar(driftspont(order,:))

figure
bar(OSI(order,:))

used=peak(:,1)>3 & peak(:,2)>2;
figure
bar(OSI(used,:))

used=peak(:,1)>2 & peak(:,2)>2 & OSI(:,1)>0.4 & OSI(:,2)>0.4;
figure
bar(width(used,:)*180/pi)

used=peak(:,1)>2 & peak(:,2)>2 & (driftwbw(1:49,1))>0 & (driftwbw(1:49,2))>0;
figure
bar(driftwbw(used,:))



