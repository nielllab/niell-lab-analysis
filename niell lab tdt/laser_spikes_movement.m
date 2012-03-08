clear all
Tank_Name = '022811_pptg_linear'
Block_Name = 'bars16d1'

%%% read laser data
TTX = openTTX(Tank_Name,Block_Name); % to initialize
[laserT laserTTL] = read_laser(TTX);
invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');


smoothwindow_secs = 2;
dt = laserT(2)-laserT(1);
smoothwindow= ceil(smoothwindow_secs/dt)
lasersmooth = zeros(size(laserTTL));
for i = smoothwindow+1:length(laserTTL);
    lasersmooth(i)= mean(laserTTL(i-smoothwindow:i));
end
lasersmooth = lasersmooth/5;  %%% to calculate duty cycle (since ttl=5V)


%%% read mouse data
[tsamp vsmooth] = getBlockVelocity(Tank_Name,Block_Name,1);


lasersamp = interp1(laserT,lasersmooth,tsamp);


figure
plot(lasersamp,vsmooth,'o');

%%% read spike data
w =1
R = zeros(4,length(tsamp));
for ch = 1:4
    ch
    tic
    T = getSpikeTimes(Tank_Name,Block_Name,(ch-1)*4+1);
    multi_times{ch} = T;
    R(ch,:) = boxcarRate(tsamp, T, w);
    toc
end;


figure
plot(R');
for i = 1:4
    figure
    plot(vsmooth,R(i,:),'o');
    figure
    plot(lasersamp,R(i,:),'o');
end

figure
plot(laserT,lasersmooth*10)
hold on
plot(tsamp,vsmooth,'g')


[fname, pname] = uigetfile('*.mat','analysis data')
load(fullfile(pname,fname));


% [fname, pname] = uiputfile('*.mat','spike laser info')
% save(fullfile(pname,fname),'tsamp','vsmooth','cluster_times','laserT','laserTTL','cells')


for i = 1:length(cells);
    w =1
    Rsingle(i,:) = boxcarRate(tsamp,cluster_times{i},w);
    
    figure
    subplot(1,3,1)
    plot(vsmooth,Rsingle(i,:),'o');
    hold on
    [xbin ybin] = binplot(vsmooth,Rsingle(i,:),[-1:2:max(vsmooth)]);
    plot(xbin,ybin,'g','Linewidth',2);
    title(sprintf('ch:%d cl:%d',cells(i,1),cells(i,2)));
    xlabel('velocity (cm/sec)');
    ylabel('Firing rate (sp/sec)');

%     figure
%     plot(lasersamp,Rsingle(i,:),'o');
%     hold on
%     [xbin ybin] = binplot(lasersamp,Rsingle(i,:),[-.001:0.002:max(lasersamp)]);
%     plot(xbin,ybin,'g','Linewidth',2);
    
    %figure
    subplot(1,3,2)
    plot(tsamp,vsmooth,'g');
    hold on
    plot(tsamp,Rsingle(i,:));
    title(sprintf('ch:%d cl:%d',cells(i,1),cells(i,2)));
    legend('velocity','firing rate')
    xlabel('secs')
    ylabel('cm/sec  sp/sec')
    
    histbins = -0.5:0.01:0.5;
    edges = laserT(find(diff(laserTTL)>2.5));
    psth=zeros(length(cells),length(hist(0,histbins)));    
    for t = 1:length(edges);
        tdiff = squeeze(cluster_times{i})-edges(t);
        tdiff = tdiff(abs(tdiff)<max(histbins));
        psth(i,:) = psth(i,:) + hist(tdiff,histbins)/(length(edges)*(histbins(2)-histbins(1)));
    end
    
   % figure
    subplot(1,3,3)
    [n histcent] = hist(tdiff,histbins)
    plot(histcent*1000,squeeze(psth(i,:)));
    title(sprintf('laser triggered psth ch:%d cl:%d',cells(i,1),cells(i,2)));
      xlabel('msec') 
      ylabel('Rate (sp/sec)')
end



histbins = -0.5:0.01:0.5;
edges = laserT(find(diff(laserTTL)>2.5));
psth=zeros(length(cells),length(hist(0,histbins)));
for c=1:length(cells);
    for t = 1:length(edges);
        tdiff = squeeze(cluster_times{c})-edges(t);
        tdiff = tdiff(abs(tdiff)<max(histbins));
        psth(c,:) = psth(c,:) + hist(tdiff,histbins);
    end
end

figure
plot(psth')