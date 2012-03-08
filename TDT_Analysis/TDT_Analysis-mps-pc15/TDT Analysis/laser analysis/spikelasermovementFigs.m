for i = 1:length(cells);
    w =1
    Rsingle(i,:) = boxcarRate(tsamp,cluster_times{i},w);
    
    figure
    plot(vsmooth,Rsingle(i,:),'o');
    hold on
    [xbin ybin] = binplot(vsmooth,Rsingle(i,:),[-1:2:max(vsmooth)]);
    plot(xbin,ybin,'g','Linewidth',2);
    title(sprintf('ch:%d cl:%d',cells(i,1),cells(i,2)));
    xlabel('velocity (cm/sec)');
    ylabel('Firing rate (sp/sec)');

%%% plot laser vs firing rate
%%% not as useful as laser-triggered psth (below)
%     figure
%     plot(lasersamp,Rsingle(i,:),'o');
%     hold on
%     [xbin ybin] = binplot(lasersamp,Rsingle(i,:),[-.001:0.002:max(lasersamp)]);
%     plot(xbin,ybin,'g','Linewidth',2);
    
    figure
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
    
    figure
    [n histcent] = hist(tdiff,histbins)
    plot(histcent*1000,squeeze(psth(i,:)));
    title(sprintf('laser triggered psth ch:%d cl:%d',cells(i,1),cells(i,2)));
      xlabel('msec') 
      ylabel('Rate (sp/sec)')
end
