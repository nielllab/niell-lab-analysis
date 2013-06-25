
for c = 1:length(spikeT);
   % sp{c} = wn_movement(c).spikes;
   s=spikeT{c};
   sp{c} = s(s>2*10^5 & s<2.1*10^5)-2*10^5;
end

for c = 1:length(wn_movement)
    plot(sp{c}',c,'k.');
    hold on
    
    pre_isi = diff(sp{c}(1:end-1));
    post_isi = diff(sp{c}(2:end));
    b = find(pre_isi>0.1 & post_isi<0.01);
    burst_times{c} = sp{c}(b+1);
    plot(burst_times{c},c,'go');
    
end
xlim([0 120]);
ylim([0 c+1]);

for c1 = 1:length(wn_movement)
    spikehist{c1} = hist(sp{c1},0:0.001:1000);
    bursthist{c1} = hist(burst_times{c1},0:0.001:1000);
    spikehist{c1} = spikehist{c1}(1:end-1);
    bursthist{c1} = bursthist{c1}(1:end-1);
end


spike_fig = figure;
burst_fig = figure;

clear xc burst_xc

n_cells = length(sp);
lag_range=-50:50;
for c1=1:n_cells
    for c2 = 1:n_cells
        for lag =lag_range;
                      
            if c1==c2 & lag ==0
                xc(c1,c2,lag-lag_range(1)+1)=0;
                burst_xc(c1,c2,lag-lag_range(1)+1)=0;
            else
            xc(c1,c2,lag-lag_range(1)+1)=spikehist{c1}*circshift(spikehist{c2},[0 lag])';
            burst_xc(c1,c2,lag-lag_range(1)+1)=bursthist{c1}*circshift(bursthist{c2},[0 lag])';
            end
        end
        figure(spike_fig);
        subplot(n_cells,n_cells,c1+(c2-1)*n_cells);
        plot(squeeze(xc(c1,c2,:)));
        xlim([1 length(lag_range)])
       set(gca,'xtick',[]);set(gca,'ytick',[]);
        
          figure(burst_fig);
        subplot(n_cells,n_cells,c1+(c2-1)*n_cells);
        plot(squeeze(burst_xc(c1,c2,:)))  
          xlim([1 length(lag_range)])
       set(gca,'xtick',[]);set(gca,'ytick',[]);
    end
end