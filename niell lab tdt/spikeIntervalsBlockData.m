close all

nblocks = round(max(spikeT{1}/10^5))+1;

for block= 1:1

clear sp
    close all
for c = 1:length(blockdata.spikes);
    sp{c} = blockdata.spikes{c};
    
        spikespeed{c} = interp1(blockdata.mouseT,blockdata.mouseV,sp{c});
   
  % s=spikeT{c};
   %sp{c} = s(s>(block-1)*10^5 & s<(block-1 + 0.01)*10^5)-(block-1)*10^5;
   
   lfp_avg_range = -0.25:0.005:0.25;
   figure
   set(gcf,'Name',sprintf('cell %d',c))
   for rep = 1:3
   lfp_avg = zeros(size(lfp_avg_range));
   if rep==1
       mvspikes=sp{c}(spikespeed{c}<1);
   elseif rep==2
        mvspikes=sp{c}(spikespeed{c}>1);
   else
       mvspikes = (0.05 + 0.9*rand(size(mvspikes)))*max(mvspikes);
       
   end
   
   ns = min(length(mvspikes),1000);
   for s = 1:ns
       s
       lfp_avg=lfp_avg+interp1(blockdata.lfpT{c},blockdata.lfpV{c},mvspikes(s)+lfp_avg_range);
   end
      
   cmap='rgk'
   subplot(1,2,1)
 hold on
   plot(lfp_avg_range,lfp_avg/ns,cmap(rep)); 
   xlim([lfp_avg_range(1) lfp_avg_range(end)])
   subplot(1,2,2)
hold on
   plot(lfp_avg_range,fftshift(abs(fft(lfp_avg)))/ns,cmap(rep));
   xlim([lfp_avg_range(1) lfp_avg_range(end)])
   getframe(gcf)
   
   end
   
end

figure

for c = 1:length(sp)
   if ~isempty(sp{c})
       plot(sp{c}',c,'k.');
    hold on
    
    pre_isi = diff(sp{c}(1:end-1));
    post_isi = diff(sp{c}(2:end));
    b = find(pre_isi>0.1 & post_isi<0.01);
    burst_times{c} = sp{c}(b+1);
    if ~isempty(burst_times{c})
        plot(burst_times{c},c,'go');
    end
   end
end
xlim([0 120]);
ylim([0 c+1]);

hist_int = 0.004;
for c1 = 1:length(sp)
    spikehist{c1} = hist(sp{c1}(spikespeed{c1}<1),0:hist_int:1000);
    bursthist{c1} = hist(sp{c1}(spikespeed{c1}>1),0:hist_int:1000);
    spikehist{c1} = spikehist{c1}(1:end-1);
    bursthist{c1} = bursthist{c1}(1:end-1);
end


spike_fig = figure;
set(gcf,'Name',sprintf('spike block %d',block))

set(gcf,'Name','stationary')
burst_fig = figure;
set(gcf,'Name','moving')

clear xc burst_xc



n_cells = length(sp);
range = 0.08;
lag_range=-(range/hist_int):(range/hist_int);

xc = zeros(n_cells,n_cells,length(lag_range));
burst_xc=xc;
for c1=1:n_cells
    for c2 = c1:n_cells
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

end