afile = {'D:\Moses_ephys_data\083112_awake_BF\wn1a\analysis3.mat','D:\Moses_ephys_data\083112_awake_BF\wn2a\analysis.mat'}

close all
lfp_channel = [  ];
n=0
frame_duration = 1/30;

for i = 1:length(afile)
    load(afile{i});
    cell_range = n+(1:size(cells,1));
    n=cell_range(end);
    
    site(cell_range)=i;
    cell_id(cell_range,:) = cells;
    wvform(cell_range,:) = wv';
    peaksite(cell_range)=peakchan;
    
    laser_speed_all(i,:)=laserspeed;
    laser_speed_all_std(i,:)=laserspeed_std;
    
    for c=cell_range
        ch=ceil(peaksite(c)/4)*4 -3;
        laser_lfp_all{c} = laserlfp(ch,:,:);
        laser_lfp_freqs{c}=freqs{ch}
        
        
        for rep=1:3
            if rep==1
                data=Rcyclerep;  %%% laser off/on
            elseif rep==2
                data = movRcyclerep; %%%stationary vs moving
            elseif rep==3
                data=statlaserRcyclerep;    %%%% laser off/on, but only stationary
            end          
            
            for state=1:2
                d = condenseData(data{find(cell_range==c),state}',15)/frame_duration;
                d = (d(1:10) + d(20:-1:11))/2;               
                spont(c,rep,state) = mean(d(1:2));
                evoked(c,rep,state)=mean(d(8:10))-mean(d(1:2));
            end
        end
    end
end
    

clear gamma alpha

for i= 1:n

    lfp = squeeze(laser_lfp_all{i});
   f=laser_lfp_freqs{i};
   close
   figure
hold on
   plot(f(f<59),lfp(1,f<59),'b');
      plot(f(f<59),lfp(2,f<59),'r');
   gammaF = find(f>50 & f<58);
   alphaF = find(f>1 & f<8);
   gamma(i,:)= mean(lfp(:,gammaF),2);
   alpha(i,:) = mean(lfp(:,alphaF),2);
end

figure
plot(gamma(:,1),gamma(:,2),'o');
hold on
title('gamma')
xlabel('laser off');
ylabel('laser on');
axis equal
plot([0 10^3],[0 10^3])

figure
plot(alpha(:,1),alpha(:,2),'o');
hold on
title('alpha')
xlabel('laser off');
ylabel('laser on');
axis equal
plot([0 10^3],[0 10^3])
axis equal

figure
barwitherr([std(alpha(:,1))/sqrt(length(gamma)) std(alpha(:,2))/sqrt(length(gamma)); std(gamma(:,1))/sqrt(length(gamma)) std(gamma(:,2))/sqrt(length(gamma))], ...
    [mean(alpha(:,1)) mean(alpha(:,2)) ; mean(gamma(:,1)) mean(gamma(:,2))]);





figure
plot(laser_speed_all(:,1),laser_speed_all(:,2),'o');
hold on
title('laser induced movement');
xlabel('laser off cm/sec');
ylabel('laser on cm/sec');
axis equal
axis square
plot([0 10],[0 10])


figure
barwitherr([std(laser_speed_all(:,1)) std(laser_speed_all(:,2))]/sqrt(length(laser_speed_all)), ...
   [ mean(laser_speed_all(:,1)) mean(laser_speed_all(:,2))]);
legend({'laser off','laser on'});
title_str = {'laser','movement','laser no move'};


used = find(evoked(:,2,2)>2);
sprintf('used fraction = %d / %d  = %f',length(used),length(evoked),length(used)/length(evoked))

for rep=1:3;
   
    figure
    plot(spont(used,rep,1),spont(used,rep,2),'o');
    axis square; axis equal;hold on
    plot([0 max(spont(:))], [0 max(spont(:))])
    title(sprintf('spont %s',title_str{rep}));
    
    figure
    plot(evoked(used,rep,1),evoked(used,rep,2),'o');
     axis square; axis equal; hold on
     plot([0 max(evoked(:))], [0 max(evoked(:))])
    title(sprintf('evoked %s',title_str{rep}));
  
    figure
    barwitherr( [std(spont(used,rep,1),1) std(spont(used,rep,2),1) ; std(evoked(used,rep,1),1) std(evoked(used,rep,2),1)]/sqrt(n),...
        [nanmean(spont(used,rep,1),1) nanmean(spont(used,rep,2),1) ; nanmean(evoked(used,rep,1),1) nanmean(evoked(used,rep,2),1)] ); 
    title(title_str{rep});
end
