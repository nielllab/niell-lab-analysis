afile = {...
%    'C:\data\ephys matlab data\021412_awake_chr2\wndrift1\analysis.mat', ...     %%LFP problems
%     'C:\data\ephys matlab data\021512_awake_pptg\wn1b\analysis.mat', ...        %%LFP problems
    'D:\Jen_ephys_data\FragileX\1_23_13\Record_1_800uM\analysis_1.mat',...
             'D:\Jen_ephys_data\FragileX\1_23_13\Record_2_930uM\analysis_drift_2 and wn_2_1_23_13.mat',...
             'D:\Jen_ephys_data\FragileX\01_24_13\B\analysis_01_24_13_cluster_drift_1_B and wn_1_B.mat',...
             'D:\Jen_ephys_data\FragileX\1_11_13\record_1\analysis_1.mat',...
           'D:\Jen_ephys_data\FragileX\1_11_13\record_2\analysis_2.mat',...
           'D:\Jen_ephys_data\FragileX\1_11_13\record_3\analysis_3.mat',...
           'D:\Jen_ephys_data\FragileX\1_14_13\record_1\analysis_1.mat',...
          'D:\Jen_ephys_data\FragileX\1_14_13\record_2\analysis_2.mat'}

lfp_channel = [  ];

n=0
frame_duration = 1/30;

for i = 1:length(afile)                             %%%Loop to go through each exp (i-th exp)
    load(afile{i});
    cell_range = n+(1:size(cells,1));            %%%Creates index number for each unit in 
    n=cell_range(end);                           %%%the i-th experiment so all data can be in one vector
    
    %Unit Data
    site(cell_range)=i;                          %%Pulls out experiment num-id for each unit in i-th exp
    cell_id(cell_range,:) = cells;               %%Pulls out tetrode id and unit-id for each unit in it-th exp
    wvform(cell_range,:) = wv';                  %%Pulls out waveform for each unit in i-th exp
    peaksite(cell_range)=peakchan;
    
%     %Speed Data
%     laser_speed_all(i,:)=laserspeed;
%     laser_speed_all_std(i,:)=laserspeed_std;
%     
    %Supra/Infra
%     display(afile{i})  
%     LayerFour = input('Which channel is Layer 4? : ');
%     infra(cell_range)= (cell_id(cell_range,1)>LayerFour);
 
    
    for c=cell_range                                %%%Loop to go through each unit (c) in i-th exp
        ch=ceil(peaksite(c)/4)*4 -3;                %%%Pulls out channel/tetrode for LFPs
        laser_lfp_all{c} = laserlfp(ch,:,:);
        laser_lfp_freqs{c}=freqs{ch};
        
        
        for rep=1:3                         %%Types of Comparisons
            if rep==1
                data=Rcyclerep;              %%% laser off/on
            elseif rep==2
                data = movRcyclerep;         %%% stationary vs moving
            elseif rep==3
                data=statlaserRcyclerep;     %%% laser off/on, but only stationary
            end          
            
            for state=1:2
                d = condenseData(data{find(cell_range==c),state}',15)/frame_duration;
                d = (d(1:10) + d(20:-1:11))/2;               
                spont(c,rep,state) = mean(d(1:2));
                evoked(c,rep,state)= mean(d(8:10))-mean(d(1:2));
                if state==1
                    fr1 = d(1:10);
                elseif state==2
                    fr2 = d(1:10);
                end
            end
            
            %%% put fit of Rm = alpha*(Rs-R0) + beta here
                Bs = polyfit(fr1,fr2,1);
                B1(c,rep) = Bs(2);
                B0(c,rep) = Bs(1);
            
        end
    end
end
    

clear gamma alpha

for i= 1:n

   lfp = squeeze(laser_lfp_all{i});
   f=laser_lfp_freqs{i};
   f = f(f<80);
   
   noiserange = (f>56 & f<63) | (f>49 &f<51) | (f>39 & f<41);
   
%    figure
%    hold on
%    plot(f( ~noiserange),lfp(1,~noiserange),'b');
%    plot(f(~noiserange),lfp(2,~noiserange),'r');
   
   gammaF = find(f>40 & f<70 & ~noiserange);
   alphaF = find(f>10 & f<30);
   gamma(i,:)= mean(lfp(:,gammaF),2);
   [peakgamma(i,:) freqs] = max(lfp(:,gammaF),[],2);
   peakgammaF(i,:) = f(gammaF(freqs));
   alpha(i,:) = mean(lfp(:,alphaF),2);
   title(sprintf('peak = %0.1f %0.1f; peakF = %0.1f %0.1f area = %0.1f %0.1f',peakgamma(i,1),peakgamma(i,2),peakgammaF(i,1),peakgammaF(i,2),gamma(i,1),gamma(i,2)))
end

figure
plot(peakgamma(:,1),peakgamma(:,2),'o');
title('peak value');

figure
plot(peakgammaF(:,1),peakgammaF(:,2),'o');
title('freq of peak value');

figure
plot(peakgamma(:,1),peakgamma(:,2),'o');
title('peak value');


figure
plot(gamma(:,1),gamma(:,2),'o');
hold on
title('gamma')
xlabel('laser off');
ylabel('laser on');
axis equal
plot([0 10^4],[0 10^4])

figure
plot(alpha(:,1),alpha(:,2),'o');
hold on
title('alpha')
xlabel('laser off');
ylabel('laser on');
axis equal
plot([0 10^4],[0 10^4])
axis equal


%plot for laser data

% figure
% plot(laser_speed_all(:,1),laser_speed_all(:,2),'o');
% hold on
% title('laser induced movement');
% xlabel('laser off cm/sec');
% ylabel('laser on cm/sec');
% axis equal
% axis square
% plot([0 10],[0 10])

%plot for laser data
% figure
% barwitherr([std(laser_speed_all(:,1)) std(laser_speed_all(:,2))]/sqrt(length(laser_speed_all)), ...
%    [ mean(laser_speed_all(:,1)) mean(laser_speed_all(:,2))]);
% legend({'laser off','laser on'});
% title_str = {'laser','movement','laser no move'};

used= find((evoked(:,2,2)>2));
usedInfra = find((evoked(:,2,2)>2) & (spont(:,2,1)>=2));                                                                    
usedSupra = find((evoked(:,2,2)>2) & (spont(:,2,1)<2));    

sprintf('used fraction = %d / %d  = %f',length(used),length(evoked),length(used)/length(evoked))    %Prints fraction of units used
sprintf('number of infra = %d',length(usedInfra))    
sprintf('number of supra = %d',length(usedSupra))    

for rep=1:3;
   
    figure
    plot(spont(used,rep,1),spont(used,rep,2),'ko');hold on;
    plot(spont(usedInfra,rep,1),spont(usedInfra,rep,2),'ro');hold on;
    plot(spont(usedSupra,rep,1),spont(usedSupra,rep,2),'bo');hold on;
    axis square; axis equal;hold on
    plot([0 max(spont(used))], [0 max(spont(used))]);hold on;
    title(sprintf('spont %s',title_str{rep}));
    
   figure
    plot(B0(used,rep),B1(used,rep),'ko');hold on;
    plot(B0(usedInfra,rep),B1(usedInfra,rep),'ro');hold on;
    plot(B0(usedSupra,rep),B1(usedSupra,rep),'bo');hold on;
    plot([1 1],[-6 10],'k-'); hold on;
    plot([.5 5],[0 0],'k-'); hold on;
%     plot([0 max(Betas(2)(used))], [0 max(Betas(2)(used))]);hold on;
    title(sprintf('spont B %s',title_str{rep}));
    
    figure
    plot(evoked(used,rep,1),evoked(used,rep,2),'ko');hold on;
    plot(evoked(usedInfra,rep,1),evoked(usedInfra,rep,2),'ro');hold on;
    plot(evoked(usedSupra,rep,1),evoked(usedSupra,rep,2),'bo');hold on;
    axis square; axis equal; hold on
    plot([0 max(evoked(used))], [0 max(evoked(used))]); hold on;
    title(sprintf('evoked %s',title_str{rep}));
  
    figure
    barwitherr( [std(spont(used,rep,1),1) std(spont(used,rep,2),1) ; std(evoked(used,rep,1),1) std(evoked(used,rep,2),1)]/sqrt(n),...
        [nanmean(spont(used,rep,1),1) nanmean(spont(used,rep,2),1) ; nanmean(evoked(used,rep,1),1) nanmean(evoked(used,rep,2),1)] ); 
    title(title_str{rep});
    [p t] = ranksum(spont(used,rep,1),spont(used,rep,2));
    [psr t] = signrank(spont(used,rep,1),spont(used,rep,2));
    sprintf('%s spont p = %f p(signrank) = %f',title_str{rep},p,psr)
     [p t] = ranksum(evoked(used,rep,1),evoked(used,rep,2));
     [psr t] = signrank(evoked(used,rep,1),evoked(used,rep,2));
    sprintf('%s evoked p = %f p(signrank) = %f',title_str{rep},p,psr)
end