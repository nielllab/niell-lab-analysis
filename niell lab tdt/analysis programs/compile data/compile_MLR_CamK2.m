clear all
afile = {...
%    'C:\data\ephys matlab data\021412_awake_chr2\wndrift1\analysis.mat', ...     %%LFP problems
%     'C:\data\ephys matlab data\021512_awake_pptg\wn1b\analysis.mat', ...        %%LFP problems
     'C:\data\ephys matlab data\021612_awake\wn1drift1\analysis.mat', ...         % MLRstim #1; great desynchronization
    'C:\data\ephys matlab data\021612_awake\wn2b\analysis.mat',...                % MLRstim #2; good desychronization
    'C:\data\ephys matlab data\021612_awake\wn3\analysis.mat',...                 % MLRstim #3; good desychronization at 50Hz
    'C:\data\ephys matlab data\021612_awake\wn4bdrift5\analysis.mat',...          % MLRstim #4    
    'C:\data\ephys matlab data\021812_awake_pptg\wn2\analysis.mat',...            % MLRstim #5    
    'C:\data\ephys matlab data\021812_awake_pptg\wn3e_drift2\analysis.mat',...    % MLRstim #6 % bad lfp
    'C:\data\ephys matlab data\021812_awake_pptg\wn5\analysis.mat',...            % MLRstim #7
    'C:\data\ephys matlab data\070112_awake_mlr\wn1fg\analysis.mat',...           % MLRstim #8  % bad lfp
    'C:\data\ephys matlab data\070112_awake_mlr\wn2a\analysis.mat',...            % MLRstim #9  % bad lfp
    'C:\data\ephys matlab data\070112_awake_mlr\wn3d\analysis.mat',...            % MLRstim #10
    'C:\data\ephys matlab data\070112_awake_mlr\wn4\analysis.mat',...             % MLRstim #11  %%% bad lfp
    'C:\data\ephys matlab data\070112_awake_mlr\wn5b\analysis.mat',...            % MLRstim #12
    'C:\data\ephys matlab data\070312_awake_mlr\wn1h\analysis.mat',...            % MLRstim #13
    'C:\data\ephys matlab data\070312_awake_mlr\wn2e\analysis.mat',...            % MLRstim #14
    'C:\data\ephys matlab data\070312_awake_mlr\wn3a\analysis.mat',...            % MLRstim #15
   'C:\data\ephys matlab data\070312_awake_mlr\wn4c\analysis.mat'...              % MLRstim #16
}
lfp_channel = [  ];

n=0
frame_duration = 1/30;

for i = 1:length(afile)                          %%%Loop to go through each exp (i-th exp)
%for i =1:8
    load(afile{i});
    cell_range = n+(1:size(cells,1));            %%%Creates index number for each unit in 
    n=cell_range(end);                           %%%the i-th experiment so all data can be in one vector
    
    %Unit Data
    site(cell_range)=i;                          %%Pulls out experiment num-id for each unit in i-th exp
    cell_id(cell_range,:) = cells;               %%Pulls out tetrode id and unit-id for each unit in it-th exp
    wvform(cell_range,:) = wv';                  %%Pulls out waveform for each unit in i-th exp
    peaksite(cell_range)=peakchan;
    
    %Speed Data
    laser_speed_all(i,:)=laserspeed;
    laser_speed_all_std(i,:)=laserspeed_std;
    
    %Supra/Infra
%     display(afile{i})  
%     LayerFour = input('Which channel is Layer 4? : ');
%     infra(cell_range)= (cell_id(cell_range,1)>LayerFour);
 
    
    for c=cell_range                                %%%Loop to go through each unit (c) in i-th exp
        ch=ceil(peaksite(c)/4)*4 -3;                %%%Pulls out channel/tetrode for LFPs
        laser_lfp_all{c} = laserlfp(ch,:,:);
        laser_lfp_freqs{c}=freqs{ch};
        
        
        for rep=1:5                       %%Types of Comparisons
            if rep==1
                data=Rcyclerep;              %%% laser off/on
            elseif rep==2
                data = movRcyclerep;         %%% stationary vs moving
            elseif rep==3
                data=statlaserRcyclerep;     %%% laser off/on, but only stationary
            elseif rep==4
            data = movnolaserRcyclerep;     %%% stationary vs moving, laser off
            elseif rep==5
                data=movinglaserRcyclerep;     %%% laser off/on, but only moving
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
badnoise = zeros(1,n);
close all
for i= 25:25
%for i = 83:93
   lfp = squeeze(laser_lfp_all{i});
   f=laser_lfp_freqs{i};
   f = f(f<80);
   
   noise(i) = max(lfp(:))
   if max(lfp(:))>10^6
       badnoise(i)=1;
   else 
       badnoise(i)=0;
   end
   
   figure
   plot(f,lfp(:,f<80)')
   ylim([0 10^4])
   noiserange = (f>57 & f<63) | (f>49 &f<51) | (f>39 & f<41);
   keyboard
   
%    figure
%    hold on
%    plot(f( ~noiserange),lfp(1,~noiserange),'b');
%    plot(f(~noiserange),lfp(2,~noiserange),'r');
   
   gammaF = find(f>45 & f<60 & ~noiserange); %%% was 40 - 70 
   alphaF = find(f>2 & f<12);
   gamma(i,:)= mean(lfp(:,gammaF),2);
   [peakgamma(i,:) freqs] = max(lfp(:,gammaF),[],2);
   peakgammaF(i,:) = f(gammaF(freqs));
   
   [peakalpha(i,:) freqs] = max(lfp(:,alphaF),[],2);
   peakalphaF(i,:) = f(alphaF(freqs));
   
   alpha(i,:) = mean(lfp(:,alphaF),2);
   title(sprintf('%d peak = %0.1f %0.1f; peakF = %0.1f %0.1f area = %0.1f %0.1f',i,peakgamma(i,1),peakgamma(i,2),peakgammaF(i,1),peakgammaF(i,2),gamma(i,1),gamma(i,2)))
end

figure
plot(peakgamma(~badnoise,1),peakgamma(~badnoise,2),'o');
title('peak value');

figure
plot(peakgammaF(:,1),peakgammaF(:,2),'o');
title('freq of peak value');


figure
plot(peakalpha(~badnoise,1),peakalpha(~badnoise,2),'o');
title('peak value');

figure
plot(peakalphaF(~badnoise,1),peakalphaF(~badnoise,2),'o');
title('freq of peak value');


sprintf('number of good lfp recordings %d',sum(~badnoise))
keyboard

clear normgamma normalpha
for rep =2:size(gamma,2);
    
%normgamma(:,rep-1) = gamma(:,rep)./gamma(:,1)
%normalpha(:,rep-1) = alpha(:,rep)./alpha(:,1)


 normgamma(:,rep-1) = peakgamma(:,rep)./peakgamma(:,1);
 normalpha(:,rep-1) = peakalpha(:,rep)./peakalpha(:,1);
end

meangamma = mean(normgamma,2);
 
normgamma = normgamma(~badnoise &~isnan(meangamma'),:);
normgamma = unique(normgamma,'rows');

normalpha = normalpha(~badnoise&~isnan(meangamma'),:);
normalpha = unique(normalpha,'rows');

figure
% barwitherr([0 nanstd(normalpha,[],1)./sqrt(sum(~isnan(normalpha))); 0 nanstd(normgamma,[],1)./sqrt(sum(~isnan(normgamma)))],...
%     [1 nanmedian(normalpha,1); 1 nanmedian(normgamma,1)]);

barwitherr([0 nanstd(bootstrp(1000,@(x) nanmedian(x,1),normalpha)); 0 nanstd(bootstrp(1000,@(x) nanmedian(x,1),normgamma))],...
    [1 nanmedian(normalpha,1); 1 nanmedian(normgamma,1)]);


legend({'laser off stat','laser on stat','laser off mov','laser on mov'});
set(gca,'xTickLabel',{'alpha','gamma'});

sprintf('lfp N = %d',length(normgamma))

gammaAll(:,2:4) = normgamma; gammaAll(:,1)=1;

alphaAll(:,2:4) = normalpha; alphaAll(:,1)=1;

for i= 1:4
    for j=1:4
        rankgamma(i,j) = ranksum(gammaAll(:,i),gammaAll(:,j));
        rankalpha(i,j) = ranksum(alphaAll(:,i),alphaAll(:,j));
    end
end

rankalpha
rankgamma
keyboard 
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
title_str = {'laser','movement','laser no move','move no laser','laser effect, moving'};

used= find((evoked(:,2,2)>1));
usedInfra = find((evoked(:,2,2)>2) & (spont(:,2,1)>=2));                                                                    
usedSupra = find((evoked(:,2,2)>2) & (spont(:,2,1)<2));    

sprintf('used fraction = %d / %d  = %f',length(used),length(evoked),length(used)/length(evoked))    %Prints fraction of units used
sprintf('number of infra = %d',length(usedInfra))    
sprintf('number of supra = %d',length(usedSupra))    

keyboard

display('mlr')
cond{1} = spont(used,3,1); cond{2} = spont(used,3,2); cond{3} = spont(used,4,2); cond{4} = spont(used,5,2);
for c = 1:4
    for d = 1:4;
        p(c,d) = signrank(cond{c},cond{d})*6;
    end
end
sprintf('p values spont')
p

cond{1} = evoked(used,3,1); cond{2} = evoked(used,3,2); cond{3} = evoked(used,4,2); cond{4} = evoked(used,5,2);
for c = 1:4
    for d = 1:4;
        p(c,d) = signrank(cond{c},cond{d})*6;
    end
end
sprintf('p values evoked')
p


%%% plot four conditions (no stim, stim, move, move+stim)
figure
barwitherr([nanstd(bootstrp(1000,@(x) nanmedian(x,1),spont(used,3,1))) nanstd(bootstrp(1000,@(x) nanmedian(x,1),spont(used,3,2))) nanstd(bootstrp(1000,@(x) nanmedian(x,1),spont(used,4,2))) nanstd(bootstrp(1000,@(x) nanmedian(x,1),spont(used,5,2))); ...
        nanstd(bootstrp(1000,@(x) nanmedian(x,1),evoked(used,3,1))) nanstd(bootstrp(1000,@(x) nanmedian(x,1),evoked(used,3,2))) nanstd(bootstrp(1000,@(x) nanmedian(x,1),evoked(used,4,2))) nanstd(bootstrp(1000,@(x) nanmedian(x,1),evoked(used,5,2)))], ...
        [nanmedian(spont(used,3,1)) nanmedian(spont(used,3,2)) nanmedian(spont(used,4,2)) nanmedian(spont(used,5,2)); ...
        nanmedian(evoked(used,3,1)) nanmedian(evoked(used,3,2)) nanmedian(evoked(used,4,2)) nanmedian(evoked(used,5,2))]);
        
    legend({'no stim','stim - no move','move - no stim','stim + move'})
    set(gca,'Xticklabel',{'spont','evoked'})
    keyboard
for rep=1:5;
   
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
    plot(mean(B0(used,rep)),mean(B1(used,rep)),'g+');hold on;
    plot(mean(B0(used,rep)),mean(B1(used,rep))+std(B1(used,rep)/sqrt(length(B1(used,rep)))),'g+');hold on;
    plot(mean(B0(used,rep)),mean(B1(used,rep))-std(B1(used,rep)/sqrt(length(B1(used,rep)))),'g+');hold on;
    plot(mean(B0(used,rep))+std(B0(used,rep)/sqrt(length(B0(used,rep)))),mean(B1(used,rep)),'g+');hold on;
    plot(mean(B0(used,rep))-std(B0(used,rep)/sqrt(length(B0(used,rep)))),mean(B1(used,rep)),'g+');hold on;
    plot([1 1],[-6 10],'k-'); hold on;
    plot([.5 5],[0 0],'k-'); hold on;
%     plot([0 max(Betas(2)(used))], [0 max(Betas(2)(used))]);hold on;
    title(sprintf('B %s',title_str{rep}));
    
    figure;
    offset=(B1(used,rep)+(B0(used,rep)-1).*spont(used,rep,1));
    plot(B0(used,rep),offset,'ko');hold on;
    plot([1 1],[-4 10],'k-'); hold on;
    plot([-1 3.5],[0 0],'k-'); hold on;
    plot(nanmedian(B0(used,rep)),nanmedian(offset),'go');hold on;
    plot([nanmedian(B0(used,rep)) nanmedian(B0(used,rep))],[nanmedian(offset)+nanstd(offset)/sqrt(length(offset)) nanmedian(offset)-nanstd(offset)/sqrt(length(offset))],'g');hold on;
    plot([nanmedian(B0(used,rep))+nanstd(B0(used,rep))/sqrt(length(offset)) nanmedian(B0(used,rep))-nanstd(B0(used,rep))/sqrt(length(offset))],[nanmedian(offset) nanmedian(offset)],'g');

%     plot([0 max(Betas(2)(used))], [0 max(Betas(2)(used))]);hold on;
    axis([-1 3.5 -4 10]);
    xlabel('Gain');
    ylabel('Offset (Hz)');
    title(sprintf('Gain vs Offset %s',title_str{rep}));
  
    figure
    plot(evoked(used,rep,1),evoked(used,rep,2),'ko');hold on;
    plot(evoked(usedInfra,rep,1),evoked(usedInfra,rep,2),'ro');hold on;
    plot(evoked(usedSupra,rep,1),evoked(usedSupra,rep,2),'bo');hold on;
    axis square; axis equal; hold on
    plot([0 max(evoked(used))], [0 max(evoked(used))]); hold on;
    title(sprintf('evoked %s',title_str{rep}));
  
    figure
%     barwitherr( [std(spont(used,rep,1),1) std(spont(used,rep,2),1) ; std(evoked(used,rep,1),1) std(evoked(used,rep,2),1)]/sqrt(n),...
%         [nanmedian(spont(used,rep,1),1) nanmedian(spont(used,rep,2),1) ; nanmedian(evoked(used,rep,1),1) nanmedian(evoked(used,rep,2),1)] ); 

    barwitherr([nanstd(bootstrp(1000,@(x) nanmedian(x,1),spont(used,rep,1))) nanstd(bootstrp(1000,@(x) nanmedian(x,1),spont(used,rep,2))) ; nanstd(bootstrp(1000,@(x) nanmedian(x,1),evoked(used,rep,1))) nanstd(bootstrp(1000,@(x) nanmedian(x,1),evoked(used,rep,2)))],...
        [nanmedian(spont(used,rep,1),1) nanmedian(spont(used,rep,2),1) ; nanmedian(evoked(used,rep,1),1) nanmedian(evoked(used,rep,2),1)] ); 
    
    
    title(title_str{rep});
    set(gca,'Xticklabel',{'spont','evoked'});
    
    [p t] = ranksum(spont(used,rep,1),spont(used,rep,2));
    [psr t] = signrank(spont(used,rep,1),spont(used,rep,2));
    sprintf('%s spont p = %f p(signrank) = %f',title_str{rep},p,psr)
     [p t] = ranksum(evoked(used,rep,1),evoked(used,rep,2));
     [psr t] = signrank(evoked(used,rep,1),evoked(used,rep,2));
    sprintf('%s evoked p = %f p(signrank) = %f',title_str{rep},p,psr)
end

figure;
offsets=(B1(used,:)+(B0(used,:)-1).*spont(used,:,1));
bar([nanmean(B0(used,3)) nanmean(B0(used,4)) nanmean(B0(used,5))]);
hold on;
errorbar([nanmean(B0(used,3)) nanmean(B0(used,4)) nanmean(B0(used,5))],[nanstd(B0(used,3))/sqrt(length(B0(used,3))) nanstd(B0(used,4))/sqrt(length(B0(used,3))) nanstd(B0(used,5))/sqrt(length(B0(used,3)))])

figure;
offsets=(B1(used,:)+(B0(used,:)-1).*spont(used,:,1));
bar([nanmean(offsets(:,3)) nanmean(offsets(:,4)) nanmean(offsets(:,5))]-1);
hold on;
errorbar([nanmean(offsets(:,3)) nanmean(offsets(:,4)) nanmean(offsets(:,5))],[nanstd(offsets(:,3))/sqrt(length(offsets(:,3))) nanstd(offsets(:,4))/sqrt(length(offsets(:,4))) nanstd(offsets(:,5))/sqrt(length(offsets(:,5)))])

figure;
offsets=(B1(used,:)+(B0(used,:)-1).*spont(used,:,1));
[bootstat1,bootsam1]=bootstrp(1000,@nanmedian,[B0(used,3) B0(used,4) B0(used,5)]-1);
se1=std(bootstat1);

bar([nanmean(B0(used,3)) nanmean(B0(used,4)) nanmean(B0(used,5))]-1);
hold on;
errorbar([nanmean(B0(used,3)) nanmean(B0(used,4)) nanmean(B0(used,5))]-1,se1)
signrank(B0(used,3)-1)
signrank(B0(used,4)-1)
signrank(B0(used,5)-1)

figure;
offsets=(B1(used,:)+(B0(used,:)-1).*spont(used,:,1));
[bootstat2,bootsam2]=bootstrp(1000,@nanmedian,[offsets(:,3) offsets(:,4) offsets(:,5)]);
se2=std(bootstat2);

bar([nanmedian(offsets(:,3)) nanmedian(offsets(:,4)) nanmedian(offsets(:,5))]);
hold on;
errorbar([nanmedian(offsets(:,3)) nanmedian(offsets(:,4)) nanmedian(offsets(:,5))],se2);
signrank(offsets(:,3))
signrank(offsets(:,4))
signrank(offsets(:,5))
