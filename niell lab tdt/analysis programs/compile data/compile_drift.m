%function compile developmental data
clear all
close all
dbstop if error
% dbclear all % this will exit out of db mode.
%[fname pname] =uiputfile('*.ps','pdf output'); psfilename=fullfile(pname,fname);  %%% get ps filename
psfilename = 'c:/test.ps';   %%% default location
if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file


N =0; cells=0;


afiles = {'C:\data\ephys matlab data\021412_awake_chr2\wndrift1\analysis.mat' ...
    'C:\data\ephys matlab data\021812_awake_pptg\wn3e_drift2\analysis.mat' ...
    'C:\data\ephys matlab data\021812_awake_pptg\drift3\analysis.mat' ...
    'C:\data\ephys matlab data\021612_awake\drift3b\analysis.mat'...
    'C:\data\ephys matlab data\021612_awake\drift4\analysis.mat'...
    'C:\data\ephys matlab data\021612_awake\drift1\analysis.mat'...
    'C:\data\ephys matlab data\021612_awake\wn4bdrift5\analysis.mat'}  %%% more locomotion with laser off?

%     'C:\data\ephys matlab data\070112_awake_mlr\drift5\analysis.mat'
% 'C:\data\ephys matlab data\070612_awake_mlr\drift7a\analysis.mat' %%basal
% forebrain stimulation

for i = 1:length(afiles)
    
    
    clear wn wn_movement
    clear LFP_movement
    
    load(afiles{i});
    n_units = length(L_ratio);
    cellrange = N+1:N+n_units;
    N=N+n_units;
    
    
    
    driftA1(cellrange,1:2)= field2array(drift,'A1');
    driftA2(cellrange,1:2)=field2array(drift,'A2');
    driftB(cellrange,1:2)= field2array(drift,'B');
    drift_theta_w(cellrange,1:2)=field2array(drift,'thetawidth');
    drift_theta(cellrange,1:2)=field2array(drift,'theta');
    
    driftspont(cellrange,1:2) = field2array(drift,'spont');
    
    driftwpref(cellrange,1:2) = field2array(drift,'wpref');
    driftwbw(cellrange,1:2) = field2array(drift,'bw') ;
    
    driftF1F0(cellrange,1:2) = field2array(drift,'F1')./field2array(drift,'F0');
    driftF0(cellrange,1:2) = field2array(drift,'F0');
    %       driftorientfreq_all(cellrange,:)=field2array(drift, 'orientfreq_all');
    
    for j = 1:n_units;
        for m = 1:2
        tuning(cellrange(j),m,:)= drift(j,m).thetatuning;
        end
    end
    
    %size(wvform)
    wvform(cellrange,:) = wv';
    
    %get firing rate at all measured orients and SF, put into an array:
    %12 rows(orientations) by 7 columns(SpatialFreqs) for each cell
    
    expt(cellrange)=i;
    
end %%loop over cells







for i = 1:size(driftA1,1)
    for j=1:size(driftA1,2)
        driftA1(i,j);
        driftA2(i,j);
        driftB(i,j);
        drift_theta_w(i,j);
        [OSI(i,j) width(i,j) peak(i,j)] = calculate_tuning(driftA1(i,j),driftA2(i,j),driftB(i,j),drift_theta_w(i,j));
        
        
    end
end

use = peak(:,1)>0 & peak(:,2)>0;

sum(use)

median(peak(use,1))
median(peak(use,2))

mean(peak(use,1))
mean(peak(use,2))

median(peak(use,2)./peak(use,1))
mean(peak(use,2)./peak(use,1))

median(OSI(use,1))
median(OSI(use,2))

col = 'rgbykv';
col = 'kkkkkkkr';
figure
hold on
for i = 1:length(afiles)
    plot(peak(use & expt'==i,1),peak(use&expt'==i,2),'o','Color',col(i));
end

plot([0 40],[0 40])
axis([0 30 0 30])
xlabel('Peak Evoked Firing Rate (Hz) without Stim');
ylabel('Peak Evoked Firing Rate (Hz) with Stim')

%peak bar plot
[bootstat1,bootsam1]=bootstrp(1000,@median,peak(use,1));
peak_se1=std(bootstat1);

[bootstat2,bootsam2]=bootstrp(1000,@median,peak(use,2));
peak_se2=std(bootstat2);

figure;
bar([median(peak(use,1)) median(peak(use,2))]); hold on;
errorbar([median(peak(use,1)) median(peak(use,2))],[peak_se1 peak_se2]);

signrank(OSI(use,1),OSI(use,2))
%OSI plotting
figure;
plot(OSI(use,1),OSI(use,2),'o');

figure;
plot([1 2],[OSI(use,1) OSI(use,2)],'k'); hold on;
plot(ones(size(OSI(use,1))),OSI(use,1),'o','Color','k'); 
plot(2*ones(size(OSI(use,2))),OSI(use,1),'o','Color','k'); 
axis([0 3 -.2 1.2])

[bootstat1,bootsam1]=bootstrp(1000,@median,OSI(use,1));
se1=std(bootstat1);

[bootstat2,bootsam2]=bootstrp(1000,@median,OSI(use,2));
se2=std(bootstat2);

figure;
bar([median(OSI(use,1)) median(OSI(use,2))]); hold on;
errorbar([median(OSI(use,1)) median(OSI(use,2))],[se1 se2]);

signrank(OSI(use,1),OSI(use,2))


