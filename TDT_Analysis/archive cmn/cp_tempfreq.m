function cp_tempfreq
% Matlab codes for reading from TTank for drifting or counterphase gratings
% plots histgrams and rasters and performs fourier analysis
%function gratfreq_cluster
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03

%%% read in cluster data, then connect to the tank and read the block

clear all
pack
[fname, pname] = uigetfile('*.mat','cluster data');
oldpname = pname;
load(fullfile(pname,fname));   %%% need to copy pname, or it can get written over in load
pname = oldpname;
if exist('nblocks')
      block = 	4;
     Block_Name = char(Block_Name(block)) % to initialize
     
else
    block=1;  
end

    TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_Snip='Snip'
Event_Name_Wave='PDec'
Sample_Interval=0.04096 % 24414.0625Hz
Sample_Number_Snip=64
Dec_Factor=32; % 
plot_duration=2.5; %in second

hist_int = 0.05;
hist_range=[0:hist_int:plot_duration];
axis_range=[0 plot_duration 0 20];
max_events=50000;

stim_duration =2;
wait_duration = .25;  %duration after each stimulus (for calculating spontaneous rate)
blank_interval = 0.1; %% length of time after stimulus not to use in calculating spontaneous rate
fft_int = .05;  %%% size of bins to use for fourier analyis (in sec)
tempfreq = 2;    %%% temporal frequency of stimulus
blank_stim = 1; %%% is there an extra stimulus to measure spontaneous
full_field = 0;  %%% is there a full-field flicker?
tetrode_linear=0;

if tetrode_linear
    ch_map = [14 8 10 4 13 7 9 3 11 1 15 5 12 2 16 6];
else
    ch_map = 1:16;
end

[fname, pname] = uigetfile('*.mat','analysis data');
if fname~=0
  afile = fullfile(pname,fname)
  load(afile); 
  use_afile=1;
else
    use_afile=0;

    %%% select which units to analyze (channel, cluster number)
    cells_bad = [1 5; 13 5;]
    cells_incomplete = [13 3]
    %cells_mod = [1 2; 1 5; 1 4;  9 1; 9 7;  9 4; 9 2; 13 2; 13 6; 13 4];
    cells_good = [1 6; 1 7; 5 1; 5 3; 9 2 ; 9 4; 9 5; 13 2; 13 3; 13 4; 13 7];

    cells = [1 4; 5 1; 9 1; 9 6; 9 2; 13 2; 13 3; 13 5]
    cells = [1 4; 5 4; 5 7; 9 2; 9 7; 13 3; 13 5; 13 7]  %%% 011606_allstim1
    cells = [1 2; 1 4; 1 7; 9 1; 9 2; 9 4; 9 6; 13 1; 13 5; ]  %%% 011606_allstim2
    cells = [1 3; 1 6 ; 1 2; 5 7; 5 4; 5 6; 9 1; 9 7;  13 3; 13 7]

    cells = [5 4 ; 5 6];
    cells = [5 4; 13 2; 13 4; ]
    cells = [1 7; 5 7; 9 2; 13 7; 13 1] %% allstim 011806
    cells = [1 1; 1 2; 5 1; 5 2; 5 4; 5 6; 5 7; 9 4; 9 5;  9 3; 9 7; 13 4; 13 3; 13 7]; %%%011707_allstim1
    cells = [5 6]
end


%%% set time based on first and last timepoints in clustered data
invoke(TTX,'CreateEpocIndexing');
event_times_all = event_times_all-(block-1)*10^5;
times = event_times_all(event_times_all>0 & event_times_all<10^5);

 times = times(times>0);
 time1 = min(times)
 time2= max(times)
MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', time1,time2, 1000);
Screen_Plot_y = [axis_range(3):0.1:axis_range(4)]; % to mark the bar appearing and disappearing time
 
%%% set number of conditions and display setup (generally rows = orientation, columns = frequency)
n_rows=8;
n_col=4;
%orients = [0 22.5 45 67.5 90 112.5 135 157.5];
orients = [ 0 22.5 45 67.5 90 112.5 135 157.5];
tempfreqs = [1 2 4 8];

n_cond=n_rows*n_col;
if blank_stim
    n_cond=n_cond+1;
end
if full_field
    n_cond=n_cond+1;
end
if ~exist('cptf_F0')
    cptf_F0 = zeros(1,size(cells,1));
    cptf_F1 = cptf_F0;
    cptf_F2 = cptf_F0;
    cptf_tftuning = zeros(size(cells,1),4);
end

printfig = input('print figures ? ');


for cell_n = 1:size(cells,1)

   channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)
    hist_fig = figure;
    rast_fig = figure;
    fft_fig = figure;
    spont_full= figure;
  channel_times =squeeze(event_times_all(channel_no,:));

  times = channel_times(channel_times>0);
  time1 = min(times);
  time2=max(times);
  
    for cond =0:n_cond-1
       tempfreq = tempfreqs(mod(cond,n_col)+1);
       clear Epocs trial_no Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
        Epocs=find(MyEpocs(1,:)==cond+1); 
        trial_no = length(Epocs);
        Epocs_TS = MyEpocs(2,Epocs);

        invoke(TTX,'ResetFilters');
        ecode= invoke(TTX, 'StringToEvCode', 'xTrg'); % to convert string xTrig to code
        invoke(TTX,'SetFilter', ecode, 69, cond+1, 0); 
                                     % 69 means equal to: 'xTrg=orientation', the last parameter not used for '69'
        N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch_map(channel_no), 0, ...
                    max(time1,Epocs_TS(1)),time2,  'FILTERED');
        if (N==max_events)
            warning('max number of events aquired');
        end
        Spike_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp

  
      %%% match up spike times for this condition with spike times from clustering
       %%% to identify cluster index for each spike
       channel_times =squeeze(event_times_all(channel_no,:));
        ntimes = size(event_times_all,2);
        spike_no=zeros(N,1);        
        n=1;
        for i = 1:N        
            while channel_times(n)~=Spike_TS(i) & n<ntimes
                n = n+1;
            end
            if channel_times(n)==Spike_TS(i)
                spike_no(i) = n;
            else    %% if not caught, start over at beginning for next search
                n=1;
                spike_no(i)=ntimes+1;
            end
                
        end        
       idx_all(channel_no,ntimes+1)=0;  %%%catchall for spikes not identified
       clust_idx = idx_all(channel_no,spike_no);


        % Now for all events, we find out when the xTrig was:
%         index = zeros(1,N);
%         TS_xTrg=index;
%         for i = 1:size(Epocs_TS,2)-1;
%             epochSpikes = find(Spike_TS>=Epocs_TS(i) & Spike_TS<Epocs_TS(i+1));
%             index(epochSpikes)=i;
%             TS_xTrg(epochSpikes)=Epocs_TS(i);
%         end

        index = zeros(1,N);
        TS_xTrg=index;
        for i = 1:size(Epocs_TS,2);
            if i<size(Epocs_TS,2)
                epochSpikes = find(Spike_TS>=Epocs_TS(i) & Spike_TS<Epocs_TS(i+1));
            else
              epochSpikes = find(Spike_TS>=Epocs_TS(i));
            end
            index(epochSpikes)=i;
            TS_xTrg(epochSpikes)=Epocs_TS(i);
        end
        
        
        
        Spike_Timing=Spike_TS-TS_xTrg;

        Spike_Timing = Spike_Timing(find(clust_idx==clust_no));
        index=index(find(clust_idx==clust_no));
        size(index);
        
        
        title_text=['Orientation: ' num2str(cond*45)];

      %raster plot 
      if cond<n_rows*n_col;
        figure(rast_fig);     
        subplot(n_rows,n_col,cond+1); hold on; set(gca, 'yDir','reverse'); 

        plot (Spike_Timing, index, '.k', 'MarkerSize',4);
        axis([0 plot_duration 0 trial_no+1]);
                 set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
      end
       %% histograms
        rate_hist = hist(Spike_Timing, hist_range)/(hist_int*trial_no);
        if cond<n_rows*n_col
            figure(hist_fig);
            subplot(n_rows,n_col,cond+1);
            bar(hist_range, rate_hist);  hold on;
            axis(axis_range); 
             set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
        end
           if cond==0
                title_text = sprintf('%s',Tank_Name);
                text(0,20,title_text,'FontSize',8);
            end
            if cond==1
              title_text = sprintf('ch%d c%d',channel_no,clust_no);
                text(0,20,title_text,'FontSize',8);  
            end
        
        if cond>=n_rows*n_col & full_field%%% blank frame or full-field, but only plot if full-field is on
            figure(spont_full)
            subplot(2,1,cond-n_rows*n_col +1);
            bar(hist_range, rate_hist);  hold on;
            axis(axis_range); 
             set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
        end
            

        %%fourier analysis
        fft_range = 0:fft_int:stim_duration-fft_int;
        fft_range = fft_range + (fft_int/2);  %%% to center bins
        
        nyq = 0.5/fft_int; %%% nyquist interval
        freq_int = nyq / (0.5*stim_duration/fft_int);

 
        fft_data=abs(fft(hist(Spike_Timing(Spike_Timing<stim_duration),fft_range)/(fft_int*trial_no)));
        fft_data = fft_data/size(fft_data,2);
        fft_data(1) = fft_data(1)/2;
        if cond<n_rows*n_col;
            figure(fft_fig);
            subplot(n_rows,n_col,cond+1);
            plot(fft_data(1:round((size(fft_data,2)+1)/2)));
            axis([0 16 0 8])
        end
        
        R(cell_n,cond+1) = sum(Spike_Timing<=stim_duration)/(stim_duration*trial_no);
        spont(cell_n,cond+1) = sum((Spike_Timing>(stim_duration+blank_interval))&(Spike_Timing<(stim_duration+wait_duration)))/(trial_no*(wait_duration-blank_interval));
        
        F0(cell_n,cond+1) = fft_data(1);
        F1(cell_n,cond+1) = fft_data(1+tempfreq/freq_int);  %%% double to count both pos & neg frequency
        F2(cell_n,cond+1) = fft_data(1+2*tempfreq/freq_int);
        
        
    end

    title_text = sprintf('channel %d cluster %d',channel_no, clust_no);

   if blank_stim
       spont_avg = R(cell_n,n_rows*n_col+1)
   else
       spont_avg = mean(spont(cell_n,:))  %%% if no blank frame, average over inter-stimulus interval of all conditions
        
   end
    plotcolor = 'bgrcmykbgr';
   
    tuning_fig = figure;
   
    
    %for f= 1:n_col;
     for f = 1:n_col
         plot(R(cell_n,f:n_col:f+n_col*(n_rows-1))-spont_avg,plotcolor(f));
        hold on;
    end
 title(title_text);
 legend('.02cpd','.04cpd','.08cpd','.16cpd')
 %set(gca,'XTickLabel',['0' '45' '90' '135' '180' '225' '270' '315']);
 
 
 %%% calculate tuning parameters
 
 
 
 
%%% reshape, and use average across SF to determine preferred angle
    n_orient = n_rows;
    orientfreq(:,1:n_col) = reshape(R(cell_n,1:n_col*n_rows),n_col,n_rows)'-spont_avg;
    n_w = n_col;
     orient_tuning_mn = mean(orientfreq,2)';
%      figure
%      plot(orient_tuning_mn);
     theta_ind = orients*pi/180;
     [theta_pref OSI A1 A2 w B yfit] = fit_tuningcurve(orient_tuning_mn,theta_ind);
    theta_pref = mod(theta_pref,pi) %%% because only goes from 0 to 180

%      hold on
%      plot(yfit,'g')
    figure(hist_fig);
     ginput(1);
     use = input('use for analysis? 0/1');
    if use && printfig
        print(hist_fig);
    end
        %%% find preferred SF at optimal theta
         [w theta] = meshgrid(1:n_w,theta_ind);
        if theta_pref >max(max(theta))
            theta_pref = max(max(theta))
        end
         
        theta_pref_col = ones(size(1:n_w));
        theta_pref_col(:) = theta_pref;
        w_tuning_curve = interp2(w,theta,orientfreq,1:n_w,theta_pref_col);

        figure
        plot(w_tuning_curve);
        if use
            cptf_tftuning(cell_n,:) = w_tuning_curve;
            cptf_tftuning
        end
        [amp w_pref] = max(w_tuning_curve);


        %%% fit orientation tuning curve at optimal SF
        w_pref_col = ones(size(1:n_orient)).*w_pref;
        theta_tuning = interp2(w,theta,orientfreq,w_pref_col,theta_ind);
        figure
        plot(theta_tuning);
        [cptf_theta(cell_n) cptf_osi(cell_n) cptf_A1(cell_n) cptf_A2(cell_n) ...
            cptf_thetawidth(cell_n) cptf_B(cell_n) cptf_null(cell_n) yfit]= fit_tuningcurve(theta_tuning,theta_ind);
        hold on
        plot(yfit,'g');
         cptf_wpref(cell_n)=tempfreqs(w_pref)
         cptf_theta
         cptf_theta = mod(cptf_theta,pi);

w_pref
cptf_theta(cell_n)
if cptf_theta(cell_n)>max(max(theta))  %%% to avoid extrapolating past end theta
    cptf_theta(cell_n) = max(max(theta))
end

     if use
        %%% interpolate F0, F1 values
        F0all(:,1:n_col) = reshape(F0(cell_n,1:n_col*n_rows),n_col,n_rows)'-spont_avg;
         cptf_F0(cell_n) = interp2(w,theta,F0all,w_pref, cptf_theta(cell_n));
        F1all(:,1:n_col) = reshape(F1(cell_n,1:n_col*n_rows),n_col,n_rows)'
        cptf_F1(cell_n) = interp2(w,theta,F1all,w_pref, cptf_theta(cell_n))
        F2all(:,1:n_col) = reshape(F2(cell_n,1:n_col*n_rows),n_col,n_rows)'
        cptf_F2(cell_n) = interp2(w,theta,F2all,w_pref, cptf_theta(cell_n))
     end
 

    figure(rast_fig)
    title(title_text);
    figure(hist_fig);
    title(title_text);
    figure(fft_fig);
    title(title_text);
    xlabel('secs');
        saveas(tuning_fig,fullfile(pname,sprintf('grattuning%s_%d_%d',Block_Name,channel_no,clust_no)),'fig')
        saveas(rast_fig,fullfile(pname,sprintf('gratrast%s_%d_%d',Block_Name,channel_no,clust_no)),'fig')
    saveas(hist_fig,fullfile(pname,sprintf('grathist%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
    saveas(fft_fig,fullfile(pname,sprintf('gratfft%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
    cell_n
end
clear event_times_all etimes_old idx_all tm used score c_score channel_times
%save(fullfile(pname,sprintf('gratanalysis%s_%s',Tank_Name,Block_Name)));
if use_afile
    cp_tf_tempfreq_all = orientfreq;
    cptf_spont = R(:,n_cond);  

    
    save(afile, 'cp_tf_tempfreq_all', 'cptf_theta','cptf_A1','cptf_A2','cptf_thetawidth','cptf_wpref','cptf_null','cptf_tftuning', ...
        'cptf_spont', 'cptf_F0','cptf_F1','cptf_F2', '-append');
end
cptf_spont
cptf_F0
cptf_F1
cptf_F2

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');