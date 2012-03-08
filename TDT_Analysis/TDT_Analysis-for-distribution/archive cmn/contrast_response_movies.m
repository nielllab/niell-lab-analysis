%function noise_cluster

% Matlab codes for reading from TTank for movie data 
% Calculates RFs from spike triggered average, and spike triggered covariance
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03


%%% read in cluster data, then connect to the tank and read the block
clear all
pack
[fname, pname] = uigetfile('*.mat','cluster data');
oldpname = pname;
load(fullfile(pname,fname));
pname = oldpname;
Block_Name
if nblocks>1
      block = 2;
     Block_Name = char(Block_Name(block)) % to initialize
else
    block=2;
end
%load wn012cpd10Hz075sig5min_0417a   %%% movie frames
%load noisewhite012cpd3hz2std3min
%load noisewhite012cpd5Hz3min_0308

%load wn_movie_5hz_30hz_012cpd_075sig_5min_092106
%load wn012cpd75hz30hz075sig5min_101206
%load wn016cpd30hz75hz075sig3min_011607
% load wn08cpd5hz075sig30hz5min_012307
load wn012alpha1_5hz_contrast10sec5min021407   %%% the usual
%load wn012cpd5hz30hz_alph1_off3_5min_021207
%load alpha_stepcontrast
%load sf_movie
clear amp width x0 baseline
TTX = openTTX(Tank_Name,Block_Name); % to initialize

contrast_modulated = 1;
correct_spectrum = 1;

Event_Name_Snip='Snip'
Event_Name_Wave='PDec'
Sample_Interval=0.04096; % 24414.0625Hz
Sample_Number_Snip=64;
Dec_Factor=32; % 
plot_duration=36; %in second
hist_int = 0.1;
hist_range=[0:hist_int:plot_duration];
axis_range=[0 plot_duration 0 20];
max_events=100000;
n_frames = 5399;

calculate_stc =0;    %%% whether or not to calculate STC (it's slow!)
 xrange =25:45;
yrange = 30:50;   
dx = size(xrange,2);
dy=size(yrange,2);
moviedata = double(moviedata);
m = reshape(moviedata(xrange,yrange,1:n_frames),[dx*dy n_frames])/128;
if calculate_stc
    covmat_prior=0;


    tic
    covmat_prior = cov(double(m)');
    figure
    imagesc(covmat_prior);
    toc
end

m= m(:,1:1000);
mov_var = var(squeeze(moviedata(50,50,:)))
movavg = mean(moviedata,3);

figure
imagesc(movavg-128, [-32/1.414 32/1.414])


tetrode_linear=0;
if tetrode_linear
    ch_map = [14 8 10 4 13 7 9 3 11 1 15 5 12 2 16 6];
else
    ch_map = 1:16;
end
[fname, pname] = uigetfile('*.mat','analysis data');
if fname~=0
  afile = fullfile(pname,fname);
  load(afile); 
  use_afile=1;
else
    use_afile=0;
    %%% select which units to analyze (channel, cluster number)
    cells = [7 1; 10 1; 12 1; 12 2; 14 1; 15 2 ; 15 3; 16 2]
    cells = [1 4; 5 1; 9 1; 9 6; 9 2; 13 2; 13 3; 13 5];
    cells = [1 5; 5 4; 5 6; 9 2; 9 4; 9 6; 13 2; 13 4; 13 6; 13 7 ]
    cells = [1 4; 5 1; 9 1; 9 2; 9 6; 13 1; 13 4; 13 6; 13 7];
    cells = [1 4; 5 4; 5 7; 9 2; 9 7; 13 3; 13 5; 13 7]  %%% 011606_allstim1
    cells = [1 4; 5 4; 5 6; 9 4; 9 6; 13 1; 13 4; ]
    cells = [1 2; 1 4; 1 7; 9 1; 9 2; 9 4; 9 6; 13 1; 13 5; ]  %%% 011606_allstim2
    cells = [1 1;  5 4; 9 4; 9 5; 13 1]
    cells = [1 1; 5 1; 5 2; 5 7; 9 4; 9 5;  9 3; 13 4]; %%%011707_allstim1
    cells = [5 4; 5 6];
    cells = [1 7; 5 7; 9 2; 13 7; 13 1] %% allstim 011806

    cells = [1 1; 1 2; 5 1; 5 2; 5 4; 5 6; 5 7; 9 4; 9 5;  9 3; 9 7; 13 4; 13 3; 13 7]; %%%011707_allstim1
end
%%% set time based on first and last timepoints in clustered data
invoke(TTX,'CreateEpocIndexing')

event_times_all = event_times_all-(block-1)*10^5;
times = event_times_all(event_times_all>0 & event_times_all<10^5);
 time1 = min(times);
 time2= max(times);


MyEpocs = invoke(TTX, 'GetEpocsV', 'fTrg',time1, time2, round(max_time*30));

printfig = input('print figures ? ');


for cell_n = 1:size(cells,1)
%for cell_n=4:6
    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)

    n_spikes= zeros(n_frames,1);
    
   
  channel_times =squeeze(event_times_all(channel_no,:));
  times = channel_times(channel_times>0 & channel_times<10^5);
  time1 = min(times)
  time2=max(times)
    
    sta0 = zeros(size(moviedata(:,:,1)));
    sta1 = sta0; sta2=sta0;
    
    

        clear Epocs trial_no Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
        invoke(TTX,'ResetFilters');
       N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch_map(channel_no), 0, ...
                    time1, time2, 'FILTERED')
        if (N==max_events)
            warning('max number of events aquired');
        end
        Spike_TS = invoke(TTX, 'ParseEvInfoV', 0, N, 6);  %   6  = Time Stamp
        n_reps = max(Spike_TS)/(.03*n_frames)
        
       
        %%% match up spike times for this condition with spike times from clustering
       %%% to identify cluster index for each spike     
%         channel_times =squeeze(event_times_all(channel_no,:));
%         spike_no=zeros(N,1);        
%         n=1;
%         for i = 1:N        
%             while channel_times(n)~=Spike_TS(i)
%                 n = n+1;
%             end
%             spike_no(i) = n;
%         end        
%         clust_idx = idx_all(channel_no,spike_no);

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
       
       
        %%% found out which frame each spike occured during
        
        Spike_Timing = Spike_TS(find(clust_idx==clust_no));
        
        figure
       size(Spike_Timing)
       hist(Spike_Timing,0:60:max(Spike_Timing));
        
        
        frame = zeros(size(Spike_Timing));
        for i=1:size(Spike_Timing,2)   % Future: a better way than this?
            frm = invoke(TTX, 'QryEpocAtV', 'fTrg', Spike_Timing(i), 0);
            frame(i) = frm;
            % The last param: 0 indicating the epoch's value returned;
            % 1 to get the epoch's start time or 2 to get the epoch's stop time.  
            % Pass 3 if you want the epocs filter status.
            if (frm>0) & (frm<=n_frames)
             n_spikes(frm)=n_spikes(frm)+1;
            end
        end
figure
plot(frame, Spike_Timing, '.k', 'MarkerSize',6);

%%% evaluate these before normalizing to rate
N=sum(n_spikes)
sta_N(cell_n)=N;
fano(cell_n) = var(n_spikes)/mean(n_spikes)
    
    if N==0
        break
    end
n_spikes = n_spikes/(.03*n_reps);

clear c
cyc_frames = 10*30;
if contrast_modulated
    for f= 1:n_frames;
        cos(2*pi*f/cyc_frames);
        c(f) = double(0.5- 0.5*cos(2*pi*f/cyc_frames));
     
%         frm = moviedata(:,:,f);
%         c(f) = double(std(frm(:)));
%         figure
%         hist(frm(:)-128,-200:10:200);
%         sum((frm(:)-128).^2)
%         mean(frm(:))
    end
    nf = zeros(300,1);
    contrastdata = zeros(cyc_frames,1);
    for i = 1:n_frames;
        contrastdata(mod(i,cyc_frames)+1)=contrastdata(mod(i,cyc_frames)+1)+double(c(i));
        nf(mod(i,cyc_frames)+1) = nf(mod(i,cyc_frames)+1)+1;
    end
    contrastdata = contrastdata./(nf/cyc_frames);
    fftsignal = (fft(n_spikes));
    fftsignal = fftsignal(1:round(n_frames/2));
    figure
    loglog(abs(fftsignal));
    
    sta_responsiveness(cell_n) = 2*abs(fftsignal(19))/(abs(fftsignal(1)));  %%double to normalize FFT relative to DC
    sta_contrastphase(cell_n) = angle(fftsignal(19));
    phaseangle = angle(fftsignal(19))
    cycledata = zeros(cyc_frames,1);
    for i = 1:n_frames;
        cycledata(mod(i,cyc_frames)+1)=cycledata(mod(i,cyc_frames)+1)+n_spikes(i);
    end
    cycledata = cycledata./(nf/cyc_frames);
 
    %cycledata = conv(cycledata,ones(1,30))/10;
    cycledata = condenseData(cycledata,15);
    contrastdata = condenseData(contrastdata,15);
    
     figure 
    hold on  
    plot(cycledata/max(cycledata),'b');
    plot(double(contrastdata)/max(contrastdata),'g');

    peakresp_wn(cell_n) = max(cycledata);
    p = polyfit(double(contrastdata(1:8))/max(contrastdata),cycledata(1:8)/max(cycledata(1:8)),3);
    figure
     hold on
    plot(double(contrastdata(1:8))/max(contrastdata),cycledata(1:8)/max(cycledata(1:8))); 
    plot(0:0.1:1,polyval(p,double(0:0.1:1)),'go');
   
    contrast_wn(:,cell_n) = contrastdata;
    response_wn(:,cell_n)= cycledata;
    p(4) = p(4)-0.5;
    r = roots(p);
    realroots = find(real(r)==r & r>0);
    halfcontrast_wn(cell_n) = min(r(realroots));
    plot(halfcontrast_wn(cell_n),0.5,'r*');
    halfslope_wn(cell_n) = polyval(polyder(p),halfcontrast_wn(cell_n));
    
    figure
 
    hold on
    
    meancontrast = 0.5*(contrastdata(1:10)+contrastdata(20:-1:11));
    meandata = 0.5*(cycledata(1:10)+cycledata(20:-1:11));
    plot(meancontrast/max(meancontrast),meandata/max(meandata),'g');
    title(sprintf('unit %d %d',channel_no, clust_no))
    xlabel('contrast');
    ylabel('response');
     plot(contrastdata/max(contrastdata),cycledata/max(meandata));
    
   adaptation_index(cell_n) = sum(cycledata(1:10)-cycledata(20:-1:11))/sum(cycledata(1:10)+cycledata(20:-1:11));
     
else

     sta_responsiveness(cell_n) = 0;
    sta_contrastphase(cell_n) = 0;
    halfcontrast_wn(cell_n)=0;
    halfslope_wn(cell_n)=0;
    adaptation_index(cell_n)=0;
end

%%%plot histogram
    titlestr = sprintf('channel %d cluster %d',channel_no, clust_no);
    figure
    plot(n_spikes);
    title(titlestr);
    saveas(gcf,fullfile(pname,sprintf('hist%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');

 
end    %%% cell
peakresp_wn
halfcontrast_wn
adaptation_index

if use_afile
    save(afile,'peakresp_wn','halfslope_wn', 'halfcontrast_wn','contrast_wn','response_wn','adaptation_index','-append');
end
    
invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');