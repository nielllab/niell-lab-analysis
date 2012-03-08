function shortbar_multiorient
% Matlab codes for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, then connect to the tank and read the block
clear all
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
Sample_Interval=0.04096 % 24414.0625Hz
Sample_Number_Snip=64
Dec_Factor=32; % 
plot_duration=3; %in second
hist_int = 0.2;
hist_range=[0:hist_int:3.2];
axis_range=[0 plot_duration 0 20];
max_events=10000;
deg_per_sec=30;
bar_width = 5;
bar_width_time = bar_width/deg_per_sec;
stim_duration = 3.015;
spacing =8;
nPos=9;
nOrient=4;
blankframe=1;
rf_pix_size = 8;   %%% bin size for calculating rf's (in deg)

orients = [0 45 90 135]
nCond = nPos*nOrient;
if blankframe
    nCond = nCond+1;
end

tetrode_linear=0;

if tetrode_linear
    ch_map = [14 8 10 4 13 7 9 3 11 1 15 5 12 2 16 6];
else
    ch_map = 1:16;
end

%%% select which units to analyze (channel, cluster number)
[afname, apname] = uigetfile('*.mat','analysis data');
if afname~=0
  afile = fullfile(apname,afname)
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
    cells = [ 9 2; 11 6; 12 4; 13 2; 14 5; 15 3]
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



for cell_n = 1:size(cells,1)
    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)
    hist_fig = figure;
    rast_fig = figure;
  channel_times =squeeze(event_times_all(channel_no,:));
  times = channel_times(channel_times>0);
  time1 = min(times);
  time2=max(times);
    for cond =1:nCond;
        
        clear Epocs trial_no Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
        Epocs=find(MyEpocs(1,:)==cond); 
        trial_no = length(Epocs)
        Epocs_TS = MyEpocs(2,Epocs);

        invoke(TTX,'ResetFilters');
        ecode= invoke(TTX, 'StringToEvCode', 'xTrg'); % to convert string xTrig to code
        invoke(TTX,'SetFilter', ecode, 69, cond, 0); 
                                     % 69 means equal to: 'xTrg=orientation', the last parameter not used for '69'
        N = invoke(TTX, 'ReadEventsV', max_events, Event_Name_Snip, ch_map(channel_no), 0, ...
                    max(time1,Epocs_TS(1)),time2, 'FILTERED');
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
        index = zeros(1,N);
        TS_xTrg=index;
        for i = 1:size(Epocs_TS,2)-1;
            epochSpikes = find(Spike_TS>=Epocs_TS(i) & Spike_TS<Epocs_TS(i+1));
            index(epochSpikes)=i;
            TS_xTrg(epochSpikes)=Epocs_TS(i);
        end
        Spike_Timing=Spike_TS-TS_xTrg;

        %%% only keep spikes that are in the desired cluster
        Spike_Timing = Spike_Timing(find(clust_idx==clust_no));
        index=index(find(clust_idx==clust_no));
        size(index)
        title_text=['Orientation: ' num2str(cond*45)];
    
         if cond<=nPos*nOrient
        %%%raster plot
        position = [6 9 8 7 4 1 2 3];
        figure(rast_fig);
        subplot(nOrient,nPos,cond); hold on; set(gca, 'yDir','reverse'); 
        axis([0 plot_duration 0 trial_no+1]); title (title_text);
        plot (Spike_Timing, index, '.k', 'MarkerSize',4);

       %% histograms
        figure(hist_fig);
        subplot(nOrient,nPos,cond);
        bar(hist_range, hist(Spike_Timing, hist_range)/(hist_int*trial_no));  hold on;
        axis(axis_range); title (title_text);
         Bar_Time = episostim_params (cond*90); % Bar appear and disappear time
         plot(Bar_Time(1)*ones(1, length(Screen_Plot_y)), Screen_Plot_y, 'r-'); % Bar appear
         plot(Bar_Time(2)*ones(1, length(Screen_Plot_y)), Screen_Plot_y, 'r-'); % Bar disappear
        
         rf_hist_int = rf_pix_size/deg_per_sec;
         rf_hist_range = 0:rf_hist_int:stim_duration;
         if ~isempty(Spike_Timing(Spike_Timing<stim_duration))
             rf_data(cond,:) = histc(Spike_Timing(Spike_Timing<stim_duration), rf_hist_range)/(rf_hist_int*trial_no);
         else
             rf_data(cond,:) = zeros(size(rf_hist_range));
         end
         axis off
         %%%% curve fitting
    
         fit_range = 0:0.1:Bar_Time(2)+0.5;
         Spike_Timing = Spike_Timing(find((Spike_Timing>min(fit_range))&(Spike_Timing<max(fit_range))));
        fit_int = fit_range(2)-fit_range(1);
        obs = hist(Spike_Timing, fit_range)/(fit_int*trial_no);
        [min(obs) max(obs) fit_range(find(max(obs))) Bar_Time(2)/5];
        peak_guess = median(fit_range(find(obs> 0.75*max(obs))));
        fit_coeff = nlinfit(fit_range,obs,@rf_fit,[min(obs) max(obs) peak_guess Bar_Time(2)/10]);
        baseline(cell_n,cond) = fit_coeff(1);
        amp(cell_n,cond) = fit_coeff(2);
       if cond<4   
        x0(cell_n,cond) = fit_coeff(3) + bar_width_time/2;
       else
        x0(cell_n,cond) = stim_duration-fit_coeff(3) - bar_width_time/2;
       end
        width(cell_n,cond) = abs(fit_coeff(4));

       %%% look for aberrant results, and set all values to zero
       if abs(fit_coeff(4)>(Bar_Time(2)-Bar_Time(1))) | (fit_coeff(2)<0) | ((fit_coeff(4)/fit_coeff(2))>1);
            amp(cell_n,cond)=0;
            x0(cell_n,cond) = 0;
            width(cell_n,cond) = 0;
            baseline(cell_n,cond) = mean(obs);
       end
        if baseline(cell_n,cond)<0
            baseline(cell_n,cond)=0;
        end
        if isnan(baseline(cell_n,cond))
            baseline(cell_n,cond)=0;
        end

        hold on
        plot(fit_range, rf_fit(fit_coeff,fit_range),'g','LineWidth',1.5);
      axis off
      
      
         else    %%% blank frame
             spontmean(cell_n) = sum(Spike_Timing<stim_duration)/(stim_duration*trial_no);
         end
         
         
    end   %orientation

    title_text = sprintf('channel %d cluster %d',channel_no, clust_no);
    figure(rast_fig)
    title(title_text);
    figure(hist_fig);
     title(title_text);
%     xlabel('secs');
%     ylabel('Hz');
    
    if ~blankframe
        spontmean(cell_n) = mean(baseline(cell_n,:));
    end
    
 

   for orient = 1:nOrient
      
  
    rf = rf_data((orient-1)*nPos +1 : orient*nPos,1:size(rf_data,2)-1);
     thresh = spontmean(cell_n) + 0.3*(max(max(rf_data))-spontmean(cell_n))  %%% do this more statistically later, like 1std
     rfthresh = rf - thresh;
    rfthresh = rfthresh.*(rfthresh>0);
    
        n= 0 ; nx = 0; nx2 = 0;
        for i = 1:nPos
         %  ns = sum(rf(i,:))-thresh*size(rf,2);
         ns = sum(rfthresh(i,:));
         
           if ns<0
               ns=0;
           end

           n = n+ns;
            nx = nx+ns*i;
            nx2 = nx2+ns*i*i;
        end;
        centroid(cell_n,orient) = nx/n;
        dispersion(cell_n,orient) = abs(sqrt((nx2/n) - (nx/n)^2))+0.5 ;
        centroid(cell_n,orient)
        dispersion(cell_n,orient)

       rf = imresize(rf,[round((nPos)*spacing) round(size(rf,2)*rf_hist_int*deg_per_sec)],'bilinear');       

        figure
        imagesc(imrotate(rf-spontmean(cell_n),orients(orient),'bilinear','crop'));
      % imagesc(rfthresh) 
        axis equal
        axis([0 90 0 90])
        colorbar
           title(title_text);
        saveas(gcf,fullfile(pname,sprintf('rf%d%s_%d_%d',orient,Block_Name,channel_no,clust_no)),'fig');
 
   end
   
   for orient = 1:nOrient
      rf = rf_data((orient-1)*nPos +1 : orient*nPos,1:size(rf_data,2)-1);
    
       if orient ==1
            rf_obs =rf;
        else
            rf_obs = rf';
       end
%           figure
%        imagesc(rf_obs);
%        title('rf data') 
       
        obs = rf_obs(:);
        [xgrid ygrid]= meshgrid(1:size(rf_obs,2),1:size(rf_obs,1));
    x(:,1) = xgrid(:);
    x(:,2)=ygrid(:);
    [spontmean(cell_n) max(obs)-spontmean(cell_n) ...
        centroid(cell_n,2) dispersion(cell_n,2) centroid(cell_n,1) dispersion(cell_n,1)]
    coeff = nlinfit(x,obs,@rf2fit,[spontmean(cell_n) max(obs)-spontmean(cell_n) ...
        centroid(cell_n,2) dispersion(cell_n,2) centroid(cell_n,1) dispersion(cell_n,1)])
    
   est = rf2fit(coeff,x);
%    figure
%    imagesc(reshape(est,[size(rf_obs,1) size(rf_obs,2)]));
%    title('rf fit');
%     axis equal
       rf_amp(cell_n,orient) = coeff(2);
   spont(cell_n,orient) = coeff(1);
   x0(cell_n,orient) = coeff(3);
   wx0(cell_n,orient) = abs(coeff(4));
   y0(cell_n,orient) = coeff(5);
   wy0(cell_n,orient) = abs(coeff(6));
       
   end
 
    
    
    
    
  
    saveas(rast_fig,fullfile(pname,sprintf('rast%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
    saveas(hist_fig,fullfile(pname,sprintf('hist%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
    
end  %%% cell

   dpos=(size(rf_obs,1)-size(rf_obs,2))/2 %%% how much to shift centers to account fo different size each direction
    y0(:,1) = y0(:,1)+dpos;
    x0(:,2) = x0(:,2)+dpos;
    ampsum = rf_amp(:,1) + rf_amp(:,2);
    xmean = (x0(:,1).*rf_amp(:,1) + x0(:,2).*rf_amp(:,2))./ampsum;
    ymean = (y0(:,1).*rf_amp(:,1) + y0(:,2).*rf_amp(:,2))./ampsum;
    wxmean = (wx0(:,1).*rf_amp(:,1) + wx0(:,2).*rf_amp(:,2))./ampsum;
    wymean = (wy0(:,1).*rf_amp(:,1) + wy0(:,2).*rf_amp(:,2))./ampsum;
    
    for c = 1:size(cells,1)
        sprintf('channel %d cluster %d',cells(c,1),cells(c,2))
        sprintf('x0 %0.2f   y0 %0.2f   wx %0.2f   wy %0.2f',xmean(c),ymean(c),wxmean(c),wymean(c))
    end

baseline;
amp;
x0=x0*deg_per_sec;
width=width*deg_per_sec;
save(fullfile(pname,sprintf('analysis%s_%s',Tank_Name,Block_Name)));

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');