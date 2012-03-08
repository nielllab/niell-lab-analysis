function shortbar_multiorient_interp
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
      block = 	3;
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

printfig = input('print figures ? ');


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




%         % Now for all events, we find out when the xTrig was:
%         index = zeros(1,N);
%         TS_xTrg=index;
%         for i = 1:size(Epocs_TS,2)-1;
%             epochSpikes = find(Spike_TS>=Epocs_TS(i) & Spike_TS<Epocs_TS(i+1));
%             index(epochSpikes)=i;
%             TS_xTrg(epochSpikes)=Epocs_TS(i);
%         end
%        
        
        
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
           if cond==1
                title_text = sprintf('%s',Tank_Name);
                text(0,20,title_text,'FontSize',8);
            end
            if cond==2
              title_text = sprintf('ch%d c%d',channel_no,clust_no);
                text(0,20,title_text,'FontSize',8);  
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
    if printfig
        print(hist_fig);
    end
%     xlabel('secs');
%     ylabel('Hz');
    
    if ~blankframe
        spontmean(cell_n) = mean(baseline(cell_n,:));
    end
    
 
   for orient = 1:nOrient

    rf = rf_data((orient-1)*nPos +1 : orient*nPos,1:size(rf_data,2)-1);    
    nx = size(rf,2);
    ny=size(rf,1);
    interp_mag=10;
    rf = imresize(rf,interp_mag,'bicubic');
    rf = imrotate(rf-spontmean(cell_n),orients(orient),'bilinear','crop');
    
        thresh = 1.5* spontmean(cell_n) + 0.2*(max(max(rf_data))-1.5*spontmean(cell_n))  %%% do this more statistically later, like 1std
     rfthresh = rf - thresh;
    rfthresh = rfthresh.*(rfthresh>0); 

    %%% initial guess at centroid, dispersion
    for ax=1:2
        n= 0 ; nsx = 0; nsx2 = 0;
      
        if ax==1
                nr = size(rfthresh,2);
            else
                nr = size(rfthresh,1);
            end
        for i = 1:nr
             %  ns = sum(rf(i,:))-thresh*size(rf,2);
             if ax ==1
                 ns = sum(rfthresh(:,i));
             else
                 ns = sum(rfthresh(i,:));
             end
               if ns<0
                   ns=0;
               end
               n = n+ns;
                nsx = nsx+ns*i;
                nsx2 = nsx2+ns*i*i;
        end;
        centroid(cell_n,ax) = nsx/n
        dispersion(cell_n,ax) = max(sqrt((nsx2/n) - (nsx/n)^2) , interp_mag* 0.5)
            
    end 
    
%     f = fspecial('gaussian',3,.5);
%     f = f/f(2,2);
%     rf = imfilter(rf,f);
    

    [xgrid ygrid]= meshgrid(1:nx*interp_mag,1:ny*interp_mag);
    obs = rf(:);
   theta = 0;
   
    x(:,1) = cos(theta)*xgrid(:) -sin(theta)*ygrid(:);
    x(:,2)= sin(theta)*xgrid(:) + cos(theta)*ygrid(:);
    
   c1 = cos(theta)*centroid(cell_n,1)-sin(theta)*centroid(cell_n,2);
   c2 = sin(theta)*centroid(cell_n,1) + cos(theta)*centroid(cell_n,2);
    
   options = statset('Robust','off');
   
   prior_coeff = [spontmean(cell_n) max(obs)-spontmean(cell_n) ...
        c1 dispersion(cell_n,1) c2 dispersion(cell_n,2)];
   coeff = nlinfit(x,obs,@rf2fit,prior_coeff,options);

      est = rf2fit(coeff,x);
  est = imresize(reshape(est,[size(rf,1) size(rf,2)]),[round(ny*spacing) round(nx*spacing)],'bilinear'); 
 
   
       rf = imresize(rf,[round(ny*spacing) round(nx*spacing)],'bilinear');       

        figure
        imagesc(rf-spontmean(cell_n));
        axis equal
        axis([0 nx*spacing 0 nx*spacing])
        colorbar
           title(title_text);
        saveas(gcf,fullfile(pname,sprintf('rf%s_%d_%d_%d',Block_Name,channel_no,clust_no,orient)),'fig');

  figure
   imagesc(est);
   title('rf fit');
   axis equal
   axis([0 nx*spacing 0 nx*spacing])
   colorbar

   sb_amp(cell_n,orient) = coeff(2);
  % sb_spont(cell_n,orient) = coeff(1);
   sb_x0(cell_n,orient) = coeff(3)*spacing/interp_mag;
   sb_wx0(cell_n,orient) = coeff(4)*spacing/interp_mag;
   sb_y0(cell_n,orient) = coeff(5)*spacing/interp_mag;
   sb_wy0(cell_n,orient) = coeff(6)*spacing/interp_mag;
   sb_allrf(cell_n,:,:) = rf-spontmean(cell_n);
   sb_spont(cell_n) = spontmean(cell_n);
   

   end %%% orientation
   done=0;
   n=0;
   clear used
   while done ==0;
        u =input('orientation to use : ');
        if u==0;
            done=1;
        else
            n=n+1;
            used(n)=u
        end
   end
    if n>0
        sb_x0mean(cell_n) = sum(sb_x0(cell_n,used).*sb_amp(cell_n,used))./sum(sb_amp(cell_n,used));
        sb_y0mean(cell_n)  = sum(sb_y0(cell_n,used).*sb_amp(cell_n,used))./sum(sb_amp(cell_n,used));
        sb_wx0mean(cell_n)  = sum(sb_wx0(cell_n,used).*sb_amp(cell_n,used))./sum(sb_amp(cell_n,used));
        sb_wy0mean(cell_n)  = sum(sb_wy0(cell_n,used).*sb_amp(cell_n,used))./sum(sb_amp(cell_n,used));
        sb_ampmean(cell_n)  = mean(sb_amp(cell_n,used));
    else
        sb_x0mean(cell_n)  = NaN;
        sb_y0mean(cell_n)  = NaN;
        sb_wx0mean(cell_n)  = NaN;
        sb_wy0mean(cell_n)  =NaN;  
        sb_ampmean(cell_n)  = 0;
    end
    saveas(rast_fig,fullfile(pname,sprintf('rast%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
    saveas(hist_fig,fullfile(pname,sprintf('hist%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
    
end  %%% cell

    if use_afile
        save(afile,'sb_amp','sb_spont','sb_x0','sb_y0','sb_wx0','sb_wy0', 'sb_x0mean','sb_y0mean','sb_wx0mean','sb_wy0mean','sb_ampmean','-append')
    end
    
sb_x0mean
sb_y0mean
sb_wx0mean
sb_wy0mean
sb_ampmean

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');