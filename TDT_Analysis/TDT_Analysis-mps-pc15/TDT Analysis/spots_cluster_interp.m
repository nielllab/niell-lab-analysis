function spots_cluster_interp

% Matlab codes for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, then connect to the tank and read the block
clear all
[fname, pname] = uigetfile('*.mat','cluster data');
load(fullfile(pname,fname));
if exist('nblocks')
    
      block =6;
     Block_Name = char(Block_Name(block)) % to initialize
     
else
    block=1;
end
Tank_Name
Block_Name
TTX = openTTX(Tank_Name,Block_Name); % to initialize

Event_Name_Snip='Snip'
Sample_Interval=0.04096 % 24414.0625Hz
Sample_Number_Snip=64
Dec_Factor=32; % 

stim_duration = 0.2;
wait_duration = 0.2;
latency = 0.025;
stim_off = 0.2;


plot_duration=stim_duration + wait_duration + latency; %in second
hist_int = 0.1;
hist_range=[0:hist_int:plot_duration];
axis_range=[0 plot_duration 0 20];
max_events=10000;
deg_per_sec=30;
bar_width = 5;
bar_width_time = bar_width/deg_per_sec;

spacing =6;
nx = 10;
ny=10;
blankframe=1;
onset_only =1;
offset_only=0;

nCond = nx*ny;
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
invoke(TTX,'CreateEpocIndexing');
event_times_all = event_times_all-(block-1)*10^5;
times = event_times_all(event_times_all>0 & event_times_all<10^5);
 times = times(times>0);
 time1 = min(times);
 time2= max(times);
MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', time1,time2, 1000);


Screen_Plot_y = [axis_range(3):0.1:axis_range(4)]; % to mark the bar appearing and disappearing time



for cell_n = 1:size(cells,1)
    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)
  
    rast_fig = figure;
  channel_times =squeeze(event_times_all(channel_no,:));
  times = channel_times(channel_times>0);
  time1 = min(times);
  time2=max(times);
    for cond =1:nCond;
        
        clear Epocs trial_no Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
        Epocs=find(MyEpocs(1,:)==cond); 
        trial_no = length(Epocs);
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
        spike_no=zeros(N,1);        
        n=1;
        for i = 1:N        
            while channel_times(n)~=Spike_TS(i)
                n = n+1;
            end
            spike_no(i) = n;
        end        
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
        
        title_text=['Orientation: ' num2str(cond*45)];
    
         if cond<=nx*ny
        %%%raster plot
        position = [6 9 8 7 4 1 2 3];
        figure(rast_fig);
        subplot(nx,ny,cond); hold on; set(gca, 'yDir','reverse'); 
        axis([0 plot_duration 0 trial_no+1]);
        %axis off
        plot (Spike_Timing, index, '.k', 'MarkerSize',4);
    set(gca,'XTickLabel',[])
    set(gca,'YTickLabel',[])

       %% histograms
%         figure(hist_fig);
%         subplot(nx,ny,cond);
%         bar(hist_range, hist(Spike_Timing, hist_range)/(hist_int*trial_no));  hold on;
%         axis(axis_range); axis off
%          Bar_Time = episostim_params (cond*90); % Bar appear and disappear time
%          plot(Bar_Time(1)*ones(1, length(Screen_Plot_y)), Screen_Plot_y, 'r-'); % Bar appear
%          plot(Bar_Time(2)*ones(1, length(Screen_Plot_y)), Screen_Plot_y, 'r-'); % Bar disappear
%         
         %rf_data(cond,:) = hist(Spike_Timing(Spike_Timing<stim_duration), hist_range)/(hist_int*trial_no);
         if onset_only
             rf_data(ceil(cond/ny), mod(cond-1,ny)+1) = sum(Spike_Timing>latency & Spike_Timing<(stim_duration))/(trial_no*(stim_duration-latency));
         elseif offset_only
              rf_data(ceil(cond/ny), mod(cond-1,ny)+1) = sum(Spike_Timing>stim_off & Spike_Timing<(stim_off+ stim_duration))/(trial_no*(stim_duration));
         else
             rf_data(ceil(cond/ny), mod(cond-1,ny)+1) = sum(Spike_Timing>latency & Spike_Timing<(stim_duration+wait_duration))/(trial_no*(stim_duration-latency+wait_duration));
         end
            
%          %%%% curve fitting
%     
%          fit_range = 0:0.1:Bar_Time(2)+0.5;
%          Spike_Timing = Spike_Timing(find((Spike_Timing>min(fit_range))&(Spike_Timing<max(fit_range))));
%         fit_int = fit_range(2)-fit_range(1);
%         obs = hist(Spike_Timing, fit_range)/(fit_int*trial_no);
%         [min(obs) max(obs) fit_range(find(max(obs))) Bar_Time(2)/5];
%         peak_guess = median(fit_range(find(obs> 0.75*max(obs))));
%         fit_coeff = nlinfit(fit_range,obs,@rf_fit,[min(obs) max(obs) peak_guess Bar_Time(2)/10]);
%         baseline(cell_n,cond) = fit_coeff(1);
%         amp(cell_n,cond) = fit_coeff(2);
%        if cond<4   
%         x0(cell_n,cond) = fit_coeff(3) + bar_width_time/2;
%        else
%         x0(cell_n,cond) = stim_duration-fit_coeff(3) - bar_width_time/2;
%        end
%         width(cell_n,cond) = abs(fit_coeff(4));
% 
%        %%% look for aberrant results, and set all values to zero
%        if abs(fit_coeff(4)>(Bar_Time(2)-Bar_Time(1))) | (fit_coeff(2)<0) | ((fit_coeff(4)/fit_coeff(2))>1);
%             amp(cell_n,cond)=0;
%             x0(cell_n,cond) = 0;
%             width(cell_n,cond) = 0;
%             baseline(cell_n,cond) = mean(obs);
%        end
%         if baseline(cell_n,cond)<0
%             baseline(cell_n,cond)=0;
%         end
%         if isnan(baseline(cell_n,cond))
%             baseline(cell_n,cond)=0;
%         end
% 
%         hold on
%         plot(fit_range, rf_fit(fit_coeff,fit_range),'g','LineWidth',1.5);

      
      
         else    %%% blank frame
%              figure
%              hold on
%              set(gca, 'yDir','reverse'); 
%         axis([0 plot_duration 0 trial_no+1]);
%         axis off
%         plot (Spike_Timing, index, '.k', 'MarkerSize',4);
             spontmean(cell_n) = sum(Spike_Timing<stim_duration)/(stim_duration*trial_no)
         end
         
         
    end   %orientation

    title_text = sprintf('channel %d cluster %d',channel_no, clust_no);
    figure(rast_fig)
    title(title_text);

%     xlabel('secs');
%     ylabel('Hz');
    
    if ~blankframe
        spontmean(cell_n) = 0;
    end
      
  
    rf = rf_data';
    
        interp_mag=10;
    rf = imresize(rf,interp_mag,'bicubic');
    
     thresh = 1.5* spontmean(cell_n) + 0.2*(max(max(rf_data))-1.5*spontmean(cell_n))  %%% do this more statistically later, like 1std
     rfthresh = rf - thresh;
    rfthresh = rfthresh.*(rfthresh>0);
%     figure
%     imagesc(rf_data')
%     title('rf_data')

    

    
    for ax=1:2
        n= 0 ; nsx = 0; nsx2 = 0;
        if ax==1
            nr = nx*interp_mag;
        else
            nr = ny*interp_mag;
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
   theta = -5 * (pi/180);
   
    x(:,1) = cos(theta)*xgrid(:) -sin(theta)*ygrid(:);
    x(:,2)= sin(theta)*xgrid(:) + cos(theta)*ygrid(:);
    
   c1 = cos(theta)*centroid(cell_n,1)-sin(theta)*centroid(cell_n,2);
   c2 = sin(theta)*centroid(cell_n,1) + cos(theta)*centroid(cell_n,2);
    
   options = statset('Robust','off');
   
   prior_coeff = [spontmean(cell_n) max(obs)-spontmean(cell_n) ...
        c1 dispersion(cell_n,1) c2 dispersion(cell_n,2)];
   coeff = nlinfit(x,obs,@rf2fit,prior_coeff,options);
   
    
   est = rf2fit(coeff,x);

   figure
   imagesc(reshape(est,[size(rf,1) size(rf,2)]));
   title('rf fit');
   
   amp(cell_n) = coeff(2);
   spont(cell_n) = coeff(1);
   x0(cell_n) = coeff(3)/interp_mag;
   wx0(cell_n) = coeff(4)/interp_mag;
   y0(cell_n) = coeff(5)/interp_mag;
   wy0(cell_n) = coeff(6)/interp_mag;
   
   
       rf = imresize(rf,[round(ny*spacing) round(nx*spacing)],'bilinear');       

        figure
        imagesc(rf-spontmean(cell_n));
      % imagesc(rfthresh) 
        axis equal
        axis([0 nx*spacing 0 nx*spacing])
        colorbar
           title(title_text);
        saveas(gcf,fullfile(pname,sprintf('rf%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');

    
    saveas(rast_fig,fullfile(pname,sprintf('rast%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');

    
end  %%% cell
wx0./wy0

x0(abs(x0)>100)=0;
y0(abs(y0)>100)=0;
wx0(abs(wx0)>100)=0;
wy0(abs(wy0)>100)=0;
x0
y0
wx0
wy0
clear event_times_all etimes_old idx_all tm used score c_score channel_times

%save(fullfile(pname,sprintf('spotanalysis%s_%s',Tank_Name,Block_Name)));

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');