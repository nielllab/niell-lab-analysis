%function generic_analysis
% Matlab codes for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, then connect to the tank and read the block
close all
clear all
SU = input('multiunit (0) or single-unit (1) : ');

cells =1;
if SU
    [fname, pname] = uigetfile('*.mat','cluster data');
    load(fullfile(pname,fname));
    for i =1:length(Block_Name);
        sprintf('%d : %s ',i,Block_Name{i})
    end
    block = input('which block to analyze ? ');
    Block_Name = Block_Name{block}
    [afname, apname] = uigetfile('*.mat','analysis data');
    
    afile = fullfile(apname,afname);
    load(afile);
    use_afile=1;
    cells
    nchan=max(cells,1)
else
    pname = uigetdir('D:\','block data')
    delims = strfind(pname,'\');
    selected_path = pname(1 :delims(length(delims))-1)
    Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
    Block_Name = pname(delims(length(delims))+1 :length(pname))
    nchan = input('number of channels : ');
    flags = struct('visStim',1,'MUspike',1);
    data = getTDTdata(Tank_Name,Block_Name,1:4:nchan,flags);
    
end

   flags = struct('mouseOn',1);
    data = getTDTdata(Tank_Name,Block_Name,1:4:nchan,flags);
    tsamp = data.mouseT;
    vsmooth = data.mouseV;
    thresh_velocity = 1.0;
figure
plot(tsamp,vsmooth);
    
stim_duration = input('duration : ');
done=0;
nEpocs = max(stimEpocs{block}(1,:))

while ~done
nrows = input('rows : ');
ncols = input('cols : ');
panels = input('panels : ');
if nrows*ncols*panels<=nEpocs
    done=1;
else
    display('this is more than # of stimuli - try again!')
end
end

plot_duration=stim_duration; %in second

hist_int = plot_duration/20;
hist_range=[0:hist_int:plot_duration];
if SU
    axis_range=[0 plot_duration 0 25];
else
    axis_range=[0 plot_duration 0 100];
end
if SU
    cell_range = 1:size(cells,1)
else
    cell_range=1:4:nchan;
end
for cell_n = cell_range;
  for runrep=1:2
      % for cell_n=9:9
    cell_n
    if SU
        channel_no = cells(cell_n,1)
        clust_no = cells(cell_n,2)
        channel_times =spikeT{cell_n} - (block-1)*10^5;
        times = channel_times(channel_times>0 & channel_times<10^5);
        times = times;
     %   hist_fig = figure('Name',sprintf('unit %d %d',channel_no,clust_no))
    else
       % hist_fig = figure('Name',sprintf('channel %d',cell_n))
    end
    for rep =1:panels
        
        if SU
            rast_fig = figure('Name',sprintf('unit %d %d rep %d',channel_no,clust_no,runrep));
           % timefig = figure('Name',sprintf('unit %d %d rep %d',channel_no,clust_no,rep));
        else
            rast_fig = figure('Name',sprintf('unit %d rep %d',cell_n,runrep));
             timefig = figure('Name',sprintf('unit %d rep %d',cell_n,runrep));
        end
        for c =1:nrows*ncols;
            orientation = (c-1)*panels+rep;
            if SU
                [Spike_Timing index numtrials Epocs_TS] = getTrialsSU(stimEpocs{block},times, orientation, stim_duration);
            else
                [Spike_Timing index numtrials Epocs_TS] = getTrialsSU(data.stimEpocs,data.MUspikeT{cell_n}, orientation, stim_duration);
            end
            
            if runrep==2
               % keyboard
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% movement specific
            clear trial_velocity newtrial
            numtrials
            for i = 1:numtrials;
                trial_velocity(i) = mean(vsmooth(find(tsamp>Epocs_TS(i) & tsamp<Epocs_TS(i)+stim_duration)));
            end
            trial_velocity
            if runrep==1
                usedtrial = trial_velocity<thresh_velocity;
            else
                usedtrial = trial_velocity>thresh_velocity;
            end
        
            trials = find(usedtrial)                   
            if isempty(trials)
                [m trials] = min(trial_velocity);
                trials
                usedtrial(trials)=1;
            end
        
            for i = 1:numtrials;
                if isempty(find(trials==i))
                    newtrial(i)=0;
                else
                    newtrial(i) = find(trials==i);
                end
            end
       
            numtrials = sum(usedtrial);
         %%% only keep spikes that are in the desired cluster

            Spike_Timing = Spike_Timing(find(usedtrial(index)));
            index=index(usedtrial(index));
            index = newtrial(index);  %%% get rid of unused trials
            
            %%%% end movememnt specific
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% calculate total spikes
           aud_duration=0.2;
           R(cell_n,orientation,runrep) = sum(Spike_Timing<stim_duration)/(stim_duration*numtrials);
            if orientation<=8
                Revoked(cell_n,orientation,runrep)= R(cell_n,orientation,runrep);
            else
                Revoked(cell_n,orientation,runrep) = sum(Spike_Timing<aud_duration)/(aud_duration*numtrials);
            end     
            title_text=['Orientation: ' num2str(orientation*22.5)];
            
            
            figure(rast_fig);
            subplot(nrows,ncols,c)
            hold on; set(gca, 'yDir','reverse');
            axis([0 plot_duration 0 numtrials+1]);
            %title (title_text);
            %
            %         plot (Spike_Timing, index, '.k', 'MarkerSize',4);
            
            plot ([Spike_Timing; Spike_Timing], [index-0.25;index+0.25], 'k', 'MarkerSize',4);
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
            
            %% histograms
%             figure(hist_fig);
%             subplot(nrows,ncols,c);
%             hold on
%             color_list = 'brg';
%         color = color_list(rep);
%             plot(hist_range, hist(Spike_Timing, hist_range)/(hist_int*numtrials),color);
%             hold on;
%             axis(axis_range);
%             set(gca,'XTickLabel',[])
%             set(gca,'YTickLabel',[])
            %             if orientation==1
            %                 title_text = sprintf('ch%d c%d rep%d',channel_no,clust_no, rep);
            %             end
            
%             figure(timefig);
%             for i=1:numtrials;
%                 r(i)=sum(index==i);
%                 subplot(nrows,ncols,c);
%                 plot(r);
%                 axis([0 length(r) 0 1+max(r)*1.1])
%             end
%             set(gca,'XTickLabel',[])
%             set(gca,'YTickLabel',[])

        end %orientation
 
        if SU
            saveas(rast_fig,fullfile(pname,sprintf('generic_rast_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
           % saveas(hist_fig,fullfile(pname,sprintf('generic_hist_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
        end
        % if use_afile
        %     barsweep_R = R;
        %     bars_x0 = x0;
        %     bars_peakwidth = width;
        %     bars_amp = amp;
        %     bars_baseline = baseline
        %     save(afile, 'barsweep_R', 'bars_x0', 'bars_peakwidth', 'bars_amp', 'bars_baseline', ...
        %         'bar_orients',  'bars_theta', 'bars_OSI', 'bars_A1', 'bars_A2', 'bars_w', 'bars_B','bars_null', 'bars_spont','rf_width','bars_dsi','-append');
        % end
        
    end
  end
         figure('Name',sprintf('unit %d %d',channel_no,clust_no));
        plot(squeeze(Revoked(cell_n,:,:)));
        hold on
        plot(1:size(Revoked,2),Revoked(cell_n,9,1)*ones(size(Revoked,2),1),'b-');
        plot(1:size(Revoked,2),Revoked(cell_n,9,2)*ones(size(Revoked,2),1),'g-');
        plot(9,Revoked(cell_n,9,1),'bo');
        plot(9,Revoked(cell_n,9,2),'go');
  legend('stationary','moving')
end

