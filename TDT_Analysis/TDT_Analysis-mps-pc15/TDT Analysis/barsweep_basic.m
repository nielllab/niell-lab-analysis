function barsweep_basic
% Matlab code for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_tetrode.m and select_units.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03
%
%%% includes 2 conditions, which can be decided by any
%%% criterion, but here is generally trackball velocity
%
%%% modified to allow negative peak amplitudes (for suppressive responses)
%
%%% set n_reps=1 for no movement
%
%%% set orientation values in bar_orients
%
%%% if no blank stimulus, set blank_stim=0

stim_duration = 3.015;
bar_orients = 0:22.5:337.5;
thresh_velocity=1;
blank_stim=1;
n_reps=1;

%%% read in cluster data, then connect to the tank and read the block
[fname, pname] = uigetfile('*.mat','cluster data');
load(fullfile(pname,fname));
for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
block = input('which block to analyze ? ');
Block_Name = Block_Name{block}

TTX = openTTX(Tank_Name,Block_Name); % to initialize

if n_reps>1
    [tsamp vsmooth groomstate] = getBlockVelocity_groom(Tank_Name,Block_Name);
end
Event_Name_Snip='Snip'
Sample_Interval=0.04096 % 24414.0625Hz
Sample_Number_Snip=64
Dec_Factor=32; %
plot_duration=stim_duration+0.5; %in second
hist_int = 0.1;
hist_range=[0:hist_int:plot_duration];
axis_range=[0 plot_duration 0 20];
max_events=200000;
deg_per_sec=30;
Bar_Time(1)=0;
Bar_Time(2) =  stim_duration+0.1;

n_orients = length(bar_orients);

ch_map =1:16;

cells=0;
[afname, apname] = uigetfile('*.mat','analysis data');
if afname~=0
    afile_bar = fullfile(apname,afname);
    save(afile_bar, 'afile_bar','apname','-append'); %%% to make sure it's not overwritten
    load(afile_bar);
    use_afile=1;
else
    use_afile=0;
    %%% select which units to analyze (channel, cluster number)
end
cells

%%% if this is a reanalysis, clear old data
clear bars_theta bars_OSI bars_A1 bars_A2 bars_w bars_B rf_width

%%% set time based on first and last timepoints in clustered data
invoke(TTX,'CreateEpocIndexing');


event_times_all = event_times_all-(block-1)*10^5;
times = event_times_all(event_times_all>0 & event_times_all<10^5);
time1 = min(times)
time2= max(times)
MyEpocs = invoke(TTX, 'GetEpocsV', 'xTrg', time1,time2, 1000);

%printfig = input('print figures? ');
printfig=0;
for cell_n = 1:size(cells,1)
    cell_n
    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)

    channel_times =squeeze(event_times_all(channel_no,:));

    times = channel_times(channel_times>0 & channel_times<10^5);
    time1 = min(times)
    time2=max(times)
    if blank_stim
        orient_list = [n_orients  0:n_orients-1]
    else
        orient_list = [0:size(bar_orients,2)-1]
    end
    for rep =1:n_reps
        hist_fig = figure('Name',sprintf('PSTH ch:%d cl:%d',channel_no,clust_no));
        rast_fig = figure( 'Name',sprintf('rasters ch:%d cl:%d',channel_no,clust_no));
        for orientation =orient_list;
            clear Epocs trial_no Epocs_TS N Spike_TS index TS_xTrg Spike_Timing;
            Epocs=find(MyEpocs(1,:)==orientation+1);
            trial_no = length(Epocs);
            Epocs_TS = MyEpocs(2,Epocs);
            invoke(TTX,'ResetFilters');
            ecode= invoke(TTX, 'StringToEvCode', 'xTrg'); % to convert string xTrig to code
            invoke(TTX,'SetFilter', ecode, 69, orientation+1, 0);
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
            numtrials = trial_no;

            if n_reps>1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% movement specific
                clear trial_velocity newtrial
                for i = 1:numtrials;
                    trial_velocity(i) = mean(vsmooth(find(tsamp>Epocs_TS(i) & tsamp<Epocs_TS(i)+stim_duration)));
                end
                trial_velocity
                if rep==1
                    usedtrial = trial_velocity<thresh_velocity;
                else
                    usedtrial = trial_velocity>thresh_velocity;
                end
                trials = find(usedtrial)
                for i = 1:numtrials;
                    if isempty(find(trials==i))
                        newtrial(i)=0;
                    else
                        newtrial(i) = find(trials==i);
                    end
                end
                numtrials = sum(usedtrial);


                %%% only keep spikes that are in the desired cluster
                Spike_Timing = Spike_Timing(find((clust_idx==clust_no)&(usedtrial(index))));
                index=index(find((clust_idx==clust_no)&(usedtrial(index))));
                index = newtrial(index);  %%% get rid of unused trials


                %%%% end movememnt specific

            else  %%no movement
                Spike_Timing = Spike_Timing(find((clust_idx==clust_no)));
                index=index(find((clust_idx==clust_no)));
            end
            %%% calculate total spikes
            R(cell_n,orientation+1) = sum(Spike_Timing<stim_duration)/(stim_duration*numtrials);
            bars_R(cell_n,rep,orientation+1)=R(cell_n,orientation+1);


            if orientation<n_orients
                title_text=['Orientation: ' num2str(bar_orients(orientation+1))];
                %%%raster plot
                place = 1:16;  %%lets you change order of subplots
                %place = [9 10 11 12 13 14 15 16 1 2 3 4 5 6 7 8];
                figure(rast_fig);
                subplot(4,4,place(orientation+1)); hold on; set(gca, 'yDir','reverse');
                axis([0 plot_duration 0 numtrials+1]);

                plot ([Spike_Timing; Spike_Timing], [index-0.25;index+0.25], 'k', 'MarkerSize',4);
                set(gca,'XTickLabel',[])
                set(gca,'YTickLabel',[])
                 if orientation==0
                    title_text = sprintf('%s',Tank_Name);
                    title(title_text,'FontSize',8);
                
            elseif orientation==1
                    title_text = sprintf('ch%d c%d',channel_no,clust_no);
                    title(title_text, 'FontSize',8);
                 else
                title(sprintf('%d',round(bar_orients(orientation+1))),'FontSize',8);
                 end
                %% histograms
                figure(hist_fig);
                subplot(4,4,orientation+1);
                bar(hist_range, hist(Spike_Timing, hist_range)/(hist_int*numtrials));  hold on;
                axis(axis_range);
                if orientation==0
                    title_text = sprintf('%s',Tank_Name);
                   title(title_text,'FontSize',8);
              
                elseif orientation==1
                    title_text = sprintf('ch%d c%d',channel_no,clust_no);
                    title(title_text,'FontSize',8);
                else
                 title(sprintf('%d',round(bar_orients(orientation+1))),'FontSize',8);
                end

                %%%% curve fitting

                spont = R(cell_n,size(bar_orients,2)+1);
                fit_range = 0:0.1:Bar_Time(2);
                Spike_Timing = Spike_Timing(find((Spike_Timing>min(fit_range))&(Spike_Timing<max(fit_range))));
                fit_int = fit_range(2)-fit_range(1);
                obs = hist(Spike_Timing, fit_range)/(fit_int*numtrials);

                fit_range_interp = 0:.02:Bar_Time(2);
                obs_interp = interp1(fit_range,obs,fit_range_interp);

                fit_range = fit_range_interp;
                obs =obs_interp - spont;
                if spont>4 & mean(obs)<0  %%% suppressive responses
                    obs = -obs;
                    inverted=1;
                else inverted=0;
                end

                peak_guess = median(fit_range(find(obs> 0.75*max(obs))));
                [ max(obs) peak_guess Bar_Time(2)/10]
                fit_coeff = nlinfit(fit_range,obs,@rf_fit_nobaseline,[ max(obs) peak_guess Bar_Time(2)/10]);

                if inverted
                    fit_coeff(1)=-fit_coeff(1);
                end
                fit_coeff

                if isnan(fit_coeff(1))
                    fit_coeff(1)=0;
                end
                baseline(cell_n,orientation+1) = spont;

                amp(cell_n,orientation+1) = fit_coeff(1);
                if isnan(fit_coeff(2))
                    amp(cell_n,orientation+1)=0;
                end

                width(cell_n,orientation+1) = abs(fit_coeff(3));

                %%% look for aberrant results, and set all values to zero
                if abs(fit_coeff(3)>(Bar_Time(2)-Bar_Time(1))) | fit_coeff(3)<.045 | sum(isnan(fit_coeff))>0
                    fit_coeff(2)=0;
                    amp(cell_n,orientation+1)=0;
                    x0(cell_n,orientation+1) = 0;
                    width(cell_n,orientation+1) = 0;
                    baseline(cell_n,orientation+1) = mean(obs);
                    fit_coeff(1)=mean(obs);
                end

                hold on
                plot(fit_range, rf_fit_nobaseline_neg(fit_coeff,fit_range)+spont,'g','LineWidth',1.5);
                %plot(fit_range, spont*ones(size(fit_range)),'g','LineWidth',1.5);
            end

        end %orientation
        saveas(rast_fig,fullfile(apname,sprintf('bar_rast_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
        saveas(hist_fig,fullfile(apname,sprintf('bar_hist_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');


        %%% polar plot (untested)
        %     figure
        %     rates = amp(cell_n,1:n_orients);
        %     rates(n_orients+1)=amp(cell_n,1);
        %     polar(bar_orients*pi/180,rates);
        %     title_text = sprintf('channel %d cluster %d',channel_no, clust_no);


        figure(rast_fig)
        if printfig
            print(hist_fig);
        end


        A(1:n_orients/2) = amp(cell_n,(n_orients/2+1):n_orients);
        A((n_orients/2+1):n_orients) = amp(cell_n,1:n_orients/2);

        [bars_theta(cell_n,rep) bars_OSI(cell_n,rep) bars_A1(cell_n,rep) bars_A2(cell_n,rep) bars_w(cell_n,rep) bars_B(cell_n,rep) bars_null(cell_n,rep) bars_curvefit(cell_n,rep,:)] = ...
            fit_tuningcurve_basic( amp(cell_n,1:n_orients), bar_orients*pi/180);

        theta_ind = bar_orients*pi/180
        i=sqrt(-1);
        mu = (sum(amp(cell_n,:).*exp(i*theta_ind)))/sum(amp(cell_n,:));
        if isnan(mu);
            mu=0;
        end
        bars_dsi(cell_n) = abs(mu);

        %%% find peak orientation, and one on either side (this is kind of a hack)
        %%% to average width
        cond_peak = round(bars_theta(cell_n,rep)/((bar_orients(2)-bar_orients(1))*pi/180))+1;
        if cond_peak>n_orients
            cond_peak = 1;
        end
        cond_down = cond_peak-1;
        if cond_down <1
            cond_down=n_orients;
        end
        cond_up = cond_peak+1;
        if cond_up>n_orients
            cond_up = 1;
        end
        conds = [cond_down cond_peak cond_up];

        rf_width(cell_n,rep) = sum(width(cell_n,conds).*amp(cell_n,conds))./sum(amp(cell_n,conds)) * 1.17; %%% 1.17 is ratio between sigma of gaussian, and half-width half-max

        bartuning(cell_n,rep,:)=amp(cell_n,1:n_orients);
        barspont(cell_n,rep) = spont;
        barwidth(cell_n,rep,:)=width(cell_n,1:n_orients);

        [bar_OSI(cell_n,rep) bar_orientwidth(cell_n,rep) bar_peak(cell_n,rep)] = ...
            calculate_tuning(bars_A1(cell_n,rep),bars_A2(cell_n,rep),bars_B(cell_n,rep),bars_w(cell_n,rep));

    end   %rep


    title_text = sprintf('tuning curve ch:%d cl:%d',channel_no, clust_no);
    tuningfig=figure('Name',title_text);
    title(title_text);
    hold on
    plot(bar_orients,squeeze(bartuning(cell_n,1,:))+barspont(cell_n,1),'b');
    plot(bar_orients,squeeze(bars_curvefit(cell_n,1,1:n_orients))+barspont(cell_n,1),'g:');
    plot(bar_orients',ones(size(amp,2))*barspont(cell_n,1),'b:');
 legend('data','fit','spont')
    if n_reps>1
        plot(bar_orients,squeeze(bartuning(cell_n,2,:))+barspont(cell_n,2),'r');
        plot(bar_orients,squeeze(bars_curvefit(cell_n,2,1:n_orients))+barspont(cell_n,2),'g:');
        plot(bar_orients',ones(size(amp,2))*barspont(cell_n,2),'r:');
    end
    xlabel('orientation (deg)');
    ylabel('sp / sec')
    saveas(tuningfig,fullfile(apname,sprintf('bar_tuning_move%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');

end  %%% cell




if use_afile
    barsweep_R = R;
    bars_x0 = x0;
    bars_peakwidth = width;
    bars_amp = amp;
    bars_baseline = baseline;
    save(afile_bar, 'bar_OSI','bar_orientwidth','bar_peak','bars_R','bartuning', 'barwidth', 'barspont', 'bars_theta', 'bars_A1', 'bars_A2', 'bars_w', 'bars_B','bars_null', 'rf_width','bars_dsi','-append');
end

bar_peak
barspont
bar_OSI
bar_orientwidth
rf_width

invoke(TTX, 'CloseTank');
invoke(TTX, 'ReleaseServer');