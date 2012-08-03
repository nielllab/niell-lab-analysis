function bars_SU_bw
% Matlab codes for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, then connect to the tank and read the block
clear all
cells =1;
[fname, pname] = uigetfile('*.mat','cluster data');
load(fullfile(pname,fname));
for i =1:length(Block_Name);
    sprintf('%d : %s ',i,Block_Name{i})
end
block = input('which block to analyze ? ');
Block_Name = Block_Name{block}

bw = input('white bars (0) or white/black (1)')+1;

plot_duration=4; %in second
hist_int = 0.1;
hist_range=[0:hist_int:9];
axis_range=[0 plot_duration 0 20];

deg_per_sec=24;

stim_duration = 4;
bar_orients = 0:45:315
blank_stim=1;
if blank_stim
    orient_list = [size(bar_orients,2)  0:size(bar_orients,2)-1]
    %orient_list = [size(bar_orients,2)*2  (0:size(bar_orients,2)-1)+16]
else
    orient_list = [0:size(bar_orients,2)-1]
end

[afname, apname] = uigetfile('*.mat','analysis data');
if afname~=0
    afile = fullfile(apname,afname);
    load(afile);
    use_afile=1;
else
    use_afile=0;
    %%% select which units to analyze (channel, cluster number)
    
end
cells


%printfig = input('print figures? ');
printfig=0;
for cell_n = 1:size(cells,1)
    % for cell_n=9:9
    cell_n
    channel_no = cells(cell_n,1)
    clust_no = cells(cell_n,2)
    channel_times =spikeT{cell_n} - (block-1)*10^5;
    times = channel_times(channel_times>0 & channel_times<10^5);
    
    for rep =1:2
        hist_fig = figure;
        rast_fig = figure;
        for orientation =orient_list;
            
            if rep ==1  &bw==1
                [Spike_Timing index numtrials] = getTrials(stimEpocs{block},times, orientation+1, stim_duration);
            else if rep==1 & bw==2
                    if orientation == size(bar_orients,2)
                        [Spike_Timing index numtrials] = getTrials(stimEpocs{block},times, 2*orientation+1, stim_duration);
                        %%% calculate total spikes
                        R(cell_n,orientation+1) = sum(Spike_Timing<stim_duration)/(stim_duration*numtrials);
                        
                        title_text=['Orientation: ' num2str(orientation*22.5)];
                        
                        if orientation<8
                            %%%raster plot
                            
                            place = [6 3 2 1 4 7 8 9];
                            %place = 1:16;  %%lets you change order of subplots
                            
                            figure(rast_fig);
                            subplot(3,3,place(orientation+1)); hold on; set(gca, 'yDir','reverse');
                            axis([0 plot_duration 0 numtrials+1]);
                            %title (title_text);
                            %
                            %         plot (Spike_Timing, index, '.k', 'MarkerSize',4);
                            
                            plot ([Spike_Timing; Spike_Timing], [index-0.25;index+0.25], 'k', 'MarkerSize',4);
                            set(gca,'XTickLabel',[])
                            set(gca,'YTickLabel',[])
                            
                            %% histograms
                            figure(hist_fig);
                            subplot(3,3,place(orientation+1));
                            bar(hist_range, hist(Spike_Timing, hist_range)/(hist_int*numtrials));  hold on;
                            axis(axis_range);
                            
                            if orientation==0
                                title_text = sprintf('%s',Tank_Name);
                                text(0,25,title_text,'FontSize',8);
                            end
                            if orientation==1
                                title_text = sprintf('ch%d c%d',channel_no,clust_no);
                                text(0,25,title_text,'FontSize',8);
                            end
                        end
                        
                            %%%% curve fitting
                            
                            
                            spont = R(cell_n,size(bar_orients,2)+1);
                            fit_range = 0:0.1:stim_duration;
                            Spike_Timing = Spike_Timing(find((Spike_Timing>min(fit_range))&(Spike_Timing<max(fit_range))));
                            fit_int = fit_range(2)-fit_range(1);
                            obs = hist(Spike_Timing, fit_range)/(fit_int*numtrials);
                            
                            fit_range_interp = 0:.02:Bar_Time(2)+0.5;
                            obs_interp = interp1(fit_range,obs,fit_range_interp);
                            
                            fit_range = fit_range_interp;
                            obs =obs_interp - spont;
                            [min(obs) max(obs) fit_range(find(max(obs))) Bar_Time(2)/5];
                            peak_guess = median(fit_range(find(obs> 0.75*max(obs))));
                            fit_coeff = nlinfit(fit_range,obs,@rf_fit_nobaseline,[ max(obs) peak_guess Bar_Time(2)/10]);
                            if isnan(fit_coeff(1))
                                fit_coeff(1)=0;
                            end
                            baseline(cell_n,orientation+1) = spont;
                            
                            amp(cell_n,orientation+1) = fit_coeff(1);
                            if isnan(fit_coeff(2))
                                amp(cell_n,orientation+1)=0;
                            end
                            if orientation<4
                                x0(cell_n,orientation+1) = fit_coeff(3) ;
                            else
                                x0(cell_n,orientation+1) = stim_duration-fit_coeff(3) ;
                            end
                            width(cell_n,orientation+1) = abs(fit_coeff(3));
                            
                            %%% look for aberrant results, and set all values to zero
                            if abs(fit_coeff(3)>(Bar_Time(2)-Bar_Time(1))) | (fit_coeff(1)<0) | fit_coeff(3)<.045 | sum(isnan(fit_coeff))>0
                                fit_coeff(2)=0;
                                amp(cell_n,orientation+1)=0;
                                x0(cell_n,orientation+1) = 0;
                                width(cell_n,orientation+1) = 0;
                                baseline(cell_n,orientation+1) = mean(obs);
                                fit_coeff(1)=mean(obs);
                            end
                            
                            
                            hold on
                            plot(fit_range, rf_fit_nobaseline(fit_coeff,fit_range)+spont,'g','LineWidth',1.5);
                            
                        end
                        
                    end %orientation
                    saveas(rast_fig,fullfile(pname,sprintf('bar_rast_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
                    saveas(hist_fig,fullfile(pname,sprintf('bar_hist_move%d%s_%d_%d',rep,Block_Name,channel_no,clust_no)),'fig');
                    
                    
                    %     figure
                    %     rates = amp(cell_n,1:16);
                    %     rates(17)=amp(cell_n,1);
                    %     polar([0:22.5:360]*pi/180,rates);
                    %     title_text = sprintf('channel %d cluster %d',channel_no, clust_no);
                    
                    
                    figure(rast_fig)
                    %  title(title_text);
                    if printfig
                        print(hist_fig);
                    end
                    
                    %%%  [bars_theta(cell_n) bars_OSI(cell_n) bars_A1(cell_n) bars_A2(cell_n) bars_w(cell_n) bars_B(cell_n) y] = ...
                    %%%      fit_tuningcurve(R(cell_n,1:size(bar_orients,2))-R(cell_n,size(bar_orients,2)+1), bar_orients*pi/180);
                    
                    [bars_theta(cell_n) bars_OSI(cell_n) bars_A1(cell_n) bars_A2(cell_n) bars_w(cell_n) bars_B(cell_n) bars_null(cell_n) y] = ...
                        fit_tuningcurve( amp(cell_n,1:8), bar_orients*pi/180);
                    
                    theta_ind = 0:pi/4:7*pi/4;
                    i=sqrt(-1);
                    mu = (sum(amp(cell_n,:).*exp(i*theta_ind)))/sum(amp(cell_n,:));
                    if isnan(mu);
                        mu=0;
                    end
                    bars_dsi(cell_n) = abs(mu);
                    
                    %%% find peak orientation, and one on either side (this is a total hack)
                    %%% to average width
                    cond_peak = round(bars_theta(cell_n)/(pi/4))+1
                    if cond_peak>8
                        cond_peak = 1;
                    end
                    cond_down = cond_peak-1;
                    if cond_down <1
                        cond_down=16;
                    end
                    cond_up = cond_peak+1;
                    if cond_up>16
                        cond_up = 1;
                    end
                    conds = [cond_down cond_peak cond_up];
                    
                    width(cell_n,:)
                    
                    rf_width(cell_n) = sum(width(cell_n,conds).*amp(cell_n,conds))./sum(amp(cell_n,conds))
                    
                    %     figure
                    %     plot(amp(cell_n,1:16));
                    %     %%% plot(R(cell_n,1:size(bar_orients,2))-R(cell_n,size(bar_orients,2)+1));
                    %     hold on
                    %     %plot(y,'g');
                    
                    bartuning(cell_n,rep,:)=amp(cell_n,1:8);
                    barspont(cell_n,rep) = spont;
                    barwidth(cell_n,rep,:)=width(cell_n,1:8);
                    
                end   %rep
                
                %  tuningfig=figure
                %  plot(squeeze(bartuning(cell_n,1,:))+barspont(cell_n,1));
                %  hold on
                %  plot(squeeze(bartuning(cell_n,2,:))+barspont(cell_n,2),'g');
                %  plot(ones(length(amp))*barspont(cell_n,1),'b:');
                %  plot(ones(length(amp))*barspont(cell_n,2),'g:');
                %  title_text = sprintf('amplitude channel %d cluster %d',channel_no, clust_no);
                %  title(title_text)
                %  saveas(tuningfig,fullfile(pname,sprintf('bar_tuning_move%s_%d_%d',Block_Name,channel_no,clust_no)),'fig');
                %
                %  figure
                %  plot(squeeze(bartuning(cell_n,1,:).*barwidth(cell_n,1,:)));
                %  hold on
                %  plot(squeeze(bartuning(cell_n,2,:).*barwidth(cell_n,2,:)),'g');
                %   title(title_text)
                %   title_text = sprintf('area channel %d cluster %d',channel_no, clust_no);
            end  %%% cell
            
            figure
            bar(barspont)
            
            bars_theta
            bars_OSI;
            bars_A1
            bars_A2
            bars_w
            bars_B %#ok<NOPRT>
            if blank_stim
                bars_spont=R(:,size(bar_orients,2)+1)
            else
                bars_spont = zeros(size(R(:,size(bar_orients,2))));
            end
            if use_afile
                barsweep_R = R;
                bars_x0 = x0;
                bars_peakwidth = width;
                bars_amp = amp;
                bars_baseline = baseline
                save(afile, 'barsweep_R', 'bars_x0', 'bars_peakwidth', 'bars_amp', 'bars_baseline', ...
                    'bar_orients',  'bars_theta', 'bars_OSI', 'bars_A1', 'bars_A2', 'bars_w', 'bars_B','bars_null', 'bars_spont','rf_width','bars_dsi','-append');
            end
            
            bars_A1./(bars_A1+bars_B)
            bars_OSI
