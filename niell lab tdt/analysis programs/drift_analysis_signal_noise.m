%function drift_orientfreq_laser
% Matlab codes for reading from TTank for drifting or counterphase gratings
% plots histgrams and rasters and performs fourier analysis
%function gratfreq_cluster
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03

%%% read in cluster data, then connect to the tank and read the block



clear all

SU = 1;
if SU
    [fname, pname] = uigetfile('*.mat','cluster data');
    load(fullfile(pname,fname));
    for i =1:length(Block_Name);
        sprintf('%d : %s ',i,Block_Name{i})
    end
    block = input('which block to analyze ? ');
    Block_Name = Block_Name{block}
    [afname, apname] = uigetfile('*.mat','analysis data');
    noisepname = apname;
    afile = fullfile(apname,afname);
    load(afile);
    use_afile=1;
    cells
else
    pname = uigetdir('C:\data\TDT tanks','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))
end


chans = 1:4:max(cells,1);

%laser = input('movement (0) or laser (1) : ');
laser = 0
if laser
    flags = struct('laserOn',1,'visStim',1)
    tdtData= getTDTdata(Tank_Name, Block_Name, chans, flags);
    tsamp  = tdtData.laserT;
    vsmooth = tdtData.laserTTL;
else
    flags = struct('mouseOn',1,'visStim',1)
    tdtData= getTDTdata(Tank_Name, Block_Name, chans, flags);
    tsamp = tdtData.mouseT;
    vsmooth = tdtData.mouseV;
end


% [fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);  %%% get ps filename
% %psfilename = 'c:/test.ps';   %%% default location
% if exist(psfilename,'file')==2;delete(psfilename);end %%% check for previous file

thresh_velocity = 0.7;
figure
plot(tsamp,vsmooth);


plot_duration=2.5; %in second

hist_int = 0.05;
hist_range=[0:hist_int:plot_duration];
axis_range=[0 plot_duration 0 25];
max_events=50000;

stim_duration =1.5;
wait_duration = .5;  %duration after each stimulus (for calculating spontaneous rate)
blank_interval = 0.1; %% length of time after stimulus not to use in calculating spontaneous rate
fft_int = .05;  %%% size of bins to use for fourier analyis (in sec)
tempfreq = 2;    %%% temporal frequency of stimulus
blank_stim = 1; %%% is there an extra stimulus to measure spontaneous
full_field = 1;  %%% is there a full-field flicker?



% set number of conditions and display setup (generally rows = orientation, columns = frequency)
n_rows=12;
n_col=6;
%orients = [0 45 90 135 180 225 270 315];
orients = 0:30:330;
spatfreqs = [.01 .02 .04 .08 .16 .32];
%



n_cond=n_rows*n_col;
if blank_stim
    n_cond=n_cond+1;
end
if full_field
    n_cond=n_cond+1;
end

%printfig = input('print ? ');
printfig=0;

if SU
    cell_range = 1:size(cells,1)
else
    cell_range=1:4:nchan;
end
for cell_n = cell_range;
    % for cell_n=9:9
    cell_n
    if SU
        channel_no = cells(cell_n,1)
        clust_no = cells(cell_n,2)
        channel_times =spikeT{cell_n} - (block-1)*10^5;
        times = channel_times(channel_times>0 & channel_times<10^5);
        hist_fig = figure('Name',sprintf('unit %d %d',channel_no,clust_no))
    else
        hist_fig = figure('Name',sprintf('channel %d',cell_n))
        channel_no = cell_n;
        clust_no = [];
    end
    
    for rep =1:2
        
        
        for cond =0:n_cond-1
                  if SU
                [Spike_Timing index numtrials Epocs_TS] = getTrialsSU(stimEpocs{block},times, cond+1, stim_duration);
            
            else
               [Spike_Timing index numtrials Epocs_TS] = getTrialsSU(tdtData.stimEpocs,tdtData.MUspikeT{cell_n}, cond+1, stim_duration);
                  end
           
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
            
            
            
            cond
            numtrials
            index
            R(cell_n,cond+1) = sum(Spike_Timing<=stim_duration)/(stim_duration*numtrials); %mean firing rate per condition
           
%             for i = 1:numtrials;
%                 R_trial(i) = sum(Spike_Timing(find(usedtrial(index==i))))/1.5;
%             end
%             vari=var(R_trial);
%             R_var(cell_n,cond+1)=vari;
             
            spont(cell_n,cond+1) = sum((Spike_Timing>(stim_duration+blank_interval))&(Spike_Timing<(stim_duration+wait_duration)))/(numtrials*(wait_duration-blank_interval));
            
            
            
            clear spikeR
            for i = 1:numtrials;
                spikeR(i) = sum(Spike_Timing<=stim_duration & index==i)/(stim_duration);
                
            end
            
            %%mean, var and err for each condition 
            R_mean(cell_n,cond+1)=mean(spikeR);
            R_var(cell_n,cond+1) = var(spikeR);
            R_err(cell_n,cond+1) = std(spikeR)/sqrt(numtrials -1);
            R_std(cell_n,cond+1)=std(spikeR);
            
            cond
%             mean(spikeR)
%             var(spikeR)
%             std(spikeR)
            spikeR
        end  %% cond
        
        title_text = sprintf('channel %d cluster %d',channel_no, clust_no);
        
        if blank_stim
            spont_avg = R(cell_n,n_rows*n_col+1);
            spont_var=var(R(cell_n,n_rows*n_col+1));
            both_drift_spont(cell_n,rep)=R(cell_n,n_rows*n_col+1);
            % spont_avg = 0.5*( R(cell_n,n_rows*n_col+1) +  mean(spont(cell_n,:)))
            
        else
            spont_avg = mean(spont(cell_n,:))  %%% if no blank frame, average over inter-stimulus interval of all conditions
        end
        
        
         [drift_peakR(cell_n) maxcond] = max(R_mean(cell_n,1:n_cond-2));
          
         
         drift_sig_noise(cell_n,rep).peakR= R_mean(cell_n,maxcond); 
         drift_sig_noise(cell_n,rep).peakvar= R_var(cell_n,maxcond);
         drift_sig_noise(cell_n,rep).peakerr= R_err(cell_n,maxcond);
         drift_sig_noise(cell_n,rep).peakstd=R_std(cell_n,maxcond);
         
         drift_sig_noise(cell_n,rep).signoise= R_mean(cell_n,maxcond)/R_std(cell_n,maxcond);
         drift_sig_noise(cell_n,rep).signoise_SE= R_mean(cell_n,maxcond)/R_err(cell_n,maxcond);
         
       
        
%         plotcolor = 'bgrcmykbgr';
%         
%         tuning_fig = figure;
%         
%         
%         %for f= 1:n_col;
%         for f = 1:n_col
%             plot(R(cell_n,f:n_col:f+n_col*(n_rows-1))-spont_avg,plotcolor(f));
%             hold on;
%         end
%         title(title_text);
%         legend('.01 cpd','.02cpd','.04cpd','.08cpd','.16cpd','.32cpd')
       
       
        cell_n
%     
%   
   
    end %%% rep
    
%    

 close all
end


if use_afile
    
    
    
    save(afile, 'drift_sig_noise','-append');
end

