function bars_analysis(clustfile,afile,pdfFile,Block_Name,blocknum)
% Matlab codes for reading from TTank for sweeping bars in 8 orientations
% plots histgrams and rasters and fits data to a gaussian peak
% Uses clustering information from cluster_linear.m or cluster_tetrode.m
% cmn 06-06, based on code by Jianhua Cang 06-27-03



%%% read in cluster data, then connect to the tank and read the block

if ~exist('Block_Name','var');
    SU = menu('recording type','multi-unit','single unit')-1;
    useArgin=0;
else
    SU=1;
    useArgin=1;
end

cells =1;
if SU
    if ~useArgin
        [fname, pname] = uigetfile('*.mat','cluster data');
        clustfile=fullfile(pname,fname);
    end
    load(clustfile);
    if ~useArgin
        blocknum = listdlg('ListString',Block_Name,'SelectionMode','single');
        
        
        [afname, apname] = uigetfile('*.mat','analysis data');
        noisepname = apname;
        afile = fullfile(apname,afname);
    end
    load(afile);
    [pname fname] = fileparts(afile);
    Block_Name = Block_Name{blocknum}
    use_afile=1;
    cells
else
    pname = uigetdir('C:\data\','block data')
    delims = strfind(pname,'\');
    selected_path = pname(1 :delims(length(delims))-1)
    Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
    Block_Name = pname(delims(length(delims))+1 :length(pname))
    nchan = input('number of channels : ');
    flags = struct('visStim',1,'MUspike',1);
    data = getTDTdata(Tank_Name,Block_Name,1:4:nchan,flags);
    
end

%stim_duration = input('duration [1] : ');
% nrows = input('# orients [8] : ');
% ncols = input('# sfs [6]  : ');
% tf = input('temp freqs [2 8] : ');
% latency = input('latency (0.05) : ');

prompt = {'duration','# orients','# sfs','temp freqs'};
num_lines = 1;
def = {'4','2','4','[1 2]'};
if ~useArgin
    answer = inputdlg(prompt,'grating parameters',num_lines,def);
else
    answer=def;
end
stim_duration = str2num(answer{1})
nrows = str2num(answer{2})
ncols = str2num(answer{3})
bw =  str2num(answer{4})


if useArgin
    psfilename = [pdfFile(1:end-4) 'bars.ps'];
else
    [fname pname] =uiputfile('*.ps'); psfilename=fullfile(pname,fname);
end
if exist(psfilename,'file')==2;delete(psfilename);end


panels= length(bw);

plot_duration=stim_duration; %in second

hist_int = plot_duration/20;
hist_range=[0:hist_int:plot_duration];
axis_range=[0 plot_duration 0 30];

if SU
    cell_range = 1:size(cells,1)
else
    cell_range=1:4:nchan;
end

orient_list = linspace(0,360,nrows*ncols+1);
orient_list = orient_list(1:end-1);

for cell_n = cell_range;
    %for cell_n=21:26
    cell_n
    if SU
        channel_no = cells(cell_n,1)
        clust_no = cells(cell_n,2)
        channel_times =spikeT{cell_n} - (blocknum-1)*10^5;
        times = channel_times(channel_times>0 & channel_times<10^5);
        hist_fig = figure('Name',sprintf('unit %d %d',channel_no,clust_no))
        epocs = stimEpocs{blocknum}
    else
        hist_fig = figure('Name',sprintf('channel %d',cell_n))
        epocs=data.stimEpocs;
    end
    
    %%% spont and flicker
    spontfig=figure
    emax = max(epocs(1,:));
    
    extra_range = 1:(emax-panels*nrows*ncols)
    
  
        cond = panels*nrows*ncols+1;
        
        if SU
            [Spike_Timing index numtrials] = getTrialsSU(epocs,times, cond, stim_duration);
        else
            [Spike_Timing index numtrials] = getTrialsSU(epocs,data.MUspikeT{cell_n}, cond, stim_duration);
        end
        
        spikes=Spike_Timing(:);
        spikes=spikes(spikes<stim_duration);
        
        %%% rasters
        figure(spontfig);
       
        hold on; set(gca, 'yDir','reverse');
        axis([0 plot_duration 0 numtrials+1]);
        plot ([Spike_Timing; Spike_Timing], [index-0.25;index+0.25], 'k', 'MarkerSize',4);
        set(gca,'XTickLabel',[]);     set(gca,'YTickLabel',[])
        
     for rep = 1:panels
    bars(cell_n,rep).spont = length(Spike_Timing)/(stim_duration*numtrials);
     end
    
    for rep =1:panels
        
        if SU
            rast_fig = figure('Name',sprintf('unit %d %d rep %d',channel_no,clust_no,rep));
        else
            rast_fig = figure('Name',sprintf('unit %d rep %d',cell_n));
        end
        
        for c =1:nrows*ncols;
            cond = (c-1)*panels+rep;
            
            orientation=c;
            if SU
                [Spike_Timing index numtrials] = getTrialsSU(stimEpocs{blocknum},times, cond, stim_duration);
            else
                [Spike_Timing index numtrials] = getTrialsSU(data.stimEpocs,data.MUspikeT{cell_n}, cond, stim_duration);
            end
            %%% calculate total spikes
            spikes=Spike_Timing(:);
            spikes=spikes(spikes<stim_duration);
            R = length(spikes)/(stim_duration*numtrials);
            
            bars(cell_n,rep).R(c)=R;
            
            figure(rast_fig);
            subplot(nrows,ncols,c)
            hold on; set(gca, 'yDir','reverse');
            axis([0 plot_duration 0 numtrials+1]);
            plot ([Spike_Timing; Spike_Timing], [index-0.25;index+0.25], 'k', 'MarkerSize',4);
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
            
            %% histograms
            figure(hist_fig);
            subplot(nrows,ncols,c);
            hold on
            if rep ==1
                color = 'b';
            else
                color = 'r';
            end
            histrate = hist(Spike_Timing, hist_range)/(hist_int*numtrials);
            plot(hist_range, histrate,color);
            hold on;
            plot([0 max(hist_range)], [bars(cell_n,1).spont bars(cell_n,1).spont],'g')
            axis(axis_range);
            set(gca,'XTickLabel',[])
            set(gca,'YTickLabel',[])
            
            bars(cell_n,rep).Rmax(c) = max(histrate);
            
            fit_range = 0:0.2:stim_duration;
            Spike_Timing = Spike_Timing(find((Spike_Timing>min(fit_range))&(Spike_Timing<max(fit_range))));
            fit_int = fit_range(2)-fit_range(1);
            obs = hist(Spike_Timing, fit_range)/(fit_int*numtrials);

            fit_range_interp = 0:.02:stim_duration;
            obs_interp = interp1(fit_range,obs,fit_range_interp,'pchip');

            fit_range = fit_range_interp;
            obs =obs_interp - bars(cell_n,1).spont;
            [min(obs) max(obs) fit_range(find(max(obs))) stim_duration/5];
            peak_guess = median(fit_range(find(obs> 0.75*max(obs))));
            fit_coeff = nlinfit(fit_range,obs,@rf_fit_nobaseline,[ max(obs) peak_guess stim_duration/10]);
            
            fit_coeff(3)=abs(fit_coeff(3));
            if isnan(fit_coeff(1))
                fit_coeff(1)=0;
            end
            baseline(cell_n,orientation) = bars(cell_n,1).spont;

            amp(cell_n,orientation) = fit_coeff(1);
            if isnan(fit_coeff(2))
                amp(cell_n,orientation)=0;
            end
            if orientation<4
                x0(cell_n,orientation) = fit_coeff(2) ;
            else
                x0(cell_n,orientation) = stim_duration-fit_coeff(2) ;
            end
            width(cell_n,orientation) = fit_coeff(3);

            %% look for aberrant results, and set all values to zero
            if abs(fit_coeff(3)>stim_duration) | (fit_coeff(1)<0) | fit_coeff(3)<.045 | sum(isnan(fit_coeff))>0
                fit_coeff(2)=0;
                amp(cell_n,orientation)=0;
                x0(cell_n,orientation) = 0;
                width(cell_n,orientation) = 0;
                baseline(cell_n,orientation) = mean(obs);
                fit_coeff(1)=mean(obs);
            end


            hold on
            plot(fit_range, rf_fit_nobaseline(fit_coeff,fit_range)+bars(cell_n,1).spont,'g','LineWidth',1);

            
            
        end %orientation
        
        bars(cell_n,rep).amp = squeeze(amp(cell_n,:));
        bars(cell_n,rep).x0 = squeeze(x0(cell_n,:));
        bars(cell_n,rep).width = squeeze(amp(cell_n,:));
       
        title(sprintf('ch %d cl %d',channel_no,clust_no));
        
    end  %%% panel
    
    figure(hist_fig) 
    set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    
    figure
    plot(orient_list, bars(cell_n,1).amp);
    hold on
    %plot(orient_list,bars(cell_n,1).Rmax-bars(cell_n,1).spont,'--');

    plot(orient_list,bars(cell_n,2).amp,'r');
   % plot(orient_list,bars(cell_n,2).Rmax-bars(cell_n,1).spont,'r--');

     set(gca,'XTick',orient_list)
    xlabel('orientation -deg');
    
       set(gcf, 'PaperPositionMode', 'auto');
    print('-dpsc',psfilename,'-append');
    
end

save(afile,'bars','-append');

ps2pdf('psfile', psfilename, 'pdffile', [psfilename(1:(end-2)) 'pdf']);
delete(psfilename);
