clear all

pname = uigetdir('C:\data\TDT tanks','block data')
delims = strfind(pname,'\');
selected_path = pname(1 :delims(length(delims))-1)
Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
Block_Name = pname(delims(length(delims))+1 :length(pname))

% Tank_Name='06222011_mlr_stim_awake'
% Block_Name='wn3_72hz'

nchan = input('# chans : ');
uselaser = input('laser on 0/1 : ');
chans = 1:4:nchan;

SU = input('multiunit (0) or single-unit (1) : ');
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
end

if uselaser
    flags = struct('lfpTseries',1,'lfpSpectra',1','mouseOn',1,'laserOn',1,'MUspike',1,'visStim',1)
else
    flags = struct('mouseOn',1,'MUspike',1,'visStim',1)
end

tdtData= getTDTdata(Tank_Name, Block_Name, chans, flags);
%tdtData= getTDTdata(Tank_Name, Block_Name, nChan, 0, 0, 0, 1,1, 0, 0, 0);

if uselaser
    smoothwindow_secs = 0.25;
    laserT = tdtData.laserT;
    dt = laserT(2)-laserT(1);
    smoothwindow = ceil(smoothwindow_secs/dt)
    lasersmooth = zeros(size(tdtData.laserTTL));
    %%% replace this with convolution (for speed)!
    for i = smoothwindow+1:length(tdtData.laserTTL);
        lasersmooth(i)= mean(tdtData.laserTTL(i-smoothwindow:i));
    end
    lasersmooth = lasersmooth/5;
    figure
    plot(lasersmooth)
    
    tic
    tsamp = 0;
    npulse = 0;
    onT=0;
    for i = 2:length(lasersmooth);
        if lasersmooth(i-1)==0 & lasersmooth(i)>0
            npulse = npulse+1;
            onT(npulse) = laserT(i);
        end
    end
    toc
    
    timeRange = -10:0.1:25;
    vdata=zeros(length(timeRange),npulse-1);
    for i = 2:npulse-2
        vdata(:,i) = interp1(tdtData.mouseT,tdtData.mouseV,timeRange+onT(i));
    end
    figure
    plot(timeRange,vdata)
    hold on
    plot([0 17],[40 40],'g','LineWidth',12)
    
    figure
    plot(timeRange,mean(vdata,2))
    hold on
    plot([0 17],[20 20],'g','LineWidth',12)
    
    
    
    for ch = chans;
        df = median(diff(tdtData.spectF{ch}));
        specdata = tdtData.spectData{ch};
        normalizer = 1:size(specdata,2);
        normalizer = repmat(normalizer,size(specdata,1),1);
        specdata = specdata.*normalizer;
        
        tdata = tdtData.spectT{ch};
        spect_avg = zeros(length(timeRange),size(specdata,2));
        for i = 2:npulse-1;
            for f=1:size(specdata,2);
                spect_avg(:,f) = spect_avg(:,f) + interp1(tdata,squeeze(specdata(:,f)),timeRange+onT(i))';
            end
        end
        figure
        imagesc(spect_avg',1.5*[0 prctile(spect_avg(:),95)])
        axis xy
        df = median(diff(tdtData.spectF{ch}));
        dt = median(diff(timeRange));
        set(gca,'YTick',(10:10:80)/df);
        set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
        set(gca,'XTick',(5:5:40)/dt);
        set(gca,'XTickLabel',{'-5','0','5','10','15','20','25','30'})
    end
    
end

bins = 0:0.25:max(tdtData.mouseT);
R= zeros(length(bins),length(chans));
for ch = chans;
    R(:,ch)= hist(tdtData.MUspikeT{ch},bins);
end
figure
plot(R)

f= tdtData.frameEpocs(1,:);
t= tdtData.frameEpocs(2,:);
if uselaser
    laserInterp = interp1(laserT,lasersmooth,t);
end

mouseInterp = interp1(tdtData.mouseT,tdtData.mouseV,t);
figure
hist(diff(t))

dt = (diff(t));
dt_avg = median(dt);
tmin = -10;
tmax = 22;
dt = 1;
tbins = tmin:dt:tmax-1;
cycframes = 300;
framescale = 30;
cycbins = cycframes/framescale;

if SU
    cell_range = 1:size(cells,1)
else
    cell_range=1:4:nchan
end
for cell_n = cell_range
    ch= cell_n;
    if SU
        channel_no = cells(cell_n,1)
        clust_no = cells(cell_n,2)
        channel_times =spikeT{cell_n} - (block-1)*10^5;
        times = channel_times(channel_times>0 & channel_times<10^5);
        frame_duration = median(diff(frameEpocs{block}(2,:)))
    else
        clust_no = [];
        channel_no = cell_n;
        frame_duration = median(diff(tdtData.frameEpocs(2,:)));
        times=tdtData.MUspikeT{cell_n};
    end
    Rtimecycle = zeros(length(tbins),cycbins);
    timecyc_trials = zeros(length(tbins),cycbins);
    for rep = 1:2
        R=zeros(1,max(f));
        
        
        Rcycle = zeros(1,cycframes);
        ntrial = zeros(1,max(f));
        cyctrials = zeros(1,cycframes);
        
        movRcycle = zeros(1,cycframes);
        movntrial = zeros(1,max(f));
        movcyctrials = zeros(1,cycframes);
        
        for i = 1:length(tdtData.frameEpocs)-1
            fspikes = sum(times>t(i) & times<t(i)+dt_avg);
            R(f(i)) = R(f(i))+fspikes;
            if rep ==1 & uselaser
                
                modcycle = mod(floor(f(i)/framescale),cycbins)+1;
                postlaserT = t(i) - max(onT(onT<t(i)));
                prelaserT = t(i)- min(onT(onT>t(i)));
                if postlaserT<tmax
                    tbin = max(find(tbins<postlaserT));
                    Rtimecycle(tbin,modcycle) = Rtimecycle(tbin,modcycle)+ fspikes;
                    timecyc_trials(tbin,modcycle)= timecyc_trials(tbin,modcycle)+1;
                end
                if prelaserT>tmin
                    tbin = max(find(tbins<prelaserT));
                    Rtimecycle(tbin,modcycle) = Rtimecycle(tbin,modcycle)+ fspikes;
                    timecyc_trials(tbin,modcycle)= timecyc_trials(tbin,modcycle)+1;
                end
            end
            if uselaser
                if rep ==1 & laserInterp(i)==0
                    ntrial(f(i))=ntrial(f(i))+1;
                    Rcycle(mod(f(i),cycframes)+1) = Rcycle(mod(f(i),cycframes)+1) + sum(times>t(i) & times<t(i)+dt_avg);
                    cyctrials(mod(f(i),cycframes)+1) = cyctrials(mod(f(i),cycframes)+1) +1;
                elseif rep==2 & laserInterp(i)>0.01;
                    ntrial(f(i))=ntrial(f(i))+1;
                    Rcycle(mod(f(i),cycframes)+1) = Rcycle(mod(f(i),cycframes)+1) + sum(times>t(i) & times<t(i)+dt_avg);
                    cyctrials(mod(f(i),cycframes)+1) = cyctrials(mod(f(i),cycframes)+1) +1;
                end
            end
            if rep ==1 & mouseInterp(i)<1
                movntrial(f(i))=movntrial(f(i))+1;
                movRcycle(mod(f(i),cycframes)+1) = movRcycle(mod(f(i),cycframes)+1) + sum(times>t(i) & times<t(i)+dt_avg);
                movcyctrials(mod(f(i),cycframes)+1) = movcyctrials(mod(f(i),cycframes)+1) +1;
            elseif rep==2 & mouseInterp(i)>1;
                movntrial(f(i))=movntrial(f(i))+1;
                movRcycle(mod(f(i),cycframes)+1) = movRcycle(mod(f(i),cycframes)+1) + sum(times>t(i) & times<t(i)+dt_avg);
                movcyctrials(mod(f(i),cycframes)+1) = movcyctrials(mod(f(i),cycframes)+1) +1;
            end
        end
        R = R./ntrial;
        %         figure
        %         plot(R)
        Rcyclerep{cell_n,rep} = Rcycle./cyctrials;
        movRcyclerep{cell_n,rep} = movRcycle./movcyctrials;
        %         figure
        %         plot(Rcycle)
        if rep ==1 & uselaser
            Rtc = Rtimecycle./timecyc_trials;
            figure
            imagesc(Rtc);
            figure
            plot(Rtc);
            figure
            plot(mean(Rtc(:,[1 10]),2));
            hold on
            plot(mean(Rtc(:,4:6),2));
            RtcAll{ch}= Rtc;
            
            figure
            barwitherr( (1/dt_avg)*[std(nanmean(Rtc(5:9,[1 10]),2)) std(nanmean(Rtc(20:25,[1 10]),2))   ; std(nanmean(Rtc(5:9,[4:6]),2)) std(nanmean(Rtc(20:25,[4:6]),2))], ...
                (1/dt_avg)*[mean(nanmean(Rtc(5:9,[1 10]))) mean(nanmean(Rtc(20:25,[1 10])))   ; ...
                mean(nanmean(Rtc(5:9,[4:6])))-mean(nanmean(Rtc(5:9,[1 10])))  mean(nanmean(Rtc(20:25,[4:6]))) - mean(nanmean(Rtc(20:25,[1 10]))) ]);
            title(sprintf('%s chan %d',Block_Name,ch));
            xlabel(['spont               evoked'])
            
            
        end
        
        %    figure
        %    plot(1-cos(2*pi*(1:300)/300),Rcycle)
    end
    if uselaser
        figure;
        d = condenseData(Rcyclerep{cell_n,1}',15)/(frame_duration);
        d = (d(1:10) + d(20:-1:11))/2;
        plot(0.5*(1-cos(2*pi*(1:10)/20)),d);
        hold on
        d = condenseData(Rcyclerep{cell_n,2}',15)/(frame_duration);
        d = (d(1:10) + d(20:-1:11))/2;
        plot(0.5*(1-cos(2*pi*(1:10)/20)),d,'g');
        legend('stationary','moving');
        xlabel('contrast');
        ylabel('sp/sec');
        title('laser effect');
    end
    
    figure
    d = condenseData(movRcyclerep{cell_n,1}',15)/(frame_duration);
    d = (d(1:10) + d(20:-1:11))/2;
    plot(0.5*(1-cos(2*pi*(1:10)/20)),d);
    hold on
    d = condenseData(movRcyclerep{cell_n,2}',15)/(frame_duration);
    d = (d(1:10) + d(20:-1:11))/2;
    plot(0.5*(1-cos(2*pi*(1:10)/20)),d,'g');
    legend('stationary','moving');
    xlabel('contrast');
    ylabel('sp/sec');
    title(sprintf('movement effect ch=%d',cell_n))
end

if SU
    save(fullfile(apname,afname),'Rcyclerep','movRcyclerep','RtcAll','-append');
end




