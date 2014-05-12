function pptgSpeedAnalysis(clustfile,afile,pdfFile,Block_Name,blocknum,uselaser, framerate);

close all
dbstop if error

if exist('clustfile','var');
    SU=1;
    
    load(clustfile);
    load(afile);
    Block_Name = Block_Name{blocknum}
    block=blocknum;
    if max(cells(:,1))>32
        nchan=64;
    else
        nchan=32;
    end
    chans = 1:4:nchan;
    use_afile=1;
    useArgin=1;
else
    useArgin=0;
    nchan = input('# chans : ');
    uselaser = input('laser on 0/1 : ');
    chans = 1:4:nchan;
    SU = input('multiunit (0) or single-unit (1) : ');
    framerate = input('movie framerate : ');
    
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
        pname = uigetdir('D:\','block data')
        delims = strfind(pname,'\');
        selected_path = pname(1 :delims(length(delims))-1)
        Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
        Block_Name = pname(delims(length(delims))+1 :length(pname))
    end
    
    
end


if uselaser
    flags = struct('mouseOn',1,'laserOn',1)
else
    if SU
        flags = struct('mouseOn',1,'visStim',1,'lfpTseries',1,'lfpSpectra',1);
    else
        flags = struct('lfpTseries',1,'lfpSpectra',1,'mouseOn',1,'MUspike',1,'visStim',1);
    end
end

tdtData= getTDTdata(Tank_Name, Block_Name, chans, flags);
%tdtData= getTDTdata(Tank_Name, Block_Name, nChan, 0, 0, 0, 1,1, 0, 0, 0);

mainfig= figure;
if uselaser
    nx=3; ny=2;
else
    nx=2;ny=1;
end
subplot(nx,ny,1);
plot(tdtData.mouseT,tdtData.mouseV); hold on;
xlabel('secs');
ylabel('cm/sec')

if ~SU
    bins = 0:0.25:max(tdtData.mouseT);
    R= zeros(length(bins),length(chans));
    for ch = chans;
        R(:,ch)= hist(tdtData.MUspikeT{ch},bins);
    end
    subplot(nx,ny,2)
    plot(bins(1:end-1),R(1:end-1,:))
    xlabel('secs')
    ylabel('sp/sec')
    xlim([0 max(tdtData.MUspikeT{ch})]);
end

if uselaser
    smoothwindow_secs = 0.25;
    laserT = tdtData.laserT;
    dt = laserT(2)-laserT(1);
    smoothwindow = ceil(smoothwindow_secs/dt)
    lasersmooth = zeros(size(tdtData.laserTTL));
    %%% replace this with convolution (for speed)!
    %     for i = smoothwindow+1:length(tdtData.laserTTL);
    %         lasersmooth(i)= mean(tdtData.laserTTL(i-smoothwindow:i));
    %     end
    smoothfilter=ones(smoothwindow,1);
    lasersmooth=conv(double(tdtData.laserTTL),smoothfilter,'same')/sum(smoothfilter);
    
    lasersmooth = lasersmooth;
    figure(mainfig)
    subplot(nx,ny,1)
    plot(laserT,lasersmooth,'g')
    xlim([0 max(laserT)])
    
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
    figure(mainfig)
    subplot(nx,ny,3);
    plot(timeRange,vdata)
    hold on
    plot([0 17],[40 40],'g','LineWidth',12)
    xlim([min(timeRange) max(timeRange)])
    xlabel('secs')
    ylabel('cm/sec');
    
    subplot(nx,ny,4);
    plot(timeRange,nanmean(vdata,2))
    hold on
    plot([0 17],[20 20],'g','LineWidth',12)
    xlim([min(timeRange) max(timeRange)])
    xlabel('secs')
    ylabel('avg cm/sec');
    tsamps = 0:0.5:max(laserT);
    laserdownsamp = interp1(laserT,lasersmooth,tsamps);
    
    for laser =1:2
        if laser ==1
            timepts = tsamps(laserdownsamp==0);
        else
            timepts = tsamps(laserdownsamp>0);
        end
        speeds = interp1(tdtData.mouseT,tdtData.mouseV,timepts);
        laserspeed(laser) = nanmean(speeds);
        laserspeed_std(laser) = nanstd(speeds)/sqrt(npulse);
    end
    subplot(nx,ny,6);
    barweb(laserspeed,laserspeed_std);
    legend({'off','on'})
    ylabel('cm/sec')
    
    
end


save(afile,'vdata','-append');






