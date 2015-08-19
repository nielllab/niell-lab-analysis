function pptgAnalysisLFP(clustfile,afile,pdfFile,Block_Name,blocknum,uselaser, framerate);

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
    uselaser = 1
    chans = 1:4:nchan;
    
    pname = uigetdir('D:\','block data')
    delims = strfind(pname,'\');
    selected_path = pname(1 :delims(length(delims))-1)
    Tank_Name = pname(delims(length(delims)-1)+1 :delims(length(delims))-1)
    Block_Name = pname(delims(length(delims))+1 :length(pname))
end


    flags = struct('lfpTseries',1,'lfpSpectra',1','laserOn',1,'MUspike',1,'visStim',1)


tdtData= getTDTdata(Tank_Name, Block_Name, chans, flags);
%tdtData= getTDTdata(Tank_Name, Block_Name, nChan, 0, 0, 0, 1,1, 0, 0, 0);

mainfig= figure;
nx=2; ny =2;


if uselaser
    smoothwindow_secs = 0.25;
    laserT = tdtData.laserT;
    laserRaw = tdtData.laserTTL;
    laserdt = median(diff(laserT));
    
  
    
    lfp = interp1(tdtData.lfpT{1},tdtData.lfpData{1},laserT);
      lfptimes = laserT;
    lfpdt = median(diff(lfptimes));
    
    rising = find(diff(laserRaw)>1)+1;
   trainstart = rising(find(diff(laserT(rising))>10)+1);
   
%    figure
%    plot(laserT(rising),1,'*')
%    title('all pulse start')
%    
%    figure
%    plot(diff(laserT(rising)),'*')
%    title('pre-pulse interval')
   
   laserT(trainstart)
   trainstart = trainstart(laserT(trainstart)>5 & laserT(trainstart)<(max(laserT-30)));
    rising= rising(laserT(rising)>5 & laserT(rising)<(max(laserT-30)));
    laserOn = conv(double(laserRaw),ones(ceil(0.001/laserdt),1),'same');
   lfpinterp = interp1(lfptimes(laserOn==0),lfp(laserOn==0),lfptimes);
   
   clear lfplocked lfpInterpLocked laserLock
   for i = 1:length(trainstart)
       startlfp = min(find(lfptimes > laserT(trainstart(i))-5));
       lfplocked(i,:) = lfp(startlfp : startlfp + round(30/lfpdt));
       lfpInterpLocked(i,:) = lfpinterp(startlfp : startlfp + round(30/lfpdt));
       laserLock(i,:) = laserRaw(startlfp : startlfp + round(30/lfpdt));
   end
   
figure(mainfig)
 subplot(2,2,3)
   
   mnLfpLocked = mean(lfplocked,1);
   mnLfpInterp = mean(lfpInterpLocked,1);
   plot((0:round(30/lfpdt))*lfpdt-5, mnLfpLocked);
   hold on
  % plot((0:round(30/lfpdt))*lfpdt-5, mnLfpInterp,'g');
   plot( (0:round(30/lfpdt))*lfpdt-5, mean(laserLock,1)*10,'g')
   xlim([-0.5 1])
   xlabel('secs'); %legend('lfp','laser')
   

    lfplocked = zeros(length(rising),1+round(0.5/lfpdt));
    lfpInterpLocked = lfplocked; laserLock = lfplocked;
   for i = 1:length(rising)
     % startlfp = min(find(lfptimes > laserT(rising(i))-0.1));
      startlfp = rising(i)-round(0.1/lfpdt);
      lfplocked(i,:) = lfp(startlfp : startlfp + round(0.5/lfpdt));
       lfpInterpLocked(i,:) = lfpinterp(startlfp : startlfp + round(0.5/lfpdt));
       laserLock(i,:) = laserRaw(startlfp : startlfp + round(0.5/lfpdt));
   end
   
subplot(2,2,4)
   
   mnLfpLocked = mean(lfplocked,1);
   mnLfpInterp = mean(lfpInterpLocked,1);
   plot(((1:size(lfplocked,2))*lfpdt-0.1)*1000, mnLfpLocked);
   hold on
  % plot((0:round(30/lfpdt))*lfpdt-5, mnLfpInterp,'g');
   plot( ((1:size(lfplocked,2))*lfpdt-0.1)*1000, mean(laserLock,1)*10,'g')
   xlim(1000*[-0.005 median(diff(laserT(rising)))])
   xlabel('msecs'); %legend('lfp','laser')
   
    
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
%     figure(mainfig)
%     subplot(nx,ny,1)
%     plot(laserT,lasersmooth,'g')
%     xlim([0 max(laserT)])
    
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


    tsamps = 0:0.5:max(laserT);
    laserdownsamp = interp1(laserT,lasersmooth,tsamps);
    
    for laser =1:2
        if laser ==1
            timepts = tsamps(laserdownsamp==0);
        else
            timepts = tsamps(laserdownsamp>0);
        end

    end
 

    for ch = chans;
        freqs{ch} = tdtData.spectF{ch};
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
        figure(mainfig)
        subplot(2,2,1)
        imagesc(spect_avg',1.5*[0 prctile(spect_avg(:),95)])
        axis xy
        df = median(diff(tdtData.spectF{ch}));
        dt = median(diff(timeRange));
       ylim([0 100/df])
       set(gca,'YTick',(10:10:80)/df);
        set(gca,'YTickLabel',{'10','20','30','40','50','60','70','80'})
        set(gca,'XTick',(5:5:40)/dt);
        set(gca,'XTickLabel',{'-5','0','5','10','15','20','25','30'})
        title(sprintf('ch=%d',ch));
        figure(mainfig)
        
        
        %%% selects timepoints for laser on/off and running speed

        for laser = 1:2
            
            if laser ==1
                timepts = tsamps(laserdownsamp==0 );
            elseif laser ==2
                timepts = tsamps(laserdownsamp>0 );
            end
            
            for f = 1:size(specdata,2)
                laserlfp(ch,laser,f)= nanmean(interp1(tdata,squeeze(specdata(:,f)),timepts));
            end
        end
        subplot(2,2,2)
        plot(tdtData.spectF{ch}, squeeze(laserlfp(ch,:,:)));
        ylim([0 1.5*prctile(laserlfp(ch,1,:),95)])
        xlim([0 90])
        
    end
   % legend({'laser off ','laser on '})
end



% if useArgin
%         psfname = [pdfFile(1:end-4) Block_Name '.ps'];
%     else
%         [fname pname] =uiputfile('*.ps'); psfname=fullfile(pname,fname);
% end
%     if exist(psfname,'file')==2;delete(psfname);end


% figure(mainfig)
% set(gcf, 'PaperPositionMode', 'auto');
% print('-dpsc',psfname,'-append');
%
% for i = 1:length(spikefig)
%     figure(spikefig(i))
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfname,'-append');
%
%     figure(burstfig(i))
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfname,'-append');
%
%     %       figure(lfpfig(i))
%     %     set(gcf, 'PaperPositionMode', 'auto');
%     %     print('-dpsc',psfname,'-append');
%
% end
%
% if uselaser
%     figure(spectrafig)
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfname,'-append');
%
%     figure(lfpfig(1))
%     set(gcf, 'PaperPositionMode', 'auto');
%     print('-dpsc',psfname,'-append');
%
%
% end



end






