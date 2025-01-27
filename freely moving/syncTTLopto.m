close all; clear all
pathname = {'./'};

fileList = []; paths = {};
n=0;
for i = 1:length(pathname)
    
    fs = dir([pathname{i} '*.dat']);
    for j= 1:length(fs);
        n=n+1;
        paths{n} =pathname{i};
    end
    fileList = [fileList ; fs];
end

for f =1:length(fileList);
    ttl_fname = fileList(f).name;
    topcam_fname = strrep(ttl_fname,'opto','top');
    topcam_fname = strrep(topcam_fname,'dat','csv');
    thispath = paths{f}
    
try
        
    TopTs = dlmread(fullfile(paths{f},topcam_fname));
    TopTs= TopTs(:,1)*60*60 + TopTs(:,2)*60 + TopTs(:,3);
    
    ttl = dlmread(fullfile(paths{f},ttl_fname));
    
    
    %%% get timestamps (column 1)
    ttlTs = ttl(:,1);
    ttlTs = mod(ttlTs-7*60*60,24*60*60); %%% time is elapsed secs since midnight 1904 GMT; subtract 7 hrs to get local time (but what about daylight savings change!)
    
    ttlSig = ttl(:,2);
    win = 15;
   
    ttlSigFilt = conv(ttlSig,ones(win,1)/win,'same');
    ttlOn= double(ttlSigFilt>0.25);
    ttlOn(1:win) = ttlOn(win+1);
    ttlOn(end:end-win) = ttlOn(end-(win+1));
    
%     figure
%     plot(ttlTs - ttlTs(1),ttlSig);hold on
%     plot(ttlTs - ttlTs(1),ttlSigFilt/10);
%     plot(ttlTs-ttlTs(1),ttlOn/10,'g');
%     title(ttl_fname);
%     
    ttlSync = interp1(ttlTs,ttlOn,TopTs);
    ttlPre = ttlSync;
    firstgood = min(find(~isnan(ttlSync)));
    if firstgood>1
        ttlSync(1:firstgood) = ttlSync(firstgood);
    end
    lastgood = max(find(~isnan(ttlSync)));
    if lastgood<length(ttlSync);
        ttlSync(lastgood:end)=ttlSync(lastgood);
    end
    
     if sum(~isnan(ttlSync))==0;
         ttlSync(:) = 0;
     end
    
%     figure
%     plot(ttlSync); hold on; plot(ttlPre);  
    
    
%     figure
%     plot(TopTs - TopTs(1),ttlSync)
%     title(ttl_fname);
    
    
    matName = strrep(ttl_fname,'dat','mat');
    save(matName,'ttlSync');
    
%     movName = strrep(topcam_fname,'csv','avi');
%     topMov = VideoReader(fullfile(thispath,movName));
%     nframe = floor(topMov.duration*topMov.framerate)
%     mov = zeros(topMov.Height,topMov.Width,3,nframe,'uint8');   
%     for i = 1:nframe
%         i
%         mov(:,:,:,i) = topMov.readFrame;
%     end
%     
%     nframe = size(mov,4);
%     figure
%     for i = 1:nframe;
%         subplot(3,3,1:6);
%         imshow(mov(:,:,:,i));
%         subplot(3,3,7:9)
%         plot(ttlSync); hold on;
%         plot(i,ttlSync(i),'o'); hold off
%         syncMov(i) = getframe(gcf);
%     end
%     
%     movieF = strrep(movName,'.avi','_sync.avi');
%     
%     outFileObj = VideoWriter(fullfile(thispath,movieF));
%     open(outFileObj);
%     writeVideo(outFileObj,syncMov);
%     close(outFileObj)
    
ttlAll{f} = ttlSync;

    % data = readtable([paths{f} 'states/' topcam_fname]); 
     data = readtable([paths{f} strrep(topcam_fname,'top','labeled')]); 
    times = (data{16:end,1});
        states = data{16:end,6};
    startstop = data{16:end,9};
    vidT = TopTs - TopTs(1);
    
    fps = data{16,4};
    fps = str2num(fps{1});
    for i = 1:length(times)
        stateT(i) = str2num(times{i});
        stateF(i) = round(stateT(i)*fps)+1;
    end
    
    vidDt = median(diff(vidT));
    
    chase = strcmp(states,'pursue/chasing') &strcmp(startstop,'START');
    endChase = strcmp(states,'pursue/chasing') &strcmp(startstop,'STOP');
     capture = strcmp(states,'consuming') &strcmp(startstop,'START');
     if sum(capture)==0
         capture(max(find(endChase)))=1;
     end
     
     if sum(endChase)==0;
         endChase(max(find(capture)))=1;
     end
     
     if sum(chase)==0;
         chase(max(find(endChase)))=1;
     end
     
     ttlSync = ttlSync(1:stateF(find(capture,1,'last')));
     optoDur(1,f) = sum(ttlSync==0)*vidDt;
     optoDur(2,f) = sum(ttlSync==1)*vidDt;
     nChase(1,f) = sum(ttlSync(stateF(chase))==0);
     nChase(2,f) = sum(ttlSync(stateF(chase))==1);
     nCapture(1,f) = sum(ttlSync(stateF(capture))==0);
   nCapture(2,f) = sum(ttlSync(stateF(capture))==1);
   
      figure
      t = (1:length(ttlSync))*vidDt;
    plot(t,ttlSync); ylim([-0.1 1.1]);
    hold on
    plot(stateF(chase)*vidDt,ones(size(stateF(chase)))*0.5,'g*');
    plot(stateF(endChase)*vidDt,ones(size(stateF(endChase)))*0.5,'r*');
    plot(stateF(capture)*vidDt,ones(size(stateF(capture)))*0.6,'b*');
    
    xlim([ 0 max(60, length(ttlSync)*vidDt)])
    
    title(sprintf('%s chase %0.2f %0.2f capture %0.2f %0.2f',ttl_fname,nChase(1,f),nChase(2,f),nCapture(1,f),nCapture(2,f)));
    
   
   catch
     sprintf('couldnt do %s boris',topcam_fname)
end
 
    
end

optoDur = optoDur(:,1:size(nCapture,2));


pChase = nChase./optoDur
pCapture = nCapture./optoDur


totalDur = sum(optoDur,1);
use = totalDur<600; %%% don't take > 2 cycles



display('chase')
%1./nanmean(pChase,2)
nansum(optoDur(:,use),2)./nansum(nChase(:,use),2)


display('capture')
%1./nanmean(pCapture,2)
nansum(optoDur(:,use),2)./nansum(nCapture(:,use),2)



