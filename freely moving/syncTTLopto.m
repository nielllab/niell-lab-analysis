close all; clear all
pathname = {'D:\gdrive\cohort5_100919\Cohort5\timestamp info\'};

fileList = [];
for i = 1:length(pathname)
    fileList = [fileList ; dir([pathname{i} '*.csv'])];
end

for f = 1:length(fileList);
    topcam_fname = fileList(f).name;
    ttl_fname = strrep(topcam_fname,'top','opto');
    ttl_fname = strrep(ttl_fname,'csv','dat');
    thispath = fileList(f).folder;
    
    TopTs = dlmread(fullfile(thispath,topcam_fname));
    TopTs= TopTs(:,1)*60*60 + TopTs(:,2)*60 + TopTs(:,3);
    
    ttl = dlmread(fullfile(thispath,ttl_fname));
    
    
    %%% get timestamps (column 1)
    ttlTs = ttl(:,1);
    ttlTs = mod(ttlTs-7*60*60,24*60*60); %%% time is elapsed secs since midnight 1904 GMT; subtract 7 hrs to get local time (but what about daylight savings change!)
    
    ttlSig = ttl(:,2);
    win = 150;
    ttlSigFilt = double(ttlSig<-0.1);
    ttlOn = conv(ttlSigFilt,ones(win,1)/win,'same');
    ttlOn(ttlOn>0)=1;
    ttlOn(1:win) = ttlOn(win+1);
    ttlOn(end:end-win) = ttlOn(end-(win+1));
    
    figure
    plot(ttlTs - ttlTs(1),ttlSig);hold on
    plot(ttlTs - ttlTs(1),ttlSigFilt/10);
    plot(ttlTs-ttlTs(1),ttlOn/10,'g');
    title(ttl_fname);
    
    ttlSync = interp1(ttlTs,ttlOn,TopTs);
    figure
    plot(TopTs - TopTs(1),ttlSync)
    title(ttl_fname);
    
    
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

end
