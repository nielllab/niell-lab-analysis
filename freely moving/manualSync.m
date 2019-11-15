%%% code to manually input ttl values based on observing light in videos
%%% correction for the fact that TTL didn't record properly in some prey capture optogenetic extps
%%% cmn 2019

close all; clear all
pathname = {'D:\gdrive\cohort5_100919\Cohort5\timestamp info\'};

fileList = [];
for i = 1:length(pathname)
    fileList = [fileList ; dir([pathname{i} '*.csv'])];
end

for f = 5:length(fileList);

 sprintf(fileList(f).name)
 thispath = fileList(f).folder;
 
     movName = strrep(fileList(f).name,'csv','avi');
    topMov = VideoReader(fullfile(thispath,movName));
    nframe = floor(topMov.duration*topMov.framerate)

ttlSync = zeros(nframe,1);
 done = 0;n = 0;
 while ~done
    
     n=n+1;
     frame(n)= input('frame # (-1 to quit) :');
     if frame(n)==nframe;
         done = 1;
     end
         val(n) = input('ttl value : ');
 end
 for i = 1:n-1;
     ttlSync(frame(i):frame(i+1)) = val(i);
 end
 
 figure
 plot(ttlSync);
    matName = strrep(movName,'.avi','_sync.mat');
    save(fullfile(thispath,matName),'ttlSync');
  
end
