%%% code to read in accelerometer data (from labjack) and sycn with video data (from bonsai)
close all
clear all

%%% load acc data
[f,p] = uigetfile('*.dat','acc data');
acc = dlmread(fullfile(p,f));

%%% get timestamps (column 1)
accTs = acc(:,1);
accTs = mod(accTs-7*60*60,24*60*60); %%% time is elapsed secs since midnight 1904 GMT; subtract 7 hrs to get local time (but what about daylight savings change!)

figure
plot(diff(accTs))
ylabel('secs'); title(sprintf('diff of acc timestamps; median = %0.3f',median(diff(accTs))));


%%% get data (columns 2-7)
accTrace = acc(:,2:7);
figure
plot(accTs-accTs(1),accTrace);
xlabel('secs'); ylim([0 5]);
%%
%%% get video
% [f p] = uigetfile('*.avi','video');
% topMov = VideoReader(fullfile(p,f));
% nframe = round(topMov.duration*topMov.framerate)
% nframe = min(nframe,1800) %%% don't read too many - slow!
% frm = 0;
% display('reading video data');
% 
% for frm = 1:nframe
%     if round((frm-1)/25) == (frm-1)/25
%         display(sprintf('done %d/%d frames',frm-1,nframe))
%     end
%     mov(:,:,:,frm)= topMov.readFrame;
% end
%%
%%% get movie timestampes
[f p] = uigetfile('*.csv','video timestamps');
movTs = dlmread(fullfile(p,f));
movTs = movTs(:,1)*60*60 + movTs(:,2)*60 + movTs(:,3); %%% timestamps is hrs, mins, secs
movTs = movTs(1:frm);

%%% interpolate acc data to match movie timestamps
clear accResamp
for i = 1:6
    accResamp(:,i) = interp1(accTs,accTrace(:,i),movTs);
end

%%% plot resampled data
figure
plot(movTs,accResamp);
xlabel('movie timestamp (sec)');
ylim([0 5]); ylabel('accelerometers')
%%
% %%% make movie with video and resampled data
% display('making movie')
% figure
% for i = 1:frm
%     if round((i-1)/25) == (i-1)/25
%         display(sprintf('done %d/%d frames',i-1,frm))
%     end
%     
%     %%% show video
%     subplot(3,3,1:3);
%     imshow(mov(:,:,:,i));
%     
%     %%% plot first three channels
%     subplot(3,3,4:6);
%     plot(movTs(1:i)-movTs(1),accResamp(1:i,1:3));
%     xlim([0 movTs(end)-movTs(1)]);ylim([0 5]);  ylabel('3-d acc');
%     
%     %%% plot second three channels
%     subplot(3,3,7:9)
%     plot(movTs(1:i)-movTs(1),accResamp(1:i,4:6));
%     xlim([0 movTs(end)-movTs(1)]); ylim([0 5]); ylabel('gyros');xlabel('secs');
%     
%     %%% get movie frame
%     drawnow
%     outMov(i) = getframe(gcf);
% end
% 
% %%% write movie
% [f p] = uiputfile('*.avi','video out');
% outFileObj = VideoWriter(fullfile(p,f));
% open(outFileObj);
% writeVideo(outFileObj,outMov);
% close(outFileObj)
