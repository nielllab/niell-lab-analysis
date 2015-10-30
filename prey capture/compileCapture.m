%close all
clear all

captureBatch

for f = 1:length(files)
    lighting(f) = files(f).lighting;
    contrast(f) = files(f).contrast;
    
    [xc(:,f) vxc(:,f) location(:,:,f) r{f}] = analyzeCapture([pathname files(f).trackpts]) ;
    latency(f) = length(r{f});
    dist(f) = mean(r{f});
end

figure
subplot(2,2,1);
imagesc(flipud(mean(location(:,:,lighting==0),3)),[0 0.05]); hold on; plot(21,21,'ro','Markersize',4)
title('dark')
subplot(2,2,3);
imagesc(flipud(mean(location(:,:,lighting==1),3)),[0 0.05]);hold on; plot(21,21,'ro','Markersize',4)
title('light')

subplot(2,2,2);
plot((-300:300)/60,mean(xc(:,lighting==0),2)); axis([-5 5 -1 1])
hold on; plot([0 0],[-1 1],'r:')
plot((-300:300)/60,mean(xc(:,lighting==1),2),'g');axis([-5 5 -1 1])
title('angle corr'); xlabel('secs')

subplot(2,2,4);
plot((-120:120)/60,mean(vxc(:,lighting==0),2)); axis([-2 2 -0.25 0.5])
hold on; plot([0 0],[-1 1],'r:')
plot((-120:120)/60,mean(vxc(:,lighting==1),2),'g');axis([-2 2 -0.25 0.5])
title('speed corr'); xlabel('secs')

sprintf('mean time to capture')
sprintf('light = %f sec',mean(latency(lighting==1))/60)
sprintf('dark = %f sec',mean(latency(lighting==0))/60)


