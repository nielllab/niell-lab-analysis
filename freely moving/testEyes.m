sprintf('subj %s session %s date %s clipnum %s',Data(vid).ani,Data(vid).sesssionnum, Data(vid).date,Data(vid).clipnum)

[f p] = uigetfile('*.avi','video file left');

[f p] = uigetfile('*.avi','video file right');

figure
for i= 1:nframe
    subplot(2,2,1);
    imshow(leftmov(:,:,:,i))
    
    subplot(2,2,2)
    imshow(rightmov(:,:,:,i))
    
    subplot(2,2,3)
    plot(Data(vid).Lthetaraw(1:i),plot(Data(vid).Lphiraw(1:i)));
    hold on;  plot(Data(vid).Lthetaraw(i),plot(Data(vid).Lphiraw(1:i)),'o');
    axis([-90 90 -90 90])
    
   subplot(2,2,4)
    plot(Data(vid).Rthetaraw(1:i),plot(Data(vid).Rphiraw(1:i)));
    hold on;  plot(Data(vid).Rthetaraw(i),plot(Data(vid).Rphiraw(1:i)),'o');
    axis([-90 90 -90 90])
    
end