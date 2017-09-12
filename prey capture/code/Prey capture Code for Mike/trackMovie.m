function mov = trackMovie(vid,head,targ,theta,frames,bin,figName);

clear mov
figure
set(gcf,'Name',figName)
for f =1:length(frames)-1;
    %f=23:65;session 9 of light sessions mouse prey capture
     %figure
    i = frames(f);
    hold off
    if ~isempty(vid)
        imshow(vid(:,:,i))
        hold on
    end
 
%     plot(targ(i*bin,1)/bin,targ(i*bin,2)/bin,'ro','Markersize',12)
%     hold on
%     plot(targ(frames(1)*bin:i*bin,1)/bin,targ(frames(1)*bin:i*bin,2)/bin,'k','LineWidth',2)
%     
    plot(head(i*bin,1)/bin,head(i*bin,2)/bin,'bo','Linewidth',2);
    plot(head(frames(1)*bin:i*bin,1)/bin,head(frames(1)*bin:i*bin,2)/bin,'c','LineWidth',2);
    
     rangeangle = atan2( targ(i*bin,2) - head(i*bin,2), targ(i*bin,1) - head(i*bin,1))*180/pi;
     headangle = rangeangle -theta(i*bin);
     plot([head(i*bin,1) head(i*bin,1) + cosd(headangle)*50 ]/bin, [head(i*bin,2) head(i*bin,2) + sind(headangle)*50]/bin,'b','Linewidth',2);
    
     axis([0 2050 0 1080]/bin);
     mov(f)=getframe(gcf);
     %save out each frame of tracked movie
     %saveas(gca,['Chase_seq' num2str(i) '.eps']) 
end

%saveas(gca,['Chase_seq' num2str(i) '.eps']) 

