% figure('units','normalized','outerposition',[0 0 1 .5])
% figure('visible','off');
figure;
axis tight manual
% ax1 = subplot(1,3,1);
ax = gca;
ax.NextPlot = 'replaceChildren';
% ax2 = subplot(1,3,2);
% ax2.NextPlot = 'replaceChildren';
% ax3=subplot(1,3,2);
% ax3.NextPlot = 'replaceChildren';

% F(size(EyeVid2,3)) = struct('cdata',[],'colormap',[]);
vidfile = VideoWriter('PupilPoints.mp4','MPEG-4');
% vidfile.FrameRate=10;
open(vidfile);
for  v=1:size(EyeVid2,3)
%      ax = subplot(1,3,1);
    imagesc(EyeVid2(:,:,v)); colormap gray; hold on; axis equal off;
    scatter(Pointsx(v,:),Pointsy(v,:),100,'.r'); hold off;
%     e_t = fit_ellipse(Pointsx(v,:),Pointsy(v,:),ax);
%     subplot(1,3,2)
%     imshow(FieldView(:,:,v));colormap gray; hold on; axis equal off;
%     subplot(1,3,3)
%     imshow(TopView(:,:,v));colormap gray; hold on; axis equal off;
    
    drawnow limitrate;
    F(v) = getframe(gcf); 
    writeVideo(vidfile,F(v));
    fprintf('Frame = %d \n',v);
end
close(vidfile)
%%
fig = figure('visible','off','units','normalized','outerposition',[0 0 1 .5]);

% figure;
axis tight manual
% ax1 = subplot(1,3,1);
% % ax = gca;
% ax1.NextPlot = 'replaceChildren';
% ax2 = subplot(1,3,2);
% ax2.NextPlot = 'replaceChildren';
% ax3=subplot(1,3,2);
% ax3.NextPlot = 'replaceChildren';

% F(size(EyeVid2,3)) = struct('cdata',[],'colormap',[]);
% vidfile = VideoWriter('PupilWorld.mp4','MPEG-4');
for v=1:size(Pointsx,1)
%     temp = EyeVid2(:,:,v);
    ax = subplot(1,3,1);
    imagesc(EyeVid2(:,:,v)); colormap gray; hold on; axis equal off;
    e_t = fit_ellipse(Pointsx(v,:),Pointsy(v,:),ax);
%     title(sprintf('Frame = %d',v));
    
    subplot(1,3,2)
    imshow(FieldView(:,:,v));colormap gray; hold on; axis equal off;
    subplot(1,3,3)
    imshow(TopView(:,:,v));colormap gray; hold on; axis equal off;
%     title(sprintf('Frame = %d',v));
    drawnow limitrate
%     F(v) = getframe(gcf);
    fprintf('Frame = %d \n',v);
end

clear myVideo;
myVideo = VideoWriter('PupilWorld.avi');
myVideo.FrameRate = 30;  % Default 30
% myVideo.Quality = 50;    % Default 75
open(myVideo);
writeVideo(myVideo, F);
close(myVideo)
fprintf('Done \n');
%%
%% Make Eye Top and World Videos

figure;%('units','normalized','outerposition',[0 0 1 1]); 
% subplot(1,2,1)
ax = gca;
ax.NextPlot = 'replaceChildren';
axis([1 size(EyeVid,1) 1 size(EyeVid,2)] )
% EllipseParams = zeros(size(EyeVid,3),6);
ptxtemp=[]; ptytemp=[];j=1; ptxtemp2 = []; ptytemp2=[]; k = 1;
% F(size(EyeVid2,3)) = struct('cdata',[],'colormap',[]);
clear F
int= 3000:size(EyeVid,3);
F(size(int,2)) = struct('cdata',[],'colormap',[]); frm = 1;

for v=int%size(EyeVid,3)
%     subplot(1,2,1)
%     imshow(EyeVid2(:,:,v)); colormap gray; axis equal off; hold on;    
%     if ~any((Pointsx(v,:)<10 | Pointsy(v,:)<10) | isnan(Pointsx(v,:))) % check points
%         e_t = fit_ellipse(Pointsx(v,:),Pointsy(v,:),ax);
%     end
%     hold off;
%     plot_ellipse(ExtraParams(v,:),ax);
%     ptxtemp2 = [ptxtemp2; EllipseParams(v,1)];
%     ptytemp2 = [ptytemp2; EllipseParams(v,2)];
%     plot(ptxtemp2(1:k,:),ptytemp2(1:k,:),'.-b'); k = k+1;
%     if size(ptxtemp2,1) >= 60
%         ptxtemp2(1,:) = [];
%         ptytemp2(1,:) = [];
%         k = k-1;
%     end
%     hold off;
    
%%%%% Top Camera
%     subplot(1,2,2)
    imshow(TopVid2(:,:,v)); colormap gray; axis equal off; hold on;
    scatter(PointsxT(v,1),PointsyT(v,1),100,'.g')
    scatter(PointsxT(v,2),PointsyT(v,2),100,'.r')
    scatter(PointsxT(v,3),PointsyT(v,3),100,'.b'); hold off
%     % Trail of points
%     ptxtemp = [ptxtemp; PointsxT(v,:)];
%     ptytemp = [ptytemp; PointsyT(v,:)];
%     plot(ptxtemp(1:j,1),ptytemp(1:j,1),'.-g'); 
%     plot(ptxtemp(1:j,2),ptytemp(1:j,2),'.-r');
%     plot(ptxtemp(1:j,3),ptytemp(1:j,3),'.-b'); j =j+1;    hold off;
%     if size(ptxtemp,1) >= 60
%         ptxtemp(1,:) = [];
%         ptytemp(1,:) = [];
%         j = j-1;
%     end

%%%%% World Camera
%     imshow(WorldVid(:,:,v)); colormap gray; axis equal off; hold on;
%     scatter(size(WorldVid,1)/2+shiftxW(v,1),size(WorldVid,2)/2-shiftyW(v,1),1000,'or','LineWidth',3)
%     scatter(size(WorldVid,1)/2+shiftxW(v,1),size(WorldVid,2)/2-shiftyW(v,1),500,'.r','LineWidth',3)
%     hold off;
%     title(sprintf('Frame = %d',v));
    drawnow limitrate
    F(frm) = getframe(gcf);
    frm=frm+1;
%     pause(.5)
end
clear myVideo;
myVideo = VideoWriter('Top_New4.avi');
myVideo.FrameRate = 30;  % Default 30
open(myVideo);
writeVideo(myVideo, F);
close(myVideo)
fprintf('Done \n');