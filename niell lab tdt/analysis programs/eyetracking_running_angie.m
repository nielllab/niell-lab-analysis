clear all
close all

dir = 'D:\Angie_analysis\DOI_experiments\12_08_16_eye';
% name = '120516_driftlowcontrast_eye'; %data file
name = {'12_08_16_drift_highcontrast_pre_eye',
         '12_08_16_drift_highcontrast_post_eye'};      
thresh = 0.80; %pupil threshold for binarization
puprange = [15 60];

%loop overwrites first data set...need to fix
%for i=1:length(name)
load(name{2})
data = squeeze(data);
warning off;

Tank_Name = '12_08_16_eye'
%if name{i} ==1
Block_Name = 'drift_highcontrast_post_eye'
% else
% Block_Name = 'dark_post_eye' %end;

flags = struct('mouseOn',1,'cameraOn',1);
tdtData= getTDTdata(Tank_Name, Block_Name, 1, flags);
tsampBar = tdtData.mouseT;
vsmoothBar = tdtData.mouseV;

cameraTTL = tdtData.cameraT;
figure
plot(diff(cameraTTL));
xlabel('frame #'); ylabel('secs');
ylim([ 0 0.2])

thresh = 0.8; %pupil threshold for binarization
puprange = [8 45]; %set

%user input to select center and right points
sprintf('Please select pupil center and top, eyeball top and right points, darkest part of eyeball')
h1 = figure('units','normalized','outerposition',[0 0 1 1])
imshow(data(:,:,10))
[cent] = ginput(5);
close(h1);
yc = cent(1,2); %pupil center y val
xc = cent(1,1); %pupil center x val
horiz = (cent(4,1) - xc); %1/2 x search range
vert = (yc - cent(3,2)); %1/2 y search range
puprad = yc - cent(2,2); %initial pupil radius
% puprange = [round(puprad - puprad*pupercent) round(puprad + puprad*pupercent)]; %range of pupil sizes to search over
ddata = double(data); 
binmaxx = cent(5,1);
binmaxy = cent(5,2);

for c = 1:size(data,3)
    binmax(c) = mean(mean(ddata(binmaxy-3:binmaxy+3,binmaxx-3:binmaxx+3,c)));
end
%end
clear bindata
for d = 1:size(ddata,3)
    bindata(:,:,d) = (ddata(yc-vert:yc+vert,xc-horiz:xc+horiz,d)/binmax(d) > thresh);
end

tic
centroid = nan(size(data,3),2);
rad = nan(size(data,3),1);
centroid(1,:) = [horiz vert];
rad(1,1) = puprad;
for n = 2:size(data,3)
    [center,radii,metric] = imfindcircles(bindata(:,:,n),puprange,'Sensitivity',0.995,'ObjectPolarity','dark');
    if(isempty(center))
        centroid(n,:) = [NaN NaN]; % could not find anything...
        rad(n) = NaN;
    else
        [~,idx] = max(metric); % pick the circle with best score
        centroid(n,:) = center(idx,:);
        rad(n,:) = radii(idx);
    end
end

%end
toc

startT = max(tsampBar(1),cameraTTL(1))
endT = min(tsampBar(end),cameraTTL(end))

dt = 0.5;
pts = startT:dt:endT;
vInterp = interp1(tsampBar,vsmoothBar,pts);
%rInterp = interp1(cameraTTL(1:length(rad')),rad,pts); %if cameraTTL is longer than rad (due to blinks?)
 rInterp = interp1(cameraTTL,rad(1:length(cameraTTL)),pts); %if rad is longer
figure
plot(vInterp,'g'); hold on; plot(rInterp)
legend('velocity','radius')

% h2 = figure
% hold on
% plot(0.1:0.1:size(data,3)/10,rad,'g-')
% plot(0.1:0.1:size(data,3)/10,centroid(:,1),'.b')
% plot(0.1:0.1:size(data,3)/10,centroid(:,2),'.m')
%hold off
% legend('radius','x pos','ypos')
% ylim([5 55])

% %plot cross correlation of velocity and radius
[corr_vrad lags] = xcorr((~isnan(rInterp-nanmean(rInterp))),(~isnan(vInterp-nanmean(vInterp))),60/dt,'coeff')
figure
plot(lags*dt,corr_vrad);
hold on
plot([0 0],[0 1],'g-')

%radius as a function of speed
figure
scatter(rInterp, vInterp, '.')
xlabel('radius');ylabel('cm/sec')
lsline;

% video with tracking
h4 = figure
for i = 1:size(data,3)
 
    subplot(1,2,1)
    imshow(data(yc-vert:yc+vert,xc-horiz:xc+horiz,i));
    colormap gray
    hold on
    circle(centroid(i,1),centroid(i,2),rad(i))
    drawnow
    hold off
    
    subplot(1,2,2)
    imshow(bindata(:,:,i));
    colormap gray
    hold on
    circle(centroid(i,1),centroid(i,2),rad(i))
    drawnow
    hold off 
%     mov(i) = getframe(gcf)
%     vid = VideoWriter('predoi_tracking.avi')
%     open(vid); writeVideo(vid,mov); close(vid)
end

