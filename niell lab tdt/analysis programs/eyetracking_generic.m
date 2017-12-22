%This analysis takes in a .mat file with eyetracking data in the form of a
%3-dimensional array "data(xpixel,ypixel,frame)" output from imaqtool. It
%binarizes the images and fits a circle to the pupil. Depending on the
%image levels and the pupil size in pixels you'll want to adjust the
%threshold and pupil range. Upon running the code, you'll be prompted to
%select points on an image to improve the analysis: 1) pupil center, 2) top
%edge of pupil, 3) top edge of eyeball, 4) right edge of eyeball, 5)
%darkest part of the eyeball. The analysis outputs the centroid (X/Y) of
%the pupil center and the radius of the pupil. Any frame where a circle
%could not be fit will contain NAN values. The code currently commented out
%is intended for closed loop analysis - it uses recent sizes and locations
%to guess the next ones to improve tracking under conditions with more eye
%movement and pupil size changes. It's still a bit buggy.

%This code builds off Scanbox code by Dario
%Ringarch and was developed by P.R.L. Parker and A.M. Michaiel under C.M.
%Niell in 2015.

clear all; close all

dir = 'D:\Angie_analysis\DOI_experiments\12_08_16_eye';
name = {'12_08_16_drift_highcontrast_pre_eye'}; %data file exported from matlab     
load(name)
data = squeeze(data);
warning off;


Tank_Name = '12_08_16_eye' %where ephys/running udp data is saved
Block_Name = 'drift_highcontrast_post_eye'

flags = struct('mouseOn',1,'cameraOn',1);
tdtData= getTDTdata(Tank_Name, Block_Name, 1, flags);
tsampBar = tdtData.mouseT; % time stamps from ephys acquisition
vsmoothBar = tdtData.mouseV; % running velocity

cameraTTL = tdtData.cameraT; %time stamps for TTL pulses sent from TDT to camera
figure
plot(diff(cameraTTL)); %check no dropped TTL pulses...acquisition is 10 hz
xlabel('frame #'); ylabel('secs');
ylim([ 0 0.2])

thresh = 0.8; %pupil threshold for binarization
puprange = [8 45]; %set

%user input to select center and right points
sprintf('Please select pupil center and top, eyeball top and right points, darkest part of eyeball')
h1 = figure('units','normalized','outerposition',[0 0 1 1])
imshow(data(:,:,10)) %select frame as example for selecting points for rad/centroid calculation
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
toc

startT = max(tsampBar(1),cameraTTL(1))
endT = min(tsampBar(end),cameraTTL(end))

%interpolate because ephys/running not acquired at same rate as camera
%frames
dt = 0.5;
pts = startT:dt:endT;
vInterp = interp1(tsampBar,vsmoothBar,pts);
%rInterp = interp1(cameraTTL(1:length(rad')),rad,pts); %if cameraTTL is
%longer than radius...in case triggering is not exact
rInterp = interp1(cameraTTL,rad(1:length(cameraTTL)),pts); %if rad is longer
figure
plot(vInterp,'g'); hold on; plot(rInterp)
legend('velocity','radius')


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
vidObj = VideoWriter('12_08_16_drift_highcontrast_pre_eye'); %save video with tracking as avi
open(vidObj);
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
%     open(vid); writeVideo(vid,mov); close(vid)
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end
close(vidObj);

% median filter for radius
r = medfilt1(rad,7);
figure;plot(r)

