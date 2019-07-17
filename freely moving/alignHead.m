function data = alignHead(fname, nPts,showMovies)
%%% reads in csv data from deepLabCut with head and cricket positions
%%%
%%% computes head position, even in presence of noisy/absent points, by
%%% aligning and fitting to relatively locations of points
%%%
%%% some aspects assume 8 points in standard configuration (e.g. first =
%%% noise, last = middle of head)
%%%
%%% also loads cricket position and computes relative range, azimuth
%%%
%%% input:
%%% fname = csv filename
%%% showMovies = show movies or not?
%%% returns :
% % % data.mouse_xy = mouse x,y coords (pix)
% % % data.mouseSp = smoothed mouse speed (pix/frame)
% % % data.theta = thAll;  head theta in world coordinates, radians
% % % data.dTheta = dTheta; head angular velocity, rad/frame
% % % data.crick_xy = crick;  cricket x,y position (pix)
% % % data.crickSp = crickSp  smoothed cricket speed (pix/frame)
% % % data.range = range;   distance to cricket from mouse (pix)
% % % data.az = az;    angle of cricket relative to head direction


if ~exist('fname','var')
    fname =  'top_cricket1_062519_3DeepCut_resnet50_TopVidJun25shuffle1_900000_numeric.csv';
end

if ~exist('showMovies','var')
    showMovies = 1; %%% change this to 0 after debugging
end

%%% load data
%%% note : removed top 3 rows from original csv file since they are non-numeric
%data = dlmread(fname);
data = dlmread(fname,',',3,1)

%%% likelihood threshold for including pts
p_thresh = 0.999;

nframes = 100

%%% head pts
%%% pts(npts, x/y, t);
clear p pts
npts = 8;
for i = 1:npts
    p(i,:) = data(:,i*3+1);
    pts(i,:,:) = data(:,(i-1)*3 + (2:3))';
    ptsRaw(i,:,:) = pts(i,:,:); %%% keep a clean version before inserting NaNs
    pts(i,:,p(i,:)<p_thresh)=NaN;
end

%%% find times when all pts are good
good = p>p_thresh; %%% select points over threshold
useN = sum(good,1);
figure
plot(useN); xlabel('frame'); ylabel('#of good points'); ylim([0 npts]);
use = useN==npts; %%% for now, only use times when all points are good

badFraction = 1-mean(good,2);
figure
bar(badFraction); ylabel('fraction bad timepoints'); xlabel('point #')


%%% draw all points on tracks
% figure
% hold on
% for i = 1:npts
%     plot(squeeze(pts(i,1,:)),squeeze(pts(i,2,:)));
% end

%%% calculate centroid at each timpoint (will be NaN for ones with a bad
%%% timepoint
centroid = squeeze(mean(pts,1));
% figure
% plot(centroid(1,:),centroid(2,:));
for i = 1:npts
    centered(i,:,:) = squeeze(pts(i,:,:))-centroid;
end


%%% draw centered points
figure
title('only 8 good points')
if showMovies
    for t= 1:nframes;
        if ~isnan(centroid(1,t))
            %plot(squeeze(centered(:,1,t)),squeeze(centered(:,2,t)),'o')
            hold off
            drawHead(centered(:,:,t)); axis square; axis([-70 70 -70 70])
            drawnow
        end
    end
end

%%% choose a reference image
%%% for now, this is just the first good one
%%% could do something like bootstrapping, choose 10 random ones and averge
%%% the results
refnum = min(find(use));
ref = centered(:,:,refnum);
% figure
% drawHead(ref)

%%% rotate all data to align to this one
display('aligning good data')
tic
aligned = zeros(size(centered));
for t = 1:size(centroid,2)
    
    if use(t)  %%% only align ones that are all good points
        c = centered(:,:,t);
        
        %%% loop through range of thetas, rotate image by that much, and
        %%% calculate how well it matches the ref
        theta = linspace(0,2*pi,101); theta = theta(1:end-1);
        for i = 1:length(theta)
            c_rot = c*rotmat(theta(i))';  %%% rotation
            d(i) = sum((ref(:) -c_rot(:)).^2);  %%% rms error
        end
        
        %%% find smallest error, and rotate by this amount
        [y ind] = min(d);
        th(t) = theta(ind);
        aligned(:,:,t) = c*rotmat(th(t))';
    else
        th(t) = NaN;
        aligned(:,:,t) = NaN;
    end
end
toc

%%% calculate the mean Head
meanHead = nanmean(aligned,3);

%%% rotate mean head to align to x-axis
longAxis = meanHead([8 1],:); %%% line between middle of head and nose points
longTheta = atan2(diff(longAxis(:,2)), diff(longAxis(:,1)));  %%% angle of line
headRot= rotmat(-longTheta);  %%% rotation matrix to fix this
for i = 1:size(aligned,3);
    aligned(:,:,i) = aligned(:,:,i)*headRot';
end
meanHead = nanmean(aligned,3);

%%% show aligned points
figure
hold on
for i = 1:npts
    plot(squeeze(aligned(i,1,:)),squeeze(aligned(i,2,:)),'o');
end
drawHead(meanHead); axis square; axis equal
title('alignment from times with all good points')

if showMovies
    figure
    for t= 1:nframes;
        if ~isnan(centroid(1,t))
            hold off
        drawHead(meanHead)
        drawHead(aligned(:,:,t)); axis square; axis([-70 70 -70 70])
        drawnow
        end
    end
end



%%% calculate centroid for all points
%%% this uses the fact that the centroid has defined distance from marked
%%% points in the meanHead.
%%% here we calculate the x,y centroid that best batches these distances
meanD = sqrt(meanHead(:,1).^2 + meanHead(:,2).^2);  %%% distances from centroid (0,0) for meanhead

display('getting all centroids')
tic
for t = 1:size(centroid,2);  %%% all timepoints
    c = pts(:,:,t);
    %%% make a meshgrid that covers the x,y position of all head points at
    %%% this time
    [x y] = meshgrid(floor(min(c(:,1))):ceil(max(c(:,1))), floor(min(c(:,2))):ceil(max(c(:,2))));
    
    %%% for each head point calculate how far the pixels are from it,
    %%% then calculate error of how far this is from where it should be,
    %%% then add these up
    err = 0;
    for i = 1:npts
        if ~isnan (c(i,1))
            r = sqrt((x-c(i,1)).^2  + (y-c(i,2)).^2); %%% distance
            err  = err+ (meanD(i) - r).^2;  %%% error
        end
    end
    %%% find minimum, then get x and y values and set as centroid
    [v ind] = min(err(:));
    [i j] = ind2sub(size(err),ind);
    cent(1,t) = x(i,j); cent(2,t) = y(i,j);
end
toc
%%% center all points using calculated centroid
for i = 1:npts
    centered(i,:,:) = squeeze(pts(i,:,:))-cent;
end

figure
plot(cent(1,:),cent(2,:), 'g','Linewidth',2);
hold on
plot(centroid(1,:),centroid(2,:),'k')
legend('all points','only good points')

%%% draw all centered
figure
if showMovies
    for t= 1:nframes;
        %plot(squeeze(centered(:,1,t)),squeeze(centered(:,2,t)),'o')
        hold off
        drawHead(centered(:,:,t)); axis square; axis([-70 70 -70 70])
        drawnow
    end
end


%%% now align all timepoints
%%% works similar to above, but only sum error over good marked points

display('aligning all data - slow')
tic
alignedAll = zeros(size(centered));
clear d
for t = 1:size(centroid,2)
    c = centered(:,:,t);
    
    if useN(t)>=4  %%% need 4 points to fit well
        %%% loop over thetas, rotate points, calculate error
        theta = linspace(0,2*pi,1001); theta = theta(1:end-1);
        clear d
        for i = 1:length(theta)
            c_rot = c*rotmat(theta(i))';
            d(i) = nansum((meanHead(:) -c_rot(:)).^2);  %%% nansum only includes good points
        end
        %%% find minimum and rotate points accordingly
        [y ind] = min(d);
        thAll(t) = theta(ind);
        alignedAll(:,:,t) = c*rotmat(thAll(t))';
    else
        thAll(t) = NaN;
        alignedAll(:,:,t)=NaN;
    end
end
toc

thAll = 2*pi-thAll; %%% because angle of head is actually negative of what we needed to correct it
thAll(thAll>pi) = thAll(thAll>pi)-2*pi;  %%% range = -pi : pi
figure
plot(thAll); title('final theta')
xlabel('frame'); ylabel('theta');

meanHeadAll = nanmean(alignedAll,3);
figure
hold on
for i = 1:npts
    plot(squeeze(alignedAll(i,1,:)),squeeze(alignedAll(i,2,:)),'o');
end
drawHead(meanHeadAll); axis square; axis equal
title('alignment from all timepoints')

if showMovies
    figure
    for t= 1:nframes;
        hold off
        drawHead(meanHeadAll)
        drawHead(alignedAll(:,:,t)); axis square; axis([-70 70 -70 70])
        drawnow
    end
end

%%% cricket pts
crick = data(:,npts*3 +(2:3))';
crick_p = data(:,npts*3 + 4);

figure
plot(crick_p)

crick(:,crick_p<p_thresh) = NaN;

vx = diff(crick(1,:)); vy = diff(crick(2,:));
filt = ones(3,1); filt = filt/sum(filt);
vx = conv(vx,filt,'same'); vy = conv(vy,filt,'same');
crickSp = sqrt(vx.^2 + vy.^2);

figure
subplot(2,3,1);
plot(crickSp); title('cricket speed'); xlabel('frames'); ylabel('pix / sec')

r = crick - cent;
range = sqrt(r(1,:).^2 + r(2,:).^2);

subplot(2,3,2)
plot(range);
title('cricket range'); xlabel('frame'); ylabel('pix')

crickTheta = atan2(r(2,:),r(1,:));
az = thAll - crickTheta;
subplot(2,3,3)
plot(az)
title('azimuth'); xlabel('frame'); ylabel('radians')

dTheta = diff(thAll);
dTheta(dTheta>pi) = dTheta(dTheta>pi)-2*pi;
dTheta(dTheta<-pi) = dTheta(dTheta<-pi) + 2*pi
filt = ones(3,1); filt = filt/sum(filt);
dTheta = conv(dTheta,filt,'same');
subplot(2,3,4)
plot(dTheta)
title('head angular velocity'); ylabel('rad/frame'); xlabel('frame')

vx = diff(cent(1,:)); vy = diff(cent(2,:));
filt = ones(3,1); filt =filt/sum(filt);
vx = conv(vx,filt,'same'); vy = conv(vy,filt,'same');
mouseV = sqrt(vx.^2 + vy.^2);

subplot(2,3,5)
plot(mouseV)
title('mouse speed'); xlabel('frame'); ylabel('pix / sec')

subplot(2,3,6)
hold on
for i = 1:npts
    plot(squeeze(alignedAll(i,1,:)),squeeze(alignedAll(i,2,:)),'.');
end
drawHead(meanHeadAll); axis square; axis equal

clear data
data.mouse_xy = cent;
data.mouseSp = mouseV;
data.theta = thAll;
data.dTheta = dTheta;
data.crick_xy = crick;
data.crickSp = crickSp
data.range = range;
data.az = az;



cricket.pos = crick;


sz = max(pts(:));
if showMovies
    figure
    for t = 1:size(pts,3);
        hold off
        currentHead = meanHeadAll*rotmat(thAll(t))';
        currentHead = currentHead + repmat(cent(:,t)',[8 1]);
        drawHead(currentHead,1);
        plot(squeeze(ptsRaw(:,1,t)),squeeze(ptsRaw(:,2,t)),'b.','MarkerSize',12);
        plot(crick(1,t),crick(1,t),'ro')
        axis([1 sz 1 sz])
        drawnow
    end
end

