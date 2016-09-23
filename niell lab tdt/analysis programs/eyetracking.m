clear all
close all

dir = 'D:\Angie_analysis\DOI_experiments\01_27_16';
name = '01_27_16_bars_eye'; %data file
thresh = 0.80; %pupil threshold for binarization
puprange = [15 55]; %set

%%closed loop parameters
% pupercent = 0.15; %set range pupil radius window
% pupchange = 0.25; %acceptable percent change in radius per framerange
% framerange = 1; %number of frames to smooth over

load(name);
data = squeeze(data);
warning off;

%user input to select center and right points
sprintf('Please select pupil center and top, eyeball top and right points, darkest part of eyeball')
h1 = figure('units','normalized','outerposition',[0 0 1 1])
imshow(data(:,:,100))
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
for i = 1:size(data,3)
    binmax(i) = mean(mean(mean(ddata(binmaxy-3:binmaxy+3,binmaxx-3:binmaxx+3,i))));
end
for i = 1:size(ddata,3)
    bindata(:,:,i) = (ddata(yc-vert:yc+vert,xc-horiz:xc+horiz,i)/binmax(i) > thresh);
end
figure
imshow(bindata(:,:,1400))

%convert from uint8 into doubles and threshold, then binarize

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
    %%closed loop execution
%     if n>framerange && (isnan(rad(n-1)) | isnan(rad(n))) %if it's a nan or preceeded by all nans don't change puprange
%         puprange = puprange;
%     elseif n>framerange && (abs(1 - rad(n)/nanmean(rad(n-framerange:n-1))) > pupchange) %if % change is bigger than specified don't change puprange
%         puprange = puprange;
%     elseif n>framerange && (rad(n)>nanmean(rad(n-framerange:n-1))) %if radius goes up, shift range up
%         puprange = puprange + round(rad(n) - nanmean(rad(n-1)));
%     elseif n>framerange && (rad(n)<nanmean(rad(n-framerange:n-1))) %if radius goes down, shift range down
%         puprange = puprange - round(nanmean(rad(n-1)) - rad(n));
%     else
%         puprange = puprange;
%     end
end
toc

%plot x and y position and radius across experiment
h2 = figure
hold on
plot(0.1:0.1:size(data,3)/10,rad,'b-')
plot(0.1:0.1:size(data,3)/10,centroid(:,1),'.g')
plot(0.1:0.1:size(data,3)/10,centroid(:,2),'.r')
hold off
legend('radius','x pos','ypos')
ylim([5 55])

% %
h3 = figure
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

end


save(fullfile(dir,name),'centroid','rad','h2','-append');


