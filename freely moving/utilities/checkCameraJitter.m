function [jitter trace] = checkCameraJitter(movieFilename);
%%% check for camera jitter by doing image registration on a small patch
%%% in a region of the video that shouldn't move
%%%
%%% input = .avi movie filename (including path); if not provided, will ask
%%%  re-saves movie as moviename_DeInter.avi, returned in fnew
%%%
%%% cmn 2020

%%% read in movie file if not passed as argument
if ~exist('movieFilename','var');
    [f p] = uigetfile('*.avi','world video');
    movieFilename = fullfile(p,f);
end

%%% read in movie
vid = VideoReader(movieFilename);
frame=1; k=1;
display('reading')

while hasFrame(vid)
    img(:,:,:,frame) = readFrame(vid);
    
    %%% status update
    if mod(frame,200)==0
        fprintf('frame = %d\n',frame)
    end
    frame=frame+1;   
end

%%% convert to bw double image
bw = double(squeeze(img(:,:,1,:)))/255;

% figure
% imshow(mean(bw(:,:,1:10:end),3))
% 
% figure
% imshow(bw(:,:,100))
% 
% figure
% imagesc(std(bw(:,:,1:10:end),[],3))

close all

%%% crop image
display('please crop image')
figure
[cropImg cropRect] = imcrop(mean(bw(:,:,1:10:end),3));
cropRect = round(cropRect);

imref = mean(bw(cropRect(2) + (1:cropRect(4)), cropRect(1) + (1:cropRect(3)), 1:10:end),3);

figure
imshow(imref)

%%% perform cross-correlation alignment
n = 0;
range = -10:10;
tic
clear xoffset yoffset
for f = 1:10:size(bw,3);  %%% loop over every 10th frame
   f
   n= n+1;
    for x = 1:length(range)  %%% loop over offsets and calculate difference with ref image
        for y= 1:length(range)
            thisImg = bw(cropRect(2) + range(y) + (1:cropRect(4)), cropRect(1) + range(x)+ (1:cropRect(3)), f);
            err(x,y) = sum(sum((imref - thisImg).^2));
        end
    end
    %%% find location of minimum error and convert it back to offset
    minErr = min(err(:));
    [ymin xmin] = ind2sub(size(err),find(err==minErr));
    yoffset(n) = range(ymin); xoffset(n) = range(xmin);
end
toc
    
%%% sometimes fails to align (rails at one end) due to image being
%%% obscured; set these to NaNs
yoffset(abs(yoffset)==range(end))=NaN;
xoffset(abs(xoffset)==range(end))=NaN;

%%% set data to return
trace(:,1) = xoffset;
trace(:,2) = yoffset;
jitter(1) = nanstd(xoffset);
jitter(2) = nanstd(yoffset);

%%% plot
figure
plot(xoffset); hold on; plot(yoffset+0.5);
ylim([range(1) range(end)]);
title(sprintf('x = %0.2f pix; y = %0.2f pix', jitter(1),jitter(2)));
    
            
            
    

