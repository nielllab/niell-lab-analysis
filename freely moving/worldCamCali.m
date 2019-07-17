clear; close all;

% Read in Eye Video
[f,p] = uigetfile('.avi');
% fileListE = dir([f '\Eye\' '*eye*.avi']);
frame=1;
% Read in Video
if f ~=0
TempVid = VideoReader(fullfile(p,f));

    for frame=1:100


    Calib(:,:,frame) = rgb2gray(readFrame(TempVid));
    end
end
%%
[centers,radii,~] = imfindcircles(Calib(:,:,1),[50 150],'ObjectPolarity','bright','Sensitivity',0.97);

px_deg = radii/(atand(5.08/30));
%%
f = uigetdir();
fileListE = dir([f '\World\' '*world*.avi']);

% Read in Eye Tracking Video
if f ~=0
    TempVid = VideoReader(fullfile([f '\World\'],fileListE.name));
    frame=1; 
    while hasFrame(TempVid)
        WorldVid(:,:,frame) = rgb2gray(readFrame(TempVid));
        if mod(frame,100)==0
            fprintf('frame = %d\n',frame)
        end
        frame=frame+1;
        
    end
end
fprintf('Done Reading in World Video \n');
%

shiftxW = newtheta*px_deg;


parfor i =1: size(shiftxW,1) - 1 
    EyeVid2(:,:,i) = imtranslate(EyeVid2(:,:,i),[shiftxW(i,1) shiftyW(i,1)]);
end