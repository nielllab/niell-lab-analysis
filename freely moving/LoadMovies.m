%% Load Full Movies

clear; close all;
% Choose avi file
[file,folder] = uigetfile({'*.avi';'*.mov';},'choose First Video');
% Read in Eye Tracking Video
if file ~=0
    TempVid = VideoReader(fullfile(folder,file));
    frame=1; k=1;
    while hasFrame(TempVid)
        Vid1(:,:,frame) = rgb2gray(readFrame(TempVid));
        if mod(frame,100)==0
            fprintf('frame = %d\n',frame)
        end
        frame=frame+1;
        
    end
    % Load the .csv file generated by DLC
    fileListDLC1 = dir([folder file(1:end-4) '*.csv']);
    DLCPoints = csvread(fullfile(folder,fileListDLC1.name),3,0);
    vidsize=1:size(Vid1,3); %1:1000;
    Pointsx = DLCPoints(vidsize,[2:3:end]); Pointsy = DLCPoints(vidsize,[3:3:end]);
    LH =  DLCPoints(vidsize,[4:3:end]);
end
vstart = 0; vstop = size(Vid1,3);
fprintf('Done Reading in Video \n');

%% Load only specific interval of frames
clear; close all;
vstart = 3300;  vstop = 3600; % Which frames to grab

[file,folder] = uigetfile({'*.avi';'*.mov';},'choose First Video');
if file ~=0
    % Video
    TempVid = VideoReader(fullfile(folder,file)); % Create video reader
    tvideo = (read(TempVid,[vstart vstop])); % Read in froms from vstart to vstop
    for fm=1:size(tvideo,4)
        Vid1(:,:,fm) = rgb2gray(tvideo(:,:,:,fm)); % Convert to gray scale.
    end
    % DLC Points
    fileListDLC1 = dir([folder file(1:end-4) '*.csv']);
    DLCPoints = csvread(fullfile(folder,fileListDLC1.name),3,0);
    vidsize=[vstart:vstop]; %
    Pointsx = DLCPoints(vidsize,[2:3:end]); Pointsy = DLCPoints(vidsize,[3:3:end]);
    LH =  DLCPoints(vidsize,[4:3:end]);
    
end
fprintf('Done\n')
%% Check for bad fits

j=1;
bdfit=[];
for i=1:size(LH,1)
    if isempty(find(LH(i,:)<.9))~=1  % Checking if Likelihood is < .9
        bdfit(j,1) = i; j=j+1;
    end
end
Pointsx(bdfit,:)=NaN; Pointsy(bdfit,:)=NaN; 

%% Plot video with DLC points

figure;
axis tight manual
ax = gca;
% ax.NextPlot = 'replaceChildren';
% axis([1 size(Vid,1) 1 size(EyeVid,2)] )
for fm = 1:size(Vid1,3)
    % Show Image Frame
    imagesc(Vid1(:,:,fm)); colormap gray; axis equal off; hold on;
    % Plot DLC points
    scatter(Pointsx(fm,:),Pointsy(fm,:),100,'.r'); % Uncomment if you wish plot points
    % If eye cam uncomment to plot ellipse fit. 
%     if ~any((Pointsx(fm,:)<10 | Pointsy(fm,:)<10)) % check points
%         e_t = fit_ellipse(Pointsx(fm,:),Pointsy(fm,:),ax);
%     end
    
    hold off;
    title(sprintf('Frame = %d',fm+vstart));
    drawnow limitrate
end



%% If Eye video Calculate Ellipse fit and angles

[theta,phi,EllipseParams,~] = EyeCameraCalc1(Vid1,Pointsx,Pointsy);







