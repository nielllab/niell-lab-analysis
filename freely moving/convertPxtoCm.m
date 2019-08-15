clear; close all;
% Choose avi file 
[f,p] = uigetfile({'*.avi';},'choose Left eye EyeFilt');
% Read in Eye Tracking Video
if f ~=0
    TempVid = VideoReader(fullfile(p,f));
    frame=1; k=1;
    while hasFrame(TempVid) 
        EyeVidL(:,:,frame) = rgb2gray(readFrame(TempVid));
        if mod(frame,100)==0
            fprintf('frame = %d\n',frame)
        end
        frame=frame+1;
        
    end
    
end

h1 = figure('units','normalized','outerposition',[0 0 1 1])
imshow(EyeVidL(:,:,24))
[cent] = ginput(3);
close(h1);
start=cent(1,:)
fin = cent(2,:)
fin2=cent(3,:)

pt1=cent(1,1)
pt2=cent(2,1)
pt3=cent(3,1)

cm1=pt2-pt1
cm2=pt3-pt2
cm3=(pt3-pt1)/2


%%
in1=start-fin
in2=fin-fin2
in22=start-fin2

(in1(:,2)+in2(:,2))/2




