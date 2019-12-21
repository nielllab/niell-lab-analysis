function fnew = deInterlaceVids(movieFilename);
%%% De-interlace video to double framerate and remove interlacing artifact
%%% Separates out every other line and places them into separate frames
%%% Down-size along the other axis, to keep pixels square
%%% Drops 2x in spatial resolution to gain 2x in temporal resolution
%%% (but not really a loss in spatial since the original resolution had artifact
%%% Ideal for NTSC/PAL cameras (e.g. head/eye cams)
%%%
%%% input = .avi movie filename (including path); if not provided, will ask
%%%  re-saves movie as moviename_DeInter.avi, returned in fnew
%%%
%%% cmn 2019

doubleSize=1;

%%% read in original movie file if not passed as argument
if ~exist('movieFilename','var');
    [f p] = uigetfile('*.avi','world video');
    movieFilename = fullfile(p,f);
end

%%% read in movie
TempVidT = VideoReader(movieFilename);
frame=1; k=1;
display('reading')

while hasFrame(TempVidT)
    worldLaced(:,:,:,frame) = (readFrame(TempVidT));
    
    %%% status update
    if mod(frame,500)==0
        fprintf('frame = %d\n',frame)
    end
    frame=frame+1;   
end

%%% separate out odd/even lines, and put them as consecutive frames
% world = zeros(size(worldLaced,1)/2, size(worldLaced,2)/2, 3, size(worldLaced,4)*2,'uint8');
% sz = [size(world,1) size(world,2)];
% world(:,:,:,1:2:end) = imresize(worldLaced(1:2:end,:,:,:),sz);
% world(:,:,:,2:2:end) = imresize(worldLaced(2:2:end,:,:,:),sz);

if doubleSize
world = zeros(size(worldLaced,1), size(worldLaced,2), 3, size(worldLaced,4)*2,'uint8');
else
    world = zeros(size(worldLaced,1)/2, size(worldLaced,2)/2, 3, size(worldLaced,4)*2,'uint8');
end
sz = [size(world,1) size(world,2)];
world(:,:,:,1:2:end) = imresize(worldLaced(1:2:end,:,:,:),sz);
world(:,:,:,2:2:end) = imresize(worldLaced(2:2:end,:,:,:),sz);


%%% plot as troubleshooting
% figure
% imshow(worldLaced(:,:,:,1));
% figure
% imshow(world(:,:,:,1));

%%% re-save as new movie
display('saving')
if doubleSize
    fnew = [movieFilename(1:end-4) '_DeInter2_100.avi'];
else
fnew = [movieFilename(1:end-4) '_DeInter100.avi'];
end

movObj = VideoWriter(fnew);
movObj.FrameRate = 60;
movObj.Quality=100;
open(movObj);
writeVideo(movObj,immovie(world));
close(movObj);

