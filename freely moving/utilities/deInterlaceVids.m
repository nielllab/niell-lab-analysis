function deInterlaceVids(movieFilename);
%%% De-interlace video to double framerate and remove interlacing artifact
%%% Separates out every other line and places them into separate frames
%%% Down-size along the other axis, to keep pixels square
%%% Drops 2x in spatial resolution to gain 2x in temporal resolution
%%% (but not really a loss in spatial since the original resolution had artifact
%%% Ideal for NTSC/PAL cameras (e.g. head/eye cams)
%%%
%%% input = .avi movie filename (including path); if not provided, will ask
%%%  re-saves movie as moviename_DeInter.avi
%%%
%%% cmn 2019


%%% read in original movie file
if ~exist('movieFilename','var');
[f p] = uigetfile('*.avi','world video');
movieFilename = fullfile(p,f);
end

TempVidT = VideoReader(movieFilename);
frame=1; k=1;

while hasFrame(TempVidT)
    worldLaced(:,:,:,frame) = (readFrame(TempVidT));
    
    if mod(frame,100)==0
        fprintf('frame = %d\n',frame)
    end
    frame=frame+1;
    
end

world = zeros(size(worldLaced,1)/2, size(worldLaced,2)/2, 3, size(worldLaced,4)*2,'uint8');
sz = [size(world,1) size(world,2)];
world(:,:,:,1:2:end) = imresize(worldLaced(1:2:end,:,:,:),sz);
world(:,:,:,2:2:end) = imresize(worldLaced(2:2:end,:,:,:),sz);

fnew = [movieFilename(1:end-4) '_DeInter.avi'];
movObj = VideoWriter(fnew);
movObj.FrameRate = 60;
open(movObj);
writeVideo(movObj,immovie(world));
close(movObj);

