[f p] = uiputfile('*.avi')
movdata = permute(moviedata,[2 1 3]);
movImg = mat2im(movdata(:,:,1:1:round(end/10)),gray,[0 255]);
mov = immovie(permute(movImg,[1 2 4 3]));
vid = VideoWriter(fullfile(p,f));
% mov = immovie(permute(shiftmov,[1 2 4 3]));
% vid = VideoWriter('bilateralS1.avi');
vid.FrameRate=60;
open(vid);
writeVideo(vid,mov);
close(vid)