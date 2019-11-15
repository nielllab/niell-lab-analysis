close all; clear all
pathname = {'./'};

fileList = [];
for i = 1:length(pathname)
    fileList = [fileList ; dir([pathname{i} '*.avi'])];
end

for f = 1:length(fileList);

 f
 
    movName = fileList(f).name;
    %thispath = fileList(f).folder;
    thispath = pathname{1};
    topMov = VideoReader(movName);
    nframe = floor(topMov.duration*topMov.framerate)
    mov = zeros(topMov.Height,topMov.Width,3,nframe,'uint8');   
   
    display('reading movie');
    for i = 1:nframe
        if round(i/100) == i/100
            sprintf('file %d : done %d / %d',f,i,nframe)
        end
        mov(:,:,:,i) = topMov.readFrame;
    end

    display('writing movie');
    movieF = strrep(movName,'.avi','_mat.avi');
    
    outFileObj = VideoWriter(fullfile(thispath,movieF));
    open(outFileObj);
    writeVideo(outFileObj,mov);
    close(outFileObj)
    


end
