function [ stack ] = MovToTiff( vid_file )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%cd('E:\Dropbox\mucsimol_V1\Trial videos for Tristans tracking\Iphone5\live_1_LT_1hr_postCNO\')

[path1,name1,~]=fileparts(vid_file);
workingDir = path1;
mkdir(name1); % path1,[name1 '.tiff']

liveVid = VideoReader(vid_file);

ii = 1;

while hasFrame(liveVid)
   img = readFrame(liveVid);
   filename = [sprintf('%06d',ii) '.tif'];
   fullname = fullfile(workingDir,name1,filename);
   imwrite(img,fullname)    % Write out to a tif file (img1.jpg, img2.jpg, etc.)
   ii = ii+1;
end

end

