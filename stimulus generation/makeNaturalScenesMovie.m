%% makeNaturalScenesMovie
%%%generate a movie of natural scenes using familiar and unfamiliar images
%%%for the enrichment mice (PRLP 060419)
clear all
warning off

%directory of familiar and unfamiliar stimuli
unfamdir = 'C:\src\images\naturalScenes\unfamiliar';
famdir = 'C:\src\images\naturalScenes\familiar';
prevmov = 'C:\Users\nlab\Desktop\naturalImages.avi'; % if you want to make a preview movie
nfiles = 13; %if you want only a predetermined number of images selected from the folders (each)
nreps = 12; %number of repititions of each image
previewStim = 0; %if you want to preview the stimuli

%parameters of the stimulus movie
scaledown = 8; %how much to scale the movie down from 1280*720
xsz = 1920/scaledown;ysz = xsz*(9/16); %size of the movie in pixels
framerate = 20; %display framerate (usually 60)
duration = 1; %stimulus duration
isi = 1; %interstimulus interval
blank = 1; %if you want to add an occasional blank image
randomOrder = 1; %if you want to randomize the presentation order (usually yes)

familiar=[];
imfiles = {}; %empty array for image files
totnum=0; %total number of images
dirs = {unfamdir,famdir};

for d = 1:length(dirs)
    cd(dirs{d})
    files = dir('*.png');
    if exist('nfiles','var')
        files = files(1:nfiles);
    end
    num = length(files);
    totnum = totnum + num;
end

sprintf('total duration = %d min',(duration+isi)*(totnum+blank)*nreps/60)

allims = zeros(xsz,ysz,totnum+blank,'uint8')+128;
allfiles = {};
fileidx = [];
%get the images
cnt=1;
for d = 1:length(dirs)
    cd(dirs{d})
    files = dir('*.png');
    if exist('nfiles','var')
        files = files(1:nfiles);
    end
    num = length(files);
    for i=1:num
        im=imread(files(i).name);
        allfiles{cnt} = files(i).name;
        fileidx(cnt) = cnt;
        if size(im,3)>1
            im = squeeze(im(:,:,1));
        end
        [y,x] = size(im);
        im = imresize(im,ysz/y);
        [y,x] = size(im);
        allims(round((xsz-x)/2):x+round((xsz-x)/2)-1,:,cnt) = im.';
        cnt=cnt+1;
    end
    familiar = [familiar ones(1,num)+(d-2)]; %array for familiar/unfamiliar trial IDs
end

if blank
    familiar = [familiar 2];
    fileidx = [fileidx 0];
end

if previewStim
    figure;
    for i=1:totnum+blank
        imshow(squeeze(allims(:,:,end)))
        drawnow
        pause(0.5)
    end
end

%repeat and randomize (finish this)
allims = repmat(allims,[1,1,nreps]); familiar = repmat(familiar,[1,nreps]); fileidx = repmat(fileidx,[1,nreps]);
trial = size(allims,3);
if randomOrder
    order = randperm(trial);
    allims = allims(:,:,order);familiar = familiar(order);fileidx = fileidx(order);
end

moviedata = zeros(xsz,ysz,framerate*(duration+isi)*trial,'uint8')+128;
for tr = 1:trial
    tr
    for t = 1:duration*framerate;
        frame = squeeze(allims(:,:,tr));
        moviedata(:,:,(tr-1)*duration*framerate+tr*isi*framerate+t) = frame;
    end
end

save C:\src\movies\naturalImagesEnrichment8mag648s.mat moviedata allfiles framerate duration isi blank randomOrder nreps familiar fileidx -v7.3

if previewStim
    figure
    for i = 1:length(moviedata)/10
        imshow(imresize(moviedata(:,:,i),0.2));
        drawnow
    end
end

%% make a sample avi
mov = moviedata(:,:,1:size(moviedata,3)/40);
mov = rot90(mov,3);
mov = mat2im(mov,gray);
mov = immovie(permute(mov,[1 2 4 3]));
vid = VideoWriter(prevmov);
vid.FrameRate=framerate;
open(vid);
disp('writing movie')
writeVideo(vid,mov);
close(vid)