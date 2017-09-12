%Tristan Ursell
%April 2017
%Conversion from all MATLAB readable formats to TIFF stack
%
 %video2tiff(vid_file)
 %stack_out=video2tiff(vid_file,'write');
%vid_file='live_1_LT_1hr_postCNO.mov';


%pathname = '/Users/jennifer/Dropbox/mucsimol_V1/Trial videos for Tristans tracking/Iphone5/test_green_ears/'
%cd '/Users/jennifer/Dropbox/mucsimol_V1/Trial videos for Tristans tracking/Iphone5/test_green_ears/'
% 

%Goal 1: generate way to input pathname and file


%vid_file = 'E:\Dropbox\mucsimol_V1\Trial videos for Tristans tracking\Iphone5\live_1_LT_1hr_postCNO\live_1_LT_1hr_postCNO.mov';
% files=dir(fullfile(,'*.mov'));
% vid_file=files;

function varargout=video2tiff(vid_file,varargin)
writeq=0;

if nargin>1
    if strcmp(varargin{1},'write')
        writeq=1;
        
        [path1,name1,~]=fileparts(vid_file);
        outname=fullfile(path1,[name1 '.tiff']);
    else
        error([varargin{1} ' is not a recognized input.'])
    end
end

%load video file object
vid1=VideoReader(vid_file);
vid1.CurrentTime = 0.0;

if ~writeq
    Nframe=vid1.Duration*vid1.FrameRate;
    frameDepth=vid1.BitsPerPixel/8;
    stack_out=zeros(vid1.Height,vid1.Width,frameDepth,ceil(Nframe),'uint8');
end

q=0;
while hasFrame(vid1)
    q=q+1;

    if writeq
        image_save(readFrame(vid1),outname)
    else
        stack_out(:,:,:,q)=readFrame(vid1);
    end
end

if nargout>0
    if writeq
        varargout{1}=outname;
    else
        varargout{1}=stack_out(:,:,:,1:q);
    end
end

%{
currAxes = axes;
while hasFrame(v1)
    vidFrame = readFrame(v1);
    image(vidFrame, 'Parent', currAxes);
    currAxes.Visible = 'off';
    pause(1/v1.FrameRate);
end
%}



%*************************************************************************
%*************************************************************************
%Image Stack Saver
%Tristan Ursell, February 2012
%
% image_save(Im1,basename)
% image_save(Im1,basename,fmax)
%
%For use with appended image stacks, e.g. most commonly stacked TIF files.
%Some operating systems throw an error indicating that the file is
%inaccessible at the time of write, which leads to a code fault and can be
%very frustrating.  This simple script fixes that issues, mainly on the
%Windows operating system.
%
% Im1 = the current image (matrix) you want to add to the stack
% 
% basename = is a string that specifices the file name, and potentially the
% path of the stack to save to.  Best practice is the put '.tif' at the end
% of the file name.  
%
% fmax = maximum number of allowed failures, i.e. if the script attempts to
% write to the file 'fmax' times and fails, it will give up.  This prevents
% entering an infinite loop.  The default value is fmax = 10.
%
% Any of the 'imwrite' parameters can be modified below on line 38.
%

function image_save(Im1,basename,varargin)

if ~isempty(varargin)
    fmax=varargin{1};
else
    fmax=10;
end

er1=0;
f=0;
while and(er1==0,f<fmax)
    try
        imwrite(Im1,basename,'writemode','append','compression','none');
        er1=1;
    catch
        f=f+1;
        pause(0.1*rand)
    end
end

if f==fmax
    error(['File: ' basename ', remained inaccessible after ' num2str(fmax) ' attempts.'])
end


