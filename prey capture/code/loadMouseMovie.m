function vid = loadMouseMovie(fn,binsize, startframe, endframe)

if ismac
    fn(fn=='\')='/';
end

matfile = [fn(1:end-3) 'mat'];

try exist(matfile,'file');
    load(matfile,'vid','bin');
    size(vid)
    if bin~=binsize
        error
    end
    
catch
    readerObj = VideoReader(fn);
    display('reading video ...')
    tic
    vid = read(readerObj,[startframe endframe]);
    toc
   % close(readerObj);
    vid = flipdim(imresize(squeeze(vid(:,:,1,binsize:binsize:end)),1/binsize),1);
    bin = binsize;
    save(matfile,'vid','bin','-v7.3');
end