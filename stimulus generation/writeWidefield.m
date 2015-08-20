close all

word = 'HELLO';
moviedata = zeros(128,128,600);
nc = length(word);
isi =102;
duration = (600-nc*isi)/nc;

for i = 1:length(word);
    
    figure
    imshow(zeros(256,256));
    
    text(64,128,word(i),'FontSize',96,'Color','w');
    drawnow
    im= getframe(gca);
    im = im.cdata;
    
    im = squeeze(im(1:256,1:256,1));
    im(im<255)=0;
    mass = mean(im,1);
    x = 1:length(mass');
    centroid = (min(find(mass>0)) +max(find(mass>0)))/2;
    image{i} = im(65:192,centroid-63:centroid+64);
    figure
    imshow(image{i})
    for f = (i-1)*(duration+isi)+(1:duration)
        moviedata(:,:,f) = image{i}*255;
    end
   
end

moviedata=uint8(moviedata);
moviedata(moviedata==0)=128;
moviedata(moviedata==255)=0;


for f = 1:size(moviedata,3);
    moviedata(:,:,f) = flipud(moviedata(:,:,f))';
end


moviedata = imresize(moviedata,[1.5*size(moviedata,1) size(moviedata,2)]);
filt = ones(1,8)/8;
moviedata=imfilter(moviedata,filt,'replicate');
moviedata(moviedata<128)=0;

moviedata(end+1:end+36,:,:)=128;
filt = fspecial('gaussian',25,4);

moviedata = imfilter(moviedata,filt,'replicate');
moviedata(moviedata>32)=128;


filt = fspecial('gaussian',25,2);
moviedata = imfilter(moviedata,filt,'replicate');
moviedata = imresize(moviedata,0.5);
figure

    for f = 1:600;
        imshow(moviedata(:,:,f));
        drawnow
    end
    moviedata = repmat(moviedata,[1 1 30]);
    figure
    imshow(moviedata(:,:,1))
    save hellomovie_black_short moviedata;
    