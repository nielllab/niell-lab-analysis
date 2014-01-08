%%% adds baseline and shortens movie to 10 frames
for i = 1:300:9000;
    moviedata(:,:,i:i+74)=128;
    moviedata(:,:,i+85:i+224)=128;
    moviedata(:,:,i+235:i+299)=128;
end

%%% adds baseline and shortens movie to 1 frames
for i = 1:300:9000;
    moviedata(:,:,i:i+74)=128;
    moviedata(:,:,i+76:i+224)=128;
    moviedata(:,:,i+226:i+299)=128;
end

%%% adds baseline and make static image
for i = 1:300:9000;
    moviedata(:,:,i:i+74)=128;
    for f = 76:149;
        moviedata(:,:,i+f)=moviedata(:,:,i+75);
    end
    moviedata(:,:,i+150:i+224)=128;
    for f = 226:299;
        moviedata(:,:,i+f)=moviedata(:,:,i+225);
    end
    
end

figure


for i = 1:600
    imshow(squeeze(moviedata(:,:,i)));
    getframe(gcf);
end