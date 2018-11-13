duration = 20;
drift_duration = 2;
interval=0.5;
framerate=60; period = 0.5;
wavelength = 50;
nframes = duration * framerate;
mov = zeros(256,256,nframes)+128;


[x y ] = meshgrid(1:256,1:256);
for orient = 1:8;
    theta = (orient-1)*pi/4;
    rotx = x*cos(theta) + y*sin(theta);
    roty = x*sin(theta)-y*cos(theta);
    for t = 1:drift_duration*framerate;
        mov(:,:,(orient-1)*(drift_duration+interval)*framerate +t) = 255*(cos(2*pi*rotx/wavelength + 2*pi*(t/framerate)/period)>0);
    end
end

moviedata = uint8(mov);

save fullcontrastGratings -v7.3

    

%%
for i = 1:size(moviedata,3);
    imshow(moviedata(:,:,i));
    getframe(gcf);
end