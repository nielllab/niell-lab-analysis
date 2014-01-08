
[x y] = meshgrid(1:128,1:128);
grating = cos(2*pi*x/25);

framerate=60;
tf = 1;
m= zeros(128,128,6000);
for f = 1:6000
    if mod(f-1,300)>=150
        m(:,:,f) =grating*sin(2*pi*(f/framerate)*tf);
    end
end
moviedata = uint8(128*m+128);
figure

for f = 1:600;
    imshow(moviedata(:,:,f));
    getframe(gcf);
    f
end

clear m f x y framerate