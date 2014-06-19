duration =10;
radius = 30;
framerate = 60;
szx = 256; szy = 256;
[x y] = meshgrid(1:szy, 1:szx);
moviedata = zeros(szx,szy,duration*framerate);

for f = 1:duration*framerate;
      if f>1 & f<30
         r = radius*(f)/30;
        x0 = szx/4;
        y0= szy/2;
        inside = ((x-x0).^2 + (y-y0).^2) < r^2;
        moviedata(:,:,f) = inside';
      elseif f>150 & f<180;
        r = radius*(f-150)/30;
        x0 = szx/4;
        y0= szy/2;
        inside = ((x-x0).^2 + (y-y0).^2) < r^2;
        moviedata(:,:,f) = -inside';
    elseif f>300 & f<330
        r = radius*(f-300)/30;
        x0 = 3*szx/4;
        y0= szy/2;
        inside = ((x-x0).^2 + (y-y0).^2) < r^2;
        moviedata(:,:,f) = inside';
    elseif f>450 & f<480
                r = radius*(f-450)/30;
        x0 = 3*szx/4;
        y0= szy/2;
        inside = ((x-x0).^2 + (y-y0).^2) < r^2;
        moviedata(:,:,f) = -inside';

    end;
end

moviedata = uint8(0.5* (1-moviedata) *255);
moviedata = moviedata(:,49:208,:);

figure
for i = 1:size(moviedata,3);
i
    imshow(moviedata(:,:,i))
    getframe(gcf);
end





