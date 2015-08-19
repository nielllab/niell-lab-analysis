clear all
duration = 2;
framerate = 60;
targDur = 0.016;
maskDur = 0.05;
soa = [0 16 33 66 100 150 1000]/1000;

sigma = 0.25;
binarize = 0;

targDur = round(targDur*framerate); maskDur = round(maskDur*framerate); soa = round(soa*framerate);

targtheta = [0 pi/4 pi/2 3*pi/4];
targSF = [0.04 0.16];
masktheta = [0 pi/4 pi/2 3*pi/4];
maskSF = [0.04 0.16];
totalduration = 600;
nrep = totalduration/duration;
nx=1; ny=1;

   
xsz = 128;
ysz = xsz*9/16;
dist = 25;
width = 50;
widthdeg = 2*atand(0.5*width/dist)
degperpix = widthdeg/xsz

xpos = linspace(1,xsz-1,nx+1);
xposrange = xpos(1:end-1);
xposrange = 0.25*xsz;
ypos = linspace(1,ysz-1,ny+1);
yposrange = ypos(1:end-1);
yposrange = 0.05*ysz;

blockwidth = floor(ysz/ny);


[x y] =meshgrid(1:blockwidth,1:blockwidth);
xgrid=(x-mean(x(:)))/max(x(:));; ygrid=(y-mean(y(:)))/max(y(:));
x = x*degperpix; y=y*degperpix;
gaussian = sqrt(exp(-0.5*(xgrid.^2 +ygrid.^2)/sigma^2));
if binarize
    gaussian(gaussian<exp(-1))=0;
end
r= 2*sqrt(xgrid.^2 + ygrid.^2);
mask = ones(size(xgrid));
%mask = mask.*cos((2*r-1)*pi/2);  mask(r>1)=0; mask(r<0.5)=1;
mask = mask.*cos(r*pi/2);  mask(r>=1)=0;
gaussian = gaussian.*mask;


moviedata = zeros(xsz,ysz,totalduration*framerate,'uint8')+128;
for r = 1:nrep
    sf(r,:) = [listPick(targSF) listPick(maskSF)];
    theta(r,:) = [listPick(targtheta) listPick(masktheta)];
    lag(r) = listPick(soa);
    xpos(r,:) =xposrange(1); ypos(r,:) =  yposrange(1);
    ph = (x*cos(theta(r,1)) + y*sin(theta(r,1)))*2*pi*sf(r,1) + rand*2*pi;
 
    for t = 1:targDur;
        frame = uint8(0.5*255*(cos(ph).*gaussian+1));
        moviedata(xpos(r):xpos(r)+blockwidth-1, ypos(r):ypos(r)+blockwidth-1,(r-1)*duration*framerate +t) = frame;
    end

      ph = (x*cos(theta(r,2)) + y*sin(theta(r,2)))*2*pi*sf(r,2) + rand*2*pi;

    for t = 1:maskDur;
        frame = uint8(0.5*255*(cos(ph).*gaussian+1));
        moviedata(xpos(r):xpos(r)+blockwidth-1, ypos(r):ypos(r)+blockwidth-1,(r-1)*duration*framerate +targDur +lag(r)+ t) = frame;
    end
  
end
if binarize
moviedata(moviedata>128)=255;
moviedata(moviedata<128)=0;
end



moviedata = moviedata(1:xsz,1:ysz,:);

figure
for i = 1:16;
    subplot(4,4,i);
    imagesc(moviedata(:,:,i),[0 255]); colormap gray;
end

figure
for i = 1:length(moviedata)/40
imshow(moviedata(:,:,i));
drawnow
end
save mask2sf4theta moviedata xpos ypos maskDur targDur sf lag theta framerate duration  nx ny sigma
