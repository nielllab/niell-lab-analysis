clear all
duration = 1;
framerate = 60;
targDur = 0.033;
maskDur = 0.033;
soa = [0 33 66 100]/1000;

sigma = 5; %width of targ / annulus (deg)
r = 15; % radius of annulus (deg)
radius =r;
binarize = 0;

targDur = round(targDur*framerate); maskDur = round(maskDur*framerate); soa = round(soa*framerate);

deltaTheta = [ 0 pi/2];
targSF = [0 0.04 0.16];
maskSF = [0 0.04 0.16];
phases = linspace(0,2*pi,7);
phases= phases(1:end-1);
nx=2; ny=1;

nrep = length(soa)*length(deltaTheta)*length(targSF)*length(maskSF)*length(phases)*nx

totalduration = nrep*duration
   
xsz = 128;
ysz = xsz*9/16;
dist = 25;
width = 50;
widthdeg = 2*atand(0.5*width/dist)
degperpix = widthdeg/xsz

blockwidth = ysz;

xposRange = [1 xsz - blockwidth];
yposRange  = 1;


[x y] =meshgrid(1:blockwidth,1:blockwidth);
xgrid=(x-mean(x(:)));; ygrid=(y-mean(y(:)));
x = xgrid*degperpix; y=ygrid*degperpix;
gaussian = (exp(-0.5*(x.^2 +y.^2)/sigma^2));
if binarize
    gaussian(gaussian<exp(-1))=0;
end
dist = sqrt(x.^2 + y.^2);
mask  = exp((-0.5*(dist-r).^2) / sigma^2);
figure
imagesc(gaussian,[ 0 1]);
figure
imagesc(mask,[0 1]);


figure
imagesc(max(gaussian,mask),[0 1])
moviedata = zeros(xsz,ysz,totalduration*framerate);
for r = 1:nrep
   r
   sf(r,:) = [listPick(targSF) listPick(maskSF)];
    theta(r) = rand(1)*2*pi;
    dOri(r) = listPick(deltaTheta); 
    lag(r) = listPick(soa);
    xpos(r) =listPick(xposRange);
    ypos(r) = yposRange(1);
    ph = (x*cos(theta(r)) + y*sin(theta(r)))*2*pi*sf(r,1) + rand*2*pi;
    xr = xpos(r):xpos(r)+blockwidth-1; yr = ypos(r):ypos(r)+blockwidth-1; tr = (r-1)*duration*framerate;
   if sf(r,1)>0
    for t = 1:targDur;
        frame = cos(ph).*gaussian;
        moviedata(xr, yr,tr +t) = moviedata(xr, yr,tr +t)  + frame;
    end
   end
      ph = (x*cos(theta(r)+dOri(r)) + y*sin(theta(r)+dOri(r)))*2*pi*sf(r,2) + rand*2*pi;

     if sf(r,2)>0
    for t = 1:maskDur;
        frame = cos(ph).*mask;
        moviedata(xr, yr,tr +t + lag(r)) = moviedata(xr, yr,tr +t+lag(r))  + frame;
    end
     end
end
moviedata(abs(moviedata)<(1/128))=0;
moviedata =uint8(floor(moviedata*128 + 128));



figure
for i = 1:16;
    subplot(4,4,i);
    imagesc(moviedata(:,:,i)',[0 255]); colormap gray; axis equal
end

figure
for i = 1:length(moviedata)/40
imshow(moviedata(:,:,i));
drawnow
end
save metamask2sf2theta4soa15min moviedata xpos ypos maskDur targDur sf lag theta dOri framerate duration  nx ny sigma radius
