clear all
close all

duration = 1;
framerate = 60;
isi = 1.5;
% sfrange = [0 0.04 0.16];
tfrange=0;
sfrange = [0.04 0.16];


phaserange = linspace(0, 2*pi,6);
phaserange=phaserange(1:2);
phaserange = [0 pi]
ntheta = 12;
nx = 3; ny =2;
sigma=20;
randomOrder=1;
randomTheta=0;
binarize=0;
blank=0;

totalDuration = (duration+isi)*length(sfrange)*length(tfrange)*length(phaserange)*ntheta*nx*ny/60
% clear all
% duration = 0.5;
% framerate = 60;
% isi = 0.5;
% sfrange = [0 0.02 0.04 0.16 0.32];
% tfrange =[0 1.5 8];
% phaserange = [0 pi];
% ntheta =10;
% nx = 1; ny =1;
% sigma=0.075;
% randomOrder=1;
% randomTheta=0;

% clear all
% duration = 0.5;
% framerate = 60;
% isi = 0.5;
% sfrange = [0 0.03 0.024];
% tfrange =[0 1.5 8];
% phaserange = [0 pi];
% ntheta =8;
% nx = 2; ny =1;
% sigma=0.1;
% randomOrder=1;
% randomTheta=0;



xsz = 128;
ysz = xsz*9/16;
dist = 25;
width = 50;
height=30;
heightdeg = 2*atand(0.5*height/dist)
degperpix = heightdeg/ysz


blockwidth = floor(0.9*ysz/ny);

xpos = linspace(1,xsz,nx+2);
xposrange = round(xpos(2:end-1)-blockwidth/2);
xposrange = [xposrange 10000]

ypos = linspace(1+ysz/10,0.9*ysz,ny+1);
%ypos = linspace(1,ysz,ny+1);
yposrange = round(ypos(1:end-1));





thetarange = linspace(0, pi,ntheta+1);
if blank
ntheta=ntheta+1;
else
    thetarange = thetarange(1:end-1);
end
trial=0;
for n= 1:length(thetarange);
    for i = 1:length(xposrange);
        for j = 1:length(yposrange);
            for k = 1:length(sfrange);
                for l = 1:length(tfrange);
                    for m= 1:length(phaserange);
                        
                        trial = trial+1;
                        xpos(trial) = xposrange(i); ypos(trial)=yposrange(j);
                        sf(trial)=sfrange(k); tf(trial)=tfrange(l);
                        phase(trial) = phaserange(m); theta(trial) = thetarange(n);
                        if randomTheta
                            theta(trial) = rand*2*pi;
                        end
                    end
                end
            end
        end
    end
end

if randomOrder
order = randperm(trial);
xpos = xpos(order); ypos=ypos(order); sf =sf(order); tf=tf(order); phase=phase(order); theta=theta(order);
end

trial*duration/60

[x y] =meshgrid(1:blockwidth,1:blockwidth);
xgrid=(x-mean(x(:)))/max(x(:));; ygrid=(y-mean(y(:)))/max(y(:));
x = x*degperpix; y=y*degperpix;
gaussian = (exp(-0.5*(xgrid.^2 +ygrid.^2)/sigma^2));
r= 2*sqrt(xgrid.^2 + ygrid.^2);
mask = ones(size(xgrid));
mask1 = mask.*cos(r*pi/2);  mask1(r>=1)=0;
mask2 = mask.*cos((2*r-1)*pi/2);  mask2(r>1)=0; mask2(r<0.5)=1;

figure
imagesc(mask1)
gaussian = gaussian.*mask;
figure
imagesc(gaussian,[0 1]);
if binarize
    gaussian(gaussian<exp(-1))=0;
end

[xfull yfull] = meshgrid(1:ysz,1:xsz);

moviedata = zeros(xsz,ysz,round(trial*(duration+isi)*framerate),'uint8')+128;
for tr = 1:trial
    tr
  
    bkg_ph = (xfull*cos(theta(tr)) + yfull*sin(theta(tr)))*2*pi*sf(tr) + phase(tr);
    ph = (x*cos(theta(tr)+pi/2) + y*sin(theta(tr)+pi/2))*2*pi*sf(tr) + phase(tr);
    if theta(tr)~=2*pi
    for t = 1:duration*framerate;
        tpt = round((tr-1)*duration*framerate +tr*isi*framerate+t);
        bk = cos(bkg_ph);
       if sf(tr)<0.08
           gaussian=mask1;
       else
           gaussian=mask2;
       end
       moviedata(:,:,tpt) = uint8(0.5*256*(bk+1));
     if xpos(tr)<1000
         frame = uint8(0.5*256*(cos(ph).*gaussian +bk(xpos(tr):xpos(tr)+blockwidth-1, ypos(tr):ypos(tr)+blockwidth-1).*(1-gaussian)+1));
 
        moviedata(xpos(tr):xpos(tr)+blockwidth-1, ypos(tr):ypos(tr)+blockwidth-1,tpt) = frame;
     end
    end
    end
end
if binarize
moviedata(moviedata>128)=255;
moviedata(moviedata<128)=0;
end

moviedata = moviedata(1:xsz,1:ysz,:);
mapMonitor

save background3x2y2sf_021215_16minBLANK moviedata xpos ypos tf sf phase theta framerate duration isi nx ny sigma
 
figure
for i = 1:length(moviedata)/20
imshow(moviedata(:,:,i));
drawnow
end

figure
imshow(moviedata(:,:,50))


