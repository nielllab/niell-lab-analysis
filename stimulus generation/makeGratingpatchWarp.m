clear all


duration = 2;
framerate = 60;
isi = 2;
% sfrange = [0 0.04 0.16];
% tfrange =[0 2 8];
sfrange = [ 0 0.04  0.16 ];

tfrange =[ 1 8 ];
phaserange = linspace(0, 2*pi,6);
phaserange=phaserange(1:2);
phaserange = [0 pi]
ntheta = 4;
nx = 3; ny =2;
sigma=0.5;
randomOrder=1;
randomTheta=0;
binarize=0;
blank=1;

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


blockwidth = floor(ysz/ny);

xpos = linspace(1,xsz,nx+2);
xposrange = round(xpos(2:end-1)-blockwidth/2);

ypos = linspace(1+ysz/20,ysz-ysz/10,ny+1);
%ypos = linspace(1,ysz,ny+1);
yposrange = round(ypos(1:end-1));





thetarange = linspace(0, 2*pi,ntheta+1);
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
%mask = mask.*cos((2*r-1)*pi/2);  mask(r>1)=0; mask(r<0.5)=1;
mask = mask.*cos(r*pi/2);  mask(r>=1)=0;
figure
imagesc(mask)
gaussian = gaussian.*mask;
figure
imagesc(gaussian,[0 1]);
if binarize
    gaussian(gaussian<exp(-1))=0;
end


moviedata = zeros(xsz,ysz,round(trial*(duration+isi)*framerate),'uint8')+128;
for tr = 1:trial
    tr
    ph = (x*cos(theta(tr)) + y*sin(theta(tr)))*2*pi*sf(tr) + phase(tr);
    if theta(tr)~=2*pi
    for t = 1:duration*framerate;
        frame = uint8(0.5*256*(cos(ph + 2*pi*t*tf(tr)/framerate).*gaussian+1));
        tpt = round((tr-1)*duration*framerate +tr*isi*framerate+t);
        moviedata(xpos(tr):xpos(tr)+blockwidth-1, ypos(tr):ypos(tr)+blockwidth-1,tpt) = frame;
    end
    end
end
if binarize
moviedata(moviedata>128)=255;
moviedata(moviedata<128)=0;
end

moviedata = moviedata(1:xsz,1:ysz,:);
mapMonitor


figure
for i = 1:length(moviedata)/50
imshow(moviedata(:,:,i));
drawnow
end

figure
imshow(moviedata(:,:,50))

save gratingPatch2p24min moviedata xpos ypos tf sf phase theta framerate duration isi nx ny sigma
