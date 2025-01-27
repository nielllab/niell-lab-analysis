clear all
close all

duration = 1;
framerate = 30;
isi = 1.5;
% sfrange = [0 0.04 0.16];
% tfrange =[0 2 8];

sfrange = 0.1067;
tfrange =[ 0 ];
%phaserange = [0 0 0 0 0];
phaserange = linspace(0, 2*pi,9);
phaserange=phaserange(1:2);

ntheta =2;
nx = 2; ny =1;
sigma=0.015;
randomOrder=1;
randomTheta=0;
binarize=0;
blank=0;

contrasts = [0 0 0 0 0 0 0 0 0 1]

totalduration = (duration+isi)*length(sfrange)*length(tfrange)*length(phaserange)*ntheta*nx*ny*length(contrasts)
totalduration/60

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



xsz = 256;
ysz = xsz*9/16;
dist = 25;
width = 50;
widthdeg = 2*atand(0.5*width/dist)
degperpix = widthdeg/xsz
sfrange = sfrange*degperpix;


blockwidth = 100*5/7;

xposrange = [0.33  0.66]*xsz - blockwidth/2
ypos = linspace(1,ysz-1,ny+1);
yposrange = ypos(1:end-1);
yposrange = 0.5*ysz - blockwidth/2;




thetarange = linspace(0, pi,ntheta+1);
if blank
    ntheta=ntheta+1;
else
    thetarange = thetarange(1:end-1);
end

%thetarange = [0  pi/2 ];
trial=0;
for n= 1:length(thetarange);
    for i = 1:length(xposrange);
        for j = 1:length(yposrange);
            for k = 1:length(sfrange);
                for l = 1:length(tfrange);
                    for m= 1:length(phaserange);
                        for c = 1:length(contrasts);
                            trial = trial+1;
                            xpos(trial) = xposrange(i); ypos(trial)=yposrange(j);
                            sf(trial)=sfrange(k); tf(trial)=tfrange(l);
                            phase(trial) = phaserange(m); theta(trial) = thetarange(n);
                            contrast(trial) = contrasts(c);
                            if randomTheta
                                theta(trial) = rand*2*pi;
                            end
                        end
                    end
                end
            end
        end
    end
end

if randomOrder
    order = randperm(trial);
    xpos = xpos(order); ypos=ypos(order); sf =sf(order); tf=tf(order); phase=phase(order); theta=theta(order); contrast = contrast(order);
end

trial*duration/60

[x y] =meshgrid(1:blockwidth,1:blockwidth);
xgrid=(x-mean(x(:))); ygrid=(y-mean(y(:)));
gaussian = sqrt((xgrid.^2 +ygrid.^2))<blockwidth/2;
% 
% sigma = 10
% gaussian = sqrt(exp(-0.5*(xgrid.^2 +ygrid.^2)/sigma^2));

figure
imagesc(gaussian)

%moviedata = zeros(xsz,ysz,trial*(duration+isi)*framerate,'uint8')+128;
moviedata = zeros(xsz,ysz,trial*(duration+isi)*framerate,'uint8');
for tr = 1:trial
    tr
    ph = (x*cos(theta(tr)) + y*sin(theta(tr)))*2*pi*sf(tr) + phase(tr);
    if theta(tr)~=2*pi
        for t = 1:duration*framerate;
            %frame = uint8(0.5*255*(cos(ph + 2*pi*t*tf(tr)/framerate).*gaussian*contrast(tr)+1));
           
            moviedata(:,:,(tr-1)*duration*framerate +tr*isi*framerate+t) = uint8(contrast(tr)*255);
        end
    end
end
if binarize
    moviedata(moviedata>128)=255;
    moviedata(moviedata<128)=0;
end

moviedata = moviedata(1:xsz,1:ysz,:);
% figure
% for i = 1:length(moviedata)/10
% imshow(moviedata(:,:,i));
% drawnow
% end
save flashStim moviedata xpos ypos tf sf phase theta framerate duration isi nx ny sigma contrast
size(moviedata)
length(moviedata)/(60)
figure
for rep = 0:11;
    subplot(4,3,rep+1);
    imshow(moviedata(:,:,framerate*(isi + rep*(isi+duration))+1))
end


figure
imshow(max(moviedata,[],3))

