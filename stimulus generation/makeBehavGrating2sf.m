clear all
close all

duration = 2;
framerate = 60;
isi = 1;
% sfrange = [0 0.04 0.16];
% tfrange =[0 2 8];
sfrange = [ 0.057 0.114];
sfrange = [0.1067 0.2134];
tfrange =[ 0 ];
%phaserange = [0 0 0 0 0];
phaserange = linspace(0, 2*pi,21);
phaserange=phaserange(1:20);

ntheta =2;
nx = 2; ny =1;
sigma=0.015;
randomOrder=1;
randomTheta=0;
binarize=1;
blank=0;

totalduration = (duration+isi)*length(sfrange)*length(tfrange)*length(phaserange)*ntheta*nx*ny
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
xgrid=(x-mean(x(:))); ygrid=(y-mean(y(:)));
gaussian = sqrt((xgrid.^2 +ygrid.^2))<blockwidth/2;

figure
imagesc(gaussian)

moviedata = zeros(xsz,ysz,trial*(duration+isi)*framerate,'uint8')+128;
for tr = 1:trial
    tr
    ph = (x*cos(theta(tr)) + y*sin(theta(tr)))*2*pi*sf(tr) + phase(tr);
    if theta(tr)~=2*pi
    for t = 1:duration*framerate;
        frame = uint8(0.5*255*(cos(ph + 2*pi*t*tf(tr)/framerate).*gaussian+1));
        moviedata(xpos(tr):xpos(tr)+blockwidth-1, ypos(tr):ypos(tr)+blockwidth-1,(tr-1)*duration*framerate +tr*isi*framerate+t) = frame;
    end
    end
end
if binarize
moviedata(moviedata>128)=255;
moviedata(moviedata<128)=0;
end

moviedata = moviedata(1:xsz,1:ysz,:);
% figure
% for i = 1:length(moviedata)/40
% imshow(moviedata(:,:,i));
% drawnow
% end
save behavStim2sf8min moviedata xpos ypos tf sf phase theta framerate duration isi nx ny sigma
size(moviedata)
length(moviedata)/(60)
figure
for rep = 1:12;
    subplot(4,3,rep);
    imshow(moviedata(:,:,framerate*(isi + rep*(isi+duration))+1))
end


figure
imshow(max(moviedata,[],3))

