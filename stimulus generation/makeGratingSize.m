clear all
duration = 1;
framerate = 60;
isi = 1;
% sfrange = [0 0.04 0.16];
% tfrange =[0 2 8];
sfrange = [ 0.04 0.16];
tfrange =[ 0 2];
radiusRange = [0 1 2 4 8 1000];
phaserange = [0 pi]

ntheta = 4;
nx = 2; ny =1;

randomOrder=1;
randomTheta=1;
binarize=0;
blank=0;

totalduration = length(sfrange)*length(tfrange)*length(radiusRange)*length(phaserange)*ntheta*nx*(isi+duration);
totalduration/60

xsz = 128;
ysz = xsz*9/16;
dist = 25;
width = 50;
widthdeg = 2*atand(0.5*width/dist)
degperpix = widthdeg/xsz



blockwidth = ysz;
xposrange = [1 xsz - blockwidth];
yposrange  = 1;


[x y] =meshgrid(1:blockwidth,1:blockwidth);
xgrid=(x-mean(x(:)));; ygrid=(y-mean(y(:)));
x = xgrid*degperpix; y=ygrid*degperpix;

for r=1:length(radiusRange)
g = (exp(-0.5*(x.^2 +y.^2)/radiusRange(r)^2));
g(g<1/128)=0;
gaussian{r}=g;
figure
imagesc(gaussian{r},[0 1]);
end

thetarange = linspace(0, 2*pi,ntheta+1);
    thetarange = thetarange(1:end-1);
trial=0;
for n= 1:length(thetarange);
    for i = 1:length(xposrange);
        for j = 1:length(radiusRange);
            for k = 1:length(sfrange);
                for l = 1:length(tfrange);
                    for m= 1:length(phaserange);
                        
                        trial = trial+1;
                        xpos(trial) = xposrange(i); ypos(trial)=yposrange(1); radius(trial)=j;
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
xpos = xpos(order); ypos=ypos(order); sf =sf(order); tf=tf(order); phase=phase(order); theta=theta(order);radius = radius(order);
end



moviedata = zeros(xsz,ysz,trial*(duration+isi)*framerate,'uint8')+128;
for tr = 1:trial
    tr
    ph = (x*cos(theta(tr)) + y*sin(theta(tr)))*2*pi*sf(tr) + phase(tr);
    if theta(tr)~=2*pi
    for t = 1:duration*framerate;
        frame = uint8(0.5*255*(cos(ph + 2*pi*t*tf(tr)/framerate).*gaussian{radius(tr)}+1));
        moviedata(xpos(tr):xpos(tr)+blockwidth-1, ypos(tr):ypos(tr)+blockwidth-1,(tr-1)*duration*framerate +tr*isi*framerate+t) = frame;
    end
    end
end
if binarize
moviedata(moviedata>128)=255;
moviedata(moviedata<128)=0;
end

moviedata = moviedata(1:xsz,1:ysz,:);
figure
for i = 1:length(moviedata)/10
    i
imshow(moviedata(:,:,i));
drawnow
end
save sizeSelect2sf5sz14min moviedata xpos ypos tf sf phase theta framerate duration isi nx ny radius radiusRange