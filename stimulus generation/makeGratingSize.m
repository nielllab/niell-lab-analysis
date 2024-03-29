clear all

% xsz = 128;
xsz = 128;
ysz = xsz*9/16;
dist = 20;
width = 50;
widthdeg = 2*atand(0.5*width/dist)
degperpix = widthdeg/xsz

duration = 0.5;
framerate = 60;
isi = 1.5;
nreps = 1; %number of movie repetitions
% sfrange = [0 0.04 0.16];
% tfrange =[0 2 8];
sfrange = [0.08];
tfrange = [2];
% sfrange = [0.01 0.04 0.16];
% tfrange =[1 4];
% radiusRange = [0 1 2 4 8 1000];
% sizeVals = [0 5 10 20 30 40 50 60];
% sizeVals = [30];
% sizeVals = [0 10 20 30 40 50];
% radiusRange = sizeVals/(8*degperpix);
% radiusRange = [0 2.5 5 10 15 20 25];
radiusRange = [0 1.5 2.5 5 15];
phaserange = linspace(0,2*pi,17); phaserange =phaserange(1:16);  %%%cmn
% phaserange = [0];
%contrastRange = [0.125 0.25 0.5 1];
 contrastRange = [1];


ntheta = 4; %cmn
nx = 1; ny =1;

randomOrder=1;
randomTheta=0;
randomPhase=0;
binarize=0;
blank=0;

totalduration = length(sfrange)*length(contrastRange)*length(tfrange)*length(radiusRange)*length(phaserange)*ntheta*nx*(isi+duration);
sprintf('total duration %0.2f min',totalduration/60)

blockwidth = ysz;
% xposrange = [1 xsz - blockwidth];
xposrange = xsz/4;
yposrange = 1;

[x y] =meshgrid(1:blockwidth,1:blockwidth);
xgrid=(x-mean(x(:))); ygrid=(y-mean(y(:)));
x = xgrid*degperpix; y=ygrid*degperpix;

for r=1:length(radiusRange)
    %%%for use if sizes are in pixels
    % g = (exp(-0.5*(x.^2 +y.^2)/radiusRange(r)^2));
    % g(g<1/128)=0;
    %%%for use if sizes are in degrees
%     g = (exp(-0.5*(x.^2 +y.^2)/(radiusRange(r))^2));
%     g(g<1/2)=0;


%     %alternative for disk
    g = zeros(size(x));
    g(x.^2+y.^2<radiusRange(r).^2) = 1;

    %%% option for disk rather than gaussian envelope
    %%% g = sqrt(x.^2 + y.^2)<radiusRange(r);
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
                        for o = 1:length(contrastRange)
                        
                            trial = trial+1;
                            xpos(trial) = xposrange(i); ypos(trial)=yposrange(1); radius(trial)=j;
                            sf(trial)=sfrange(k); tf(trial)=tfrange(l);
                            phase(trial) = phaserange(m); theta(trial) = thetarange(n);
                            contrasts(trial) = contrastRange(o);
                            if randomTheta
                                theta(trial) = rand*2*pi;
                            end
                            if randomPhase
                                phase(trial) = rand*pi;
                            end
                        end
                    end
                end
            end
        end
    end
end

if randomOrder
order = randperm(trial*nreps);
xpos=repmat(xpos,[1,nreps]);ypos=repmat(ypos,[1,nreps]);sf=repmat(sf,[1,nreps]);tf=repmat(tf,[1,nreps]);
phase=repmat(phase,[1,nreps]);theta=repmat(theta,[1,nreps]);radius=repmat(radius,[1,nreps]);contrasts=repmat(contrasts,[1,nreps]);
xpos = xpos(order); ypos=ypos(order); sf =sf(order); tf=tf(order); phase=phase(order); theta=theta(order);radius = radius(order); contrasts = contrasts(order);
end



moviedata = zeros(xsz,ysz,trial*(duration+isi)*framerate*nreps,'uint8')+128;
for tr = 1:trial*nreps

    ph = (x*cos(theta(tr)) + y*sin(theta(tr)))*2*pi*sf(tr) + phase(tr);
    if theta(tr)~=2*pi
    for t = 1:duration*framerate;
        frame = uint8(0.5*255*(cos(ph + 2*pi*t*tf(tr)/framerate).*gaussian{radius(tr)}*contrasts(tr)+1));
        moviedata(xpos(tr):xpos(tr)+blockwidth-1, ypos(tr):ypos(tr)+blockwidth-1,(tr-1)*duration*framerate +tr*isi*framerate+t) = frame;
    end
    end
end

if binarize
moviedata(moviedata>128)=255;
moviedata(moviedata<128)=0;
end

moviedata = moviedata(1:xsz,1:ysz,:);

% % %%%preview 1/10 of movie
% figure
% for i = 1:length(moviedata)/10
%     i
% imshow(moviedata(:,:,i));
% drawnow
% end
% 
% %%%play a single stimulus
% figure
% stim = 13;
% for i = 31+60*(stim-1):60+60*(stim-1)
%     imshow(moviedata(:,:,i))
%     drawnow
% end
% 
% %%%look at frame i - (i-1) to see change
% i=i+1;
% imshow(moviedata(:,:,i)-moviedata(:,:,i-1))
% drawnow
% 
% %%%look at single frame
% imshow(moviedata(:,:,i))
% drawnow

% save sizeSelect2sf8sz26min moviedata xpos ypos tf sf phase theta framerate duration isi nx ny radius radiusRange contrasts order totalduration sizeVals
%save sizeselect4Ctrst24min moviedata xpos ypos tf sf phase theta framerate duration isi nx ny radius radiusRange contrasts order totalduration -v7.3
save sizeselect4size_blank_05dur_15isi_10min moviedata xpos ypos tf sf phase theta framerate duration isi nx ny radius radiusRange contrasts order totalduration -v7.3
