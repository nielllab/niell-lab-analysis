clear; close all;

% xsz = 128;
xsz = 128;
ysz = xsz*9/16;
dist = 20;
width = 50;
widthdeg = 2*atand(0.5*width/dist)
degperpix = widthdeg/xsz

duration = 1;
framerate = 60;
isi = 1;

sfrange = [0.16];
tfrange = [0];
% radiusRangeCent = [0 5 10];
% radiusRangeSurr = [0  25 25 25];
% totalRangeCS = [0 10 15 10 15 10 15;0 0 0 25 25 25 25];
% CenterOn = [0 10 15 10 15 0 0];
% SurroundOn = [0 0 0 25 25 25 25];
totalRangeCS = [0 15 15 15;0 0 25 25];
CenterOn = [0 15 15 0];
SurroundOn = [0 0 25 25];
phaserange = linspace(0,2*pi,30); phaserange=phaserange(1:30);  %%%cmn
contrastRange = [1];
nthetaCent = 2;
nthetaSurr = 2;
nx = 1; ny =1;

randomOrder=1;
randomTheta=0;
randomPhase=0;
binarize=1;
blank=0;

totalduration = length(sfrange)*length(contrastRange)*length(tfrange)*length(totalRangeCS)*length(phaserange)*nthetaCent*nthetaSurr*nx*(isi+duration);
% totalduration = length(sfrange)*length(contrastRange)*length(tfrange)*length(radiusRangeCent)*length(radiusRangeSurr)*length(phaserange)*nthetaCent*nthetaSurr*nx*(isi+duration);
sprintf('total duration %0.2f min',totalduration/60)
numtr = length(sfrange)*length(contrastRange)*length(tfrange)*length(totalRangeCS)*length(phaserange)*nthetaCent*nthetaSurr*nx;

blockwidth = ysz;
% xposrange = [1 xsz - blockwidth];
xposrange = xsz/4;
yposrange = 1;

[x y] =meshgrid(1:blockwidth,1:blockwidth);
xgrid=(x-mean(x(:))); ygrid=(y-mean(y(:)));
x = xgrid*degperpix; y=ygrid*degperpix;

for r=1:length(totalRangeCS)
    %%%for use if sizes are in pixels
    % g = (exp(-0.5*(x.^2 +y.^2)/radiusRange(r)^2));
    % g(g<1/128)=0;
    %%%for use if sizes are in degrees
    %     g = (exp(-0.5*(x.^2 +y.^2)/(radiusRange(r))^2));
    %     g(g<1/2)=0;

%     %alternative for disk
    gCenter = zeros(size(x));
    gCenter(x.^2+y.^2 < CenterOn(1,r).^2 ) = 1;
    gSurr = zeros(size(x));
    gSurr(x.^2+y.^2 > totalRangeCS(1,r).^2 & x.^2+y.^2 < totalRangeCS(2,r).^2) = 1;

    %%% option for disk rather than gaussian envelope
    %%% g = sqrt(x.^2 + y.^2)<radiusRange(r);
    gaussianCent{r}=gCenter;
    gaussianSurr{r}=gSurr;
%     figure
%     imagesc(gaussianCent{r}+gaussianSurr{r},[0 1]);
end
%%
% for r=1:length(radiusRangeCent)
%     %%%for use if sizes are in pixels
%     % g = (exp(-0.5*(x.^2 +y.^2)/radiusRange(r)^2));
%     % g(g<1/128)=0;
%     %%%for use if sizes are in degrees
%     %     g = (exp(-0.5*(x.^2 +y.^2)/(radiusRange(r))^2));
%     %     g(g<1/2)=0;
% 
% %     %alternative for disk
%     g = zeros(size(x));
%     g(x.^2+y.^2<radiusRangeCent(r).^2  ) = 1;
% 
%     %%% option for disk rather than gaussian envelope
%     %%% g = sqrt(x.^2 + y.^2)<radiusRange(r);
%     gaussianCent{r}=g;
%     figure
%     imagesc(gaussianCent{r},[0 1]);
% end
% 
% for r=1:length(radiusRangeSurr)
%     %%%for use if sizes are in pixels
%     % g = (exp(-0.5*(x.^2 +y.^2)/radiusRange(r)^2));
%     % g(g<1/128)=0;
%     %%%for use if sizes are in degrees
% %     g = (exp(-0.5*(x.^2 +y.^2)/(radiusRange(r))^2));
% %     g(g<1/2)=0;
% 
% 
% %     %alternative for disk
%     g = zeros(size(x));
%     g(x.^2+y.^2 <radiusRangeSurr(r).^2 & x.^2+y.^2>radiusRangeCent(r+1).^2) = 1;
% 
%     %%% option for disk rather than gaussian envelope
%     %%% g = sqrt(x.^2 + y.^2)<radiusRange(r);
%     gaussianSurr{r}=g;
%     figure
%     imagesc(gaussianSurr{r},[0 1]);
% end

thetarangeCent = linspace(0, pi,nthetaCent+1);
thetarangeCent = thetarangeCent(1:end-1);
thetarangeSurr = linspace(0, pi,nthetaSurr+1);
thetarangeSurr = (thetarangeSurr(1:end-1));

innerSurr = zeros(1,numtr);
trial=0;
for n= 1:length(thetarangeCent)
    for i = 1:length(xposrange)
        for j = 1:length(totalRangeCS) %length(radiusRangeCent);
            for k = 1:length(sfrange)
                for l = 1:length(tfrange)
                    for m= 1:length(phaserange)
                        for o = 1:length(contrastRange)
                            for p = 1:length(thetarangeSurr)
%                                 for q = 1:length(radiusRangeSurr)
                        
                                    trial = trial+1;
                                    xpos(trial) = xposrange(i); ypos(trial)=yposrange(1); radiusCent(trial)=j;
                                    sf(trial)=sfrange(k); tf(trial)=tfrange(l);
                                    phase(trial) = phaserange(m); thetaCent(trial) = thetarangeCent(n);
                                    contrasts(trial) = contrastRange(o); radiusSurr(trial)=j;
                                    if randomTheta
                                        thetaCent(trial) = rand*2*pi;
                                        thetaSurr(trial) = rand*2*pi;
                                    end
                                    if randomPhase
                                        phase(trial) = rand*pi;
                                    end
%                                     if j < floor(length(totalRangeCS)/2)
                                        thetaSurr(trial) = thetarangeSurr(p);
%                                     else
%                                         thetaSurr(trial) = thetarangeSurr(2);
%                                     end
                                    trialID(:,trial) = [j;n;p];
%                                         
%                                 end
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
xpos = xpos(order); ypos=ypos(order); sf =sf(order); tf=tf(order); phase=phase(order); thetaCent=thetaCent(order);thetaSurr=thetaSurr(order);radiusCent = radiusCent(order);radiusSurr = radiusSurr(order); contrasts = contrasts(order);
trialID = trialID(:,order);
end

moviedata = zeros(xsz,ysz,trial*(duration+isi)*framerate,'uint8')+128;
for tr = 1:trial
    tr
    phCent = (x*cos(thetaCent(tr)) + y*sin(thetaCent(tr)))*2*pi*sf(tr) + phase(tr);
    phSurr = (x*cos(thetaSurr(tr)) + y*sin(thetaSurr(tr)))*2*pi*sf(tr) + phase(tr);
    if thetaCent(tr)~=2*pi
        for t = 1:duration*framerate;
            frame = uint8(0.5*255*(((cos(phCent + 2*pi*t*tf(tr)/framerate).*(gaussianCent{radiusCent(tr)})) + (cos(phSurr + 2*pi*t*tf(tr)/framerate).*(gaussianSurr{radiusSurr(tr)}))   )*contrasts(tr)+1));
            moviedata(xpos(tr):xpos(tr)+blockwidth-1, ypos(tr):ypos(tr)+blockwidth-1,(tr-1)*duration*framerate +tr*isi*framerate+t) = frame;
        end
    end
end

if binarize
moviedata(moviedata>128)=255;
moviedata(moviedata<128)=0;
end

moviedata = moviedata(1:xsz,1:ysz,:);
%%
%%%preview 1/10 of movie
% figure
% for i = 1:length(moviedata)/10
%     i
% imshow(moviedata(:,:,i));
% drawnow
% end
% 
%%%play a single stimulus
% figure
% stim = 10;
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
% imshow(moviedata(:,:,end)) 
% drawnow

% save sizeSelect2sf8sz26min moviedata xpos ypos tf sf phase theta framerate duration isi nx ny radius radiusRange contrasts order totalduration sizeVals
save patchonpatch16min moviedata xpos ypos tf sf phase framerate duration isi nx ny thetaCent thetaSurr radiusCent radiusSurr CenterOn SurroundOn totalRangeCS trialID contrasts order totalduration -v7.3
%%
for ii=1:size(totalRangeCS,2)
%     temp1(ii,1) = sum(trialID(1,:) == 1 & trialID(2,:) == 1 & trialID(3,:)==1);
%     temp2(ii,1) = sum(trialID(1,:) == 1 & trialID(2,:) == 1 & trialID(3,:)==2);
%     temp3(ii,1) = sum(trialID(1,:) == 1 & trialID(2,:) == 2 & trialID(3,:)==1);

    temp1(ii,:) = find(trialID(1,:)==ii & trialID(2,:)==1 & trialID(3,:)==1); %temp1=temp1(1);
    temp2(ii,:) = find(trialID(1,:)==ii & trialID(2,:)==1 & trialID(3,:)==2); %temp1=temp1(1);
    temp3(ii,:) = find(trialID(1,:)==ii & trialID(2,:)==2 & trialID(3,:)==1); %temp1=temp1(1);
    temp4(ii,:) = find(trialID(1,:)==ii & trialID(2,:)==2 & trialID(3,:)==2); %temp1=temp1(1);
end

tempT = ([temp1(:,1) temp2(:,1) temp3(:,1) temp4(:,1)]);
tempT = tempT(:);
moviedata = zeros(xsz,ysz,14,'uint8')+128; i=1;

for ii = 1:size(tempT)
    
    tr=tempT(ii);
    phCent = (x*cos(thetaCent(tr)) + y*sin(thetaCent(tr)))*2*pi*sf(tr) + phase(tr);
    phSurr = (x*cos(thetaSurr(tr)) + y*sin(thetaSurr(tr)))*2*pi*sf(tr) + phase(tr);
%     if thetaCent(tr)~=2*pi || thetaSurr(tr) ~=2*pi
        for t = 1:duration*framerate
            frame = uint8(0.5*255*(((cos(phCent + 2*pi*t*tf(tr)/framerate).*(gaussianCent{radiusCent(tr)})) + (cos(phSurr + 2*pi*t*tf(tr)/framerate).*(gaussianSurr{radiusSurr(tr)}))   )*contrasts(tr)+1));
            moviedata(xpos(tr):xpos(tr)+blockwidth-1, ypos(tr):ypos(tr)+blockwidth-1,i) = frame; 
        end
        if binarize
            moviedata(moviedata>128)=255;
            moviedata(moviedata<128)=0;
        end
%     end
    figure; 
    %imshow(moviedata(:,:,i));
    imagesc(frame);
%     title1 = ['stim' num2str(i) '.jpg'];
%     imwrite(moviedata(:,:,i),title1);
    i=i+1;
end