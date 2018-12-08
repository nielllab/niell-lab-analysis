clear all
close all
%%% makes stimuli to test occlusion

%%% original version
%%% cond 1-4 are patch gets covered by moving square
%%% cond 5-8 are moving patch goes behing stationary square
%%% cond 9 = blank

%%% cond 1 = stationary grating patch
%%% cond 2 = moving square
%%% cond 3 = moving square occludes patch
%%% cond 4 = stationary patch blinks off during time it would be occluded

%%% cond 5 = moving patch
%%% cond 6 = stationary square
%%% cond 7=  moving patch goes behind square
%%% cond 8 = moving patch blinks during time it would be behind square

%%% updated version
%%% cond 1-5 are patch gets covered by moving square
%%% cond 6-10 are moving patch goes behing stationary square
%%% cond 11 = blank

% condList{1} = 'stationary grating patch'
% condList{2}= 'moving square'
% condList{3} = 'moving square occludes patch'
% condList{4} = 'stationary patch blinks off during time it would be occluded'
% condList{5}= 'stationary patch gradually occludes, no square'

condList{1} = 'moving patch'
condList{2} = 'stationary square'
condList{3} =  'moving patch goes behind square'
condList{4} = 'moving patch blinks during time it would be behind square'
condList{5} = 'moving patch gradually occludes, no square'
condList{6} = 'blank'

duration = 4;
framerate = 60;
isi = 1;
% sfrange = [0 0.04 0.16];
% tfrange =[0 2 8];
sf = 0.1067;
tf = 2 ;
thetaRange = [pi/4 3*pi/4];
phaseRange = linspace(0,2*pi,4);
phaseRange = phaseRange(1:end-1);

condRange = 1:6;
contRange = [-1 1];

randomOrder=1;
randomTheta=0;
binarize=1;
blank=0;

%patch only, square only, square over patch, patch disappears'
% square moves vvs patch moves
% randome theta (horiz or vertical)
% 4 sec to cross screen
% in stationary patch, tf =2

totalduration = (duration+isi)*length(phaseRange)*length(condRange)*length(contRange);

totalduration/60;


xsz = 256;
ysz = xsz*9/16;
dist = 25;
width = 50;
widthdeg = 2*atand(0.5*width/dist)
degperpix = widthdeg/xsz



patchSz = ysz/4;
[x y] = meshgrid(1:patchSz,1:patchSz);
x = x-mean(x(:)); y= y-mean(y(:));
patchmask = sqrt(x.^2 + y.^2)<(patchSz/2);
occluder = ones(ysz/2,ysz/2);

patchPh = 2*pi*y*degperpix*sf;


trial=0;
   


    for c = 1:length(contRange);  %%% occluder contrast
        for t = 1:length(thetaRange);
            for p = 1:length(phaseRange);
                for i = 1:length(condRange);
                trial = trial+1;
                cond(trial)=condRange(i);
                occContrast(trial) = contRange(c);
                theta(trial) = thetaRange(t);
                phase(trial) = phaseRange(p);
            end            
        end
    end
end

if randomOrder
order = randperm(trial);
cond = cond(order); occContrast = occContrast(order); theta = theta(order);
end

occSz = size(occluder,1);
boxSpd = (xsz-patchSz-1)/(duration*framerate);
centerRangeX = (xsz/2 - patchSz/2 +1) : (xsz/2 +patchSz/2);
centerRangeY = (ysz/2 - patchSz/2 +1) : (ysz/2 +patchSz/2);
boxRangeX = 1:occSz; boxRangeY = (ysz/2 - occSz/2 + 1):(ysz/2 +occSz/2);
boxCenterX = (xsz/2 - occSz/2 +1) : (xsz/2 +occSz/2);
boxCenterY = (ysz/2 - occSz/2 +1) : (ysz/2 +occSz/2);

moviedata = zeros(xsz,ysz,trial*(duration+isi)*framerate,'uint8')+128;
for tr = 1:trial
    tr
    ph = (x*cos(theta(tr)) + y*sin(theta(tr)))*2*pi*sf*degperpix + phase(tr);
    if theta(tr)~=2*pi
        for t = 1:duration*framerate;
%             if cond(tr) ==1
%                 %%% stationary patch
%                 patch= uint8(0.5*255*(sign(cos(ph + 2*pi*t*tf/framerate)).*patchmask+1));
%                 moviedata(centerRangeX,centerRangeY,(tr-1)*duration*framerate +tr*isi*framerate+t) = patch;
%             elseif cond(tr)==2
%                 %%% moving box
%                 boxRangeX = (1:occSz) + ceil(boxSpd*t -patchSz/2);
%                 boxRangeX=boxRangeX(boxRangeX>0); boxRangeX = boxRangeX(boxRangeX<=xsz);
%                 moviedata(boxRangeX,boxRangeY,(tr-1)*duration*framerate +tr*isi*framerate+t) = 0.5*(occContrast(tr)+1)*255;
%             elseif cond(tr)==3
%                 %%% stationary patch
%                 patch= uint8(0.5*255*(sign(cos(ph + 2*pi*t*tf/framerate)).*patchmask+1));
%                 moviedata(centerRangeX,centerRangeY,(tr-1)*duration*framerate +tr*isi*framerate+t) = patch;
%                 %%% moving box
%                  boxRangeX = (1:occSz) + ceil(boxSpd*t -patchSz/2);
%                 boxRangeX=boxRangeX(boxRangeX>0); boxRangeX = boxRangeX(boxRangeX<=xsz);
%                 
%                 moviedata(boxRangeX,boxRangeY,(tr-1)*duration*framerate +tr*isi*framerate+t) = 0.5*(occContrast(tr)+1)*255;
%             elseif cond(tr) ==4;
%                 %%% stationary patch
%                 patch= uint8(0.5*255*(sign(cos(ph + 2*pi*t*tf/framerate)).*patchmask+1));
%                 %%% where the box would be
%                  boxRangeX = (1:occSz) + ceil(boxSpd*t -patchSz/2);
%                 boxRangeX=boxRangeX(boxRangeX>0); boxRangeX = boxRangeX(boxRangeX<=xsz);
%                 
%                 %%% only show stationary patch when box wouldn't be on it
%                 %%% (blinking instead of occlusion)
%                 if max(boxRangeX)<min(centerRangeX) | min(boxRangeX)>max(centerRangeX)
%                     moviedata(centerRangeX,centerRangeY,(tr-1)*duration*framerate +tr*isi*framerate+t) = patch;
%                 end
%                             elseif cond(tr) ==5;
%                       %%% stationary patch
%                 patch= uint8(0.5*255*(sign(cos(ph + 2*pi*t*tf/framerate)).*patchmask+1));
%                 moviedata(centerRangeX,centerRangeY,(tr-1)*duration*framerate +tr*isi*framerate+t) = patch;
%                 %%% moving box -> grey
%                  boxRangeX = (1:occSz) + ceil(boxSpd*t -patchSz/2);
%                 boxRangeX=boxRangeX(boxRangeX>0); boxRangeX = boxRangeX(boxRangeX<=xsz);
%                 
%                 moviedata(boxRangeX,boxRangeY,(tr-1)*duration*framerate +tr*isi*framerate+t) = 128;
%       
           if cond(tr)==1
               %%% moving patch
               patch= uint8(0.5*255*(sign(cos(ph)).*patchmask+1));
                patchRangeX = (1:patchSz) + ceil(boxSpd*t );
                patchRangeX=patchRangeX(patchRangeX>0); patchRangeX = patchRangeX(patchRangeX<=xsz);
                moviedata(patchRangeX,centerRangeY,(tr-1)*duration*framerate +tr*isi*framerate+t) = patch;
            elseif cond(tr) == 2
               %%% stationary box
               moviedata(boxCenterX,boxCenterY,(tr-1)*duration*framerate +tr*isi*framerate+t) = 0.5*(occContrast(tr)+1)*255;
            elseif cond(tr)==3
                %%% moving patch
                patch= uint8(0.5*255*(sign(cos(ph)).*patchmask+1));
                patchRangeX = (1:patchSz) + ceil(boxSpd*t );
                patchRangeX=patchRangeX(patchRangeX>0); patchRangeX = patchRangeX(patchRangeX<=xsz);
                moviedata(patchRangeX,centerRangeY,(tr-1)*duration*framerate +tr*isi*framerate+t) = patch;
               %%% stationary box
               moviedata(boxCenterX,boxCenterY,(tr-1)*duration*framerate +tr*isi*framerate+t) = 0.5*(occContrast(tr)+1)*255;
            elseif cond(tr)==4
                patch= uint8(0.5*255*(sign(cos(ph)).*patchmask+1));
              patchRangeX = (1:patchSz) + ceil(boxSpd*t );
                patchRangeX=patchRangeX(patchRangeX>0); patchRangeX = patchRangeX(patchRangeX<=xsz);
                if max(patchRangeX)<min(boxCenterX) | min(patchRangeX)>max(boxCenterX)
                    moviedata(patchRangeX,centerRangeY,(tr-1)*duration*framerate +tr*isi*framerate+t) = patch;
                end
                         elseif cond(tr)==5
                %%% moving patch
                patch= uint8(0.5*255*(sign(cos(ph)).*patchmask+1));
                patchRangeX = (1:patchSz) + ceil(boxSpd*t );
                patchRangeX=patchRangeX(patchRangeX>0); patchRangeX = patchRangeX(patchRangeX<=xsz);
                moviedata(patchRangeX,centerRangeY,(tr-1)*duration*framerate +tr*isi*framerate+t) = patch;
               %%% stationary box
               moviedata(boxCenterX,boxCenterY,(tr-1)*duration*framerate +tr*isi*framerate+t) =128;
          
            end
            
        end
    end
end


moviedata = moviedata(1:xsz,1:ysz,:);
% figure
% for i = 1:length(moviedata)/40
% imshow(moviedata(:,:,i));
% drawnow
% end



save occludeStim5cond_2theta_6min moviedata  tf sf phase theta framerate duration isi cond occContrast condList -v7.3
size(moviedata)
length(moviedata)/(60)
figure
for rep = 25:36;
    subplot(4,3,rep-24);
    imshow(moviedata(:,:,framerate*(isi + rep*(isi+duration))+1))
end

figure
for fr = 1:5:(15*duration*framerate)
   rep =ceil(fr/((duration+isi)*framerate)) 
   cond(rep)
    imshow(moviedata(:,:,fr)');
 
    title(sprintf('%d  %d  %d', fr,rep,cond(rep)))
       drawnow
end

