%PRLP 10/29/2015 Niell Lab
%This code is used to generate simple unidirectional stimuli based on
%images made in Adobe Illustrator

nstim = 1; %choose number of repeats per stimulus
MovieMag=15;                 %% magnification that movie will be played at
screenWidthPix = 1280;        %% Screen width in Pixels
screenWidthCm = 50;         %% Width in cm
screenDistanceCm = 25;      %% Distance in cm

screenWidthDeg = 2*atan(0.5*screenWidthCm/screenDistanceCm)*180/pi;
degperpix = (screenWidthDeg/screenWidthPix)*MovieMag;

%load stimuli - 4 shapes, 2 sizes each, centered
triBig = imread('geometric shapes-05.jpg');
triSml = imread('geometric shapes-40.jpg');
sqrBig = imread('geometric shapes-15.jpg');
sqrSml = imread('geometric shapes-30.jpg');
crcBig = imread('geometric shapes-10.jpg');
crcSml = imread('geometric shapes-35.jpg');
starBig = imread('geometric shapes-20.jpg');
starSml = imread('geometric shapes-25.jpg');

%load basic stimuli into one array
shapes = zeros(256,256,8);
shapes(:,:,1) = triBig;
shapes(:,:,2) = triSml;
shapes(:,:,3) = sqrBig;
shapes(:,:,4) = sqrSml;
shapes(:,:,5) = crcBig;
shapes(:,:,6) = crcSml;
shapes(:,:,7) = starBig;
shapes(:,:,8) = starSml;

scaleDownBy = 2;
shapes = shapes(1:scaleDownBy:256,1:scaleDownBy:256,:);
stimframes = 30;

%create array to hold all generated stimuli
stimlib = 127*ones(size(shapes,1),size(shapes,2),stimframes,144);
asp = stimframes/2; %used to control amount of motion/duration of frame
ISI = size(stimlib,3);
leg = zeros(144,3); %legend containing info about each stimulus for later analysis

cnt = 1;
for a = 1:size(shapes,3) %for all starting stimuli
    for color = 1:2
        if color==1 %one black and one white, both on grey
            stim = shapes(:,:,a);
            stim(stim<=127) = 127; stim(stim>127) = 0;
            for direc = 1:9 % 9 directions (5 stationary, 4 moving out from center)
                if direc==1
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = stim;
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==2;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,-asp,1);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==3;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,asp,2);
                    end        
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==4;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,asp,1);
                    end 
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==5;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,-asp,2);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==6;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,asp-n,1);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==7;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,-asp+n,2);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==8;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,-asp+n,1);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==9;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,asp-n,2);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                end
            end
        else 
            stim = shapes(:,:,a);
            stim(stim<=127) = 127; stim(stim>127) = 255;
            for direc = 1:9
                if direc==1
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = stim;
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==2;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,-asp,1);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==3;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,asp,2);
                    end        
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==4;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,asp,1);
                    end 
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==5;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,-asp,2);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==6;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,asp-n,1);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==7;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,-asp+n,2);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                 elseif direc==8;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,-asp+n,1);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                elseif direc==9;
                    for n = 1:stimframes
                        stimlib(:,:,n,cnt) = circshift(stim,asp-n,2);
                    end
                    leg(cnt,1) = cnt;
                    leg(cnt,2) = a;
                    leg(cnt,3) = color;
                    leg(cnt,4) = direc;
                    cnt = cnt+1;
                end
            end
        end
    end
end

leg = repmat(leg,nstim,1); %copy * nstim
stimlib = repmat(stimlib,1,1,1,nstim); %copy * nstim

%generate stimulus presentation order
order = randperm(size(leg,1));
ordleg = leg(order,:);
ordstimlib = stimlib(:,:,:,order);

% dir = 'C:\Users\nlab\Desktop\Stimuli';
% nam = 'GeomStim';
% save(fullfile(dir,nam),'ordstimlib','ordleg','nstim','shapes','degperpix','-v7.3');
 
