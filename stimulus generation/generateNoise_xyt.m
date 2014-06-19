function moviedata=generateNoise_xyt(maxSpatFreq,maxTempFreq,contrastSigma,duration,movtype,binarize);
%%% generates white noise movies with limited spatial and temporal
%%% frequency, via inverse fourier transform
%%% now with more options ...
%%% binarize (0/1) = convert grayscale to binary, for higher contrast
%%% movie type = allows noise patches, moving noise bars, etc.
%%% movietype 1 = standard contrast modulated noise
%%% movietype 2 = full-field step function of white noise
%%% movietype 3= moving bar of white noise
%%% movietype 4 = alternating patches of white noise along x-axis
%%% movietype 5 = alternating patches of white noise along y-axis

rand('state',sum(100*clock))
tic

%%% stimulus/display parameters

if ~exist('binarize','var') |  isempty(binarize)
    binarize=0;
end
if ~exist('movtype','var') | isempty(movtype)
    movtype=0;
end

cmod=1; step=2; bar = 3;  Xpatches=4; Ypatches=5; %%%% movietypes

imsize = 128;                %% size in pixels
framerate = 30;             %% Hz
imageMag=18;                 %% magnification that movie will be played at

screenWidthPix = 1920        %% Screen width in Pixels
screenWidthCm = 50;         %% Width in cm
screenDistanceCm = 25;      %% Distance in cm

%duration = .5;               %% duration in minutes
%maxSpatFreq = 0.12          %% spatial frequency cutoff (cpd)
%maxTempFreq = 10;          %% temporal frequency cutoff
%contrastSigma =0.5;         %% one-sigma value for contrast

%% derived parameters
screenWidthDeg = 2*atan(0.5*screenWidthCm/screenDistanceCm)*180/pi
degperpix = (screenWidthDeg/screenWidthPix)*imageMag

nframes = framerate*60*duration;

%% frequency intervals for FFT
nyq_pix = 0.5;
nyq_deg=nyq_pix/degperpix;
freqInt_deg = nyq_deg / (0.5*imsize);
freqInt_pix = nyq_pix / (0.5*imsize);
nyq = framerate/2;
tempFreq_int = nyq/(0.5*nframes)


%% cutoffs in terms of frequency intervals

tempCutoff = round(maxTempFreq/tempFreq_int);
maxFreq_pix = maxSpatFreq*degperpix;
spatCutoff = round(maxFreq_pix / freqInt_pix);


%%% generate frequency spectrum (invFFT)
alpha=-1;
offset=3;
range_mult =1;
%for noise that extends past cutoff parameter (i.e. if cutoff = 1sigma)
%range_mult=2;
spaceRange = (imsize/2 - range_mult*spatCutoff : imsize/2 + range_mult*spatCutoff)+1;
tempRange =   (nframes /2 - range_mult*tempCutoff : nframes/2 + range_mult*tempCutoff)+1;
[x y z] = meshgrid(-range_mult*spatCutoff:range_mult*spatCutoff,-range_mult*spatCutoff:range_mult*spatCutoff,-range_mult*tempCutoff:range_mult*tempCutoff);
%% can put any other function to describe frequency spectrum in here,
%% e.g. gaussian spectrum
% use = exp(-1*((0.5*x.^2/spatCutoff^2) + (0.5*y.^2/spatCutoff^2) + (0.5*z.^2/tempCutoff^2)));
%  use =single(((x.^2 + y.^2)<=(spatCutoff^2))& ((z.^2)<(tempCutoff^2)) );
use =single(((x.^2 + y.^2)<=(spatCutoff^2))& ((z.^2)<(tempCutoff^2)) ).*(sqrt(x.^2 + y.^2 +offset).^alpha);
clear x y z;


%%%
invFFT = zeros(imsize,imsize,nframes,'single');
mu = zeros(size(spaceRange,2), size(spaceRange,2), size(tempRange,2));
sig = ones(size(spaceRange,2), size(spaceRange,2), size(tempRange,2));
invFFT(spaceRange, spaceRange, tempRange) = single(use .* normrnd(mu,sig).*exp(2*pi*i*rand(size(spaceRange,2), size(spaceRange,2), size(tempRange,2))));
clear use;

%% in order to get real values for image, need to make spectrum
%% symmetric
fullspace = -range_mult*spatCutoff:range_mult*spatCutoff; halftemp = 1:range_mult*tempCutoff;
halfspace = 1:range_mult*spatCutoff;
invFFT(imsize/2 + fullspace+1, imsize/2+fullspace+1, nframes/2 + halftemp+1) = ...
    conj(invFFT(imsize/2 - fullspace+1, imsize/2-fullspace+1, nframes/2 - halftemp+1));
invFFT(imsize/2+fullspace+1, imsize/2 + halfspace+1,nframes/2+1) = ...
    conj( invFFT(imsize/2-fullspace+1, imsize/2 - halfspace+1,nframes/2+1));
invFFT(imsize/2+halfspace+1, imsize/2 +1,nframes/2+1) = ...
    conj( invFFT(imsize/2-halfspace+1, imsize/2+1,nframes/2+1));

figure
imagesc(abs(invFFT(:,:,nframes/2+1)));
figure
imagesc(angle(invFFT(:,:,nframes/2)));

pack
shiftinvFFT = ifftshift(invFFT);
clear invFFT;

%%% invert FFT and scale it to 0 -255

imraw = real(ifftn(shiftinvFFT));
clear shiftinvFFT;
immean = mean(imraw(:))
immax = std(imraw(:))/contrastSigma
immin = -1*immax
imscaled = (imraw - immin-immean) / (immax - immin);
clear imfiltered;
contrast_period =10;


%%% modify movie for different patterns (patches, bars, etc)
if binarize
    imscaled(imscaled<0.5)=0;
    imscaled(imscaled>0.5)=1;
end

if movtype==bar
    barPix=20;
    center = linspace(0-barPix/2,imsize+barPix/2,contrast_period*framerate);
    center=fliplr(center);
    loweredge = center-barPix/2; upperedge=center+barPix/2;
    loweredge(loweredge<1)=1; loweredge(loweredge>imsize)=imsize;
     upperedge(upperedge<1)=1; upperedge(upperedge>imsize)=imsize;
    loweredge = round(loweredge); upperedge=round(upperedge);
end

for f = 1:nframes
    if movtype == cmod
        imscaled(:,:,f) = (imscaled(:,:,f)-.5).*(0.5-0.5*cos(2*pi*f/(contrast_period*framerate)));
    elseif movtype ==step
        %imscaled(:,:,f) = (imscaled(:,:,f)-.5).*(sin(pi+2*pi*f/(contrast_period*framerate))>0); 
        imscaled(:,:,f) = 0.5*sign(sin(pi+2*pi*f/(contrast_period*framerate))); 

    elseif movtype ==bar
        le = loweredge(mod(f-1,contrast_period*framerate)+1);
        ue = upperedge(mod(f-1,contrast_period*framerate)+1);
        imscaled(1:le,:,f)=0.5; imscaled(ue:imsize,:,f)=0.5;
    elseif movtype==Xpatches | movtype==Ypatches
        if movtype==Xpatches
       % wx = 0.25;  xpos = 0.25; ypos = 0.5;   
       wx = 0.2; xpos = 0.4; ypos=0.5;
        widthpix = 2*round(wx*imsize/2);  yposPix = ypos*imsize;
        if mod((f-1)/(contrast_period*framerate),1)<0.5
            xposPix = xpos*imsize;
        else 
            xposPix=(1-xpos)*imsize
        end
        elseif movtype==Ypatches
          % wx = 0.25;  xpos = 0.5; ypos = 0.35;  
           wx = 0.2;  xpos = 0.5; ypos = 0.4;   
        widthpix = 2*round(wx*imsize/2);  xposPix = xpos*imsize;
        if mod((f-1)/(contrast_period*framerate),1)<0.5
            yposPix = ypos*imsize;
        else 
            yposPix=(1-ypos)*imsize
        end
        
        end
        xposPix=round(xposPix); yposPix=round(yposPix);
        imscaled(1:(xposPix-widthpix/2),:,f)=0.5; imscaled(xposPix+widthpix/2:imsize,:,f)=0.5; 
        imscaled(:,1:(yposPix-widthpix/2),f)=0.5; imscaled(:,(yposPix+widthpix/2):imsize,f)=0.5; 
    end
        
    %imscaled(:,:,f) = (imscaled(:,:,f)-.5).*(contrast(mod(f-1,1800)+1));
end
if movtype ==cmod | movtype==step
imscaled = imscaled+0.5;
end
moviedata = uint8(floor(imscaled(1:imsize,1:imsize,:)*255)+1);

%%%   to check pixel intensity distribution      (slow!)

%     pixdata = single(moviedata);
%     figure
%     hist(pixdata(:));
%     figure
%


%% to check that the spectrum is still correct
clear imscaled
c = fftn(single(moviedata)-128);
c = fftshift(c);
figure
imagesc(mean(abs(c(:,:,:)),3));

MovieRate= framerate;

MovieMag= imageMag;

[bname output_path] = uiputfile('','data folder');
save(fullfile(output_path,bname),'moviedata','MovieMag','MovieRate','maxSpatFreq','maxTempFreq','contrastSigma','alpha','offset','degperpix','screenWidthPix')
%%% to view movie
%
%     for f=1:1000
%
%         imshow(moviedata(:,:,f));
%         mov(f) = getframe(gcf);
%     end
%     toc
%     movie(mov,10,30)