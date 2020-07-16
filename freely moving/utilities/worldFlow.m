%%% read in original movie file if not passed as argument
if ~exist('movieFilename','var');
    [f p] = uigetfile('*.avi','world video');
    movieFilename = fullfile(p,f);
end

%%% read in movie
TempVidT = VideoReader(movieFilename);
frame=1; k=1;
display('reading')

while hasFrame(TempVidT)
    worldLaced(:,:,:,frame) = (readFrame(TempVidT));
    
    %%% status update
    if mod(frame,500)==0
        fprintf('frame = %d\n',frame)
    end
    frame=frame+1;   
end

doubleSize = 0;
if doubleSize
world = zeros(size(worldLaced,1), size(worldLaced,2), 3, size(worldLaced,4)*2,'uint8');
else
    world = zeros(size(worldLaced,1)/2, size(worldLaced,2)/2, 3, size(worldLaced,4)*2,'uint8');
end
sz = [size(world,1) size(world,2)];
world(:,:,:,1:2:end) = imresize(worldLaced(1:2:end,:,:,:),sz);
world(:,:,:,2:2:end) = imresize(worldLaced(2:2:end,:,:,:),sz);

worldGrey = (squeeze(mean(imresize(world,0.25),3)));

figure
opticFlow = opticalFlowHS;
for i = 1:size(worldGrey,3);

    flow(i) = estimateFlow(opticFlow,worldGrey(:,:,i));
%subplot(10,10,i);
%imagesc(flow(i).Magnitude)
%imagesc(worldGrey(:,:,i))
%axis off
flowMag(:,:,i) = flow(i).Magnitude;
flowX(:,:,i) = flow(i).Vx;
flowY(:,:,i) = flow(i).Vy;
flowTheta(:,:,i) = flow(i).Orientation;
end

figure
imagesc(mean(flowMag,3))

figure
imagesc(mean(abs(flowX),3))

figure
imagesc(mean(abs(flowY),3))


figure
imagesc(mean(flowX,3))

figure
imagesc(mean(flowY,3))

shortMov = worldGrey;
data = reshape(shortMov,size(shortMov,1)*size(shortMov,2),size(shortMov,3));

c = corrcoef(data);

figure
imagesc(c(300:600,300:600))


figure
[xc lags] = xcorr(squeeze(worldGrey(50,50,:))-92,squeeze(worldGrey(50,50,:))-92,120,'coeff');
plot(lags/60,xc)

figure
x = squeeze(worldGrey(50,50,:)-92);
fs = 60;
y = fft(x);
n = length(x);          % number of samples
f = (0:n-1)*(fs/n);     % frequency range
power = abs(y).^2/n;  

loglog(f,power);

plot(abs(fft(squeeze(worldGrey(50,50,:)-92))))

figure
quiver(flowX(2:2:end,2:2:end,50),flowY(2:2:end,2:2:end,50),2)





