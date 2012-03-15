ny = 128;
nx = 256;

    MovieMag=6;                 %% magnification that movie will be played at
    screenWidthPix = 1280        %% Screen width in Pixels
    screenWidthCm = 50;         %% Width in cm
    screenDistanceCm = 25;      %% Distance in cm

    screenWidthDeg = 2*atan(0.5*screenWidthCm/screenDistanceCm)*180/pi
    degperpix = (screenWidthDeg/screenWidthPix)*MovieMag


sz = [1 2 4 8 16];

density = 0.2;

flash_duration = 0.25;
MovieRate = 1/flash_duration
frame_duration = 1/MovieRate;
total_duration = 600;
flash_frames = total_duration/flash_duration;

for s= 1:length(sz);
    f=fspecial('disk',round(sz(s)/degperpix));
    f(f>0)=1;
    filt{s}=f;
    area(s) = sum(sum(f));
end


total_area = sum(area);
each_area = (1/length(sz))*density*nx*ny;
n = each_area./area
%n_perframe = n*frame_duration/frame_rate;
mov = zeros(nx,ny,flash_frames);
sz_mov = mov;

for f = 1:flash_frames;
    f
    mov_frm = zeros(nx,ny);
    sz_frm = zeros(nx,ny);
    for s= 1:length(sz);
       
        filt = fspecial('disk', round(sz(s)/degperpix));
%         filt = filt/(max(max(filt)));
        filt(filt>0)=1;
        frm = zeros(nx,ny);
        
        for i = 1:poissrnd(n(s));
            
            frm(ceil(nx*rand(1)),ceil(ny*rand(1))) = round(rand(1))*2 -1;
        end
        %frm_img = imdilate(sz_frm,se);
        frm_img = imfilter(frm,filt,'same');
        mov_frm(frm_img~=0) = frm_img(frm_img~=0);
        sz_frm(frm_img~=0) = sz(s);
    end
    mov(:,:,f) = mov_frm;
    sz_mov(:,:,f) = sz_frm;
end
figure
imagesc(mov(:,:,1))

figure
imagesc(sz_mov(:,:,1))
        

figure
imagesc(mean(mov,3))

figure
imagesc(mean(sz_mov,3));

nan_sz_mov= sz_mov;
nan_sz_mov(nan_sz_mov==0)= NaN;
nan_sz_mov = log2(nan_sz_mov);

n=0;
mov_avg=0;
sz_avg=0;
f=0;
for f = 1:flash_frames;
    if mov(100,100,f)>0 & sz_mov(100,100,f)==1
        n=n+1;
        f_all(n) = f;
        
    end
end
figure
imagesc(mean(mov(:,:,f_all),3));
figure
imagesc(nanmean(nan_sz_mov(:,:,f_all),3));

moviedata = uint8((mov+1)*128-1);
[fname pname]=uiputfile('*.mat');
save(fullfile(pname,fname),'moviedata','sz_mov','degperpix','MovieMag','MovieRate');

            


