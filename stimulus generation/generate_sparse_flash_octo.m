clear all
ny = 192;
nx = 256;

MovieMag=6;                 %% magnification that movie will be played at
screenWidthPix = 1280        %% Screen width in Pixels
screenWidthCm = 10  ;         %% Width in cm
screenDistanceCm = 2;      %% Distance in cm

screenWidthDeg = 2*atan(0.5*screenWidthCm/screenDistanceCm)*180/pi
degperpix = (screenWidthDeg/screenWidthPix)*MovieMag


sz = [ 2 4 8 ];
p_flash=0.01;  %%% prob of full-field flash
density = 0.1;

flash_duration = 1;
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
    
    r = rand(1);
    if r<p_flash
        mov_frm(:,:)=1;
        sz_frm=255;
    elseif r<2*p_flash
        mov_frm(:,:)=-1;
        sz_frm = 255;
    else
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
    end
    mov(:,:,f) = mov_frm;
    sz_mov(:,:,f) = sz_frm;
end
figure
imagesc(mov(:,:,1))

figure
imshow(sz_mov(:,:,1))


figure
imagesc(mean(mov,3))

figure
imagesc(median(sz_mov,3));

moviedata = uint8((mov+1)*128-1);


[fname pname]=uiputfile('*.mat');
save(fullfile(pname,fname),'moviedata','sz_mov','degperpix','MovieMag','MovieRate');

figure
imshow(moviedata(:,:,2))


