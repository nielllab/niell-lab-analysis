%[fname pname]=uiputfile('*.mat');

clear all
pack
    MovieMag=6;                 %% magnification that movie will be played at
    screenWidthPix = 1280        %% Screen width in Pixels
    screenWidthCm = 50;         %% Width in cm
    screenDistanceCm = 25;      %% Distance in cm

    screenWidthDeg = 2*atan(0.5*screenWidthCm/screenDistanceCm)*180/pi
    degperpix = (screenWidthDeg/screenWidthPix)*MovieMag


sz = [2 4 8];
speeds = [10 20 40 80 160];

ny = ceil(256 +2*max(sz)/degperpix);
nx = ceil(256+2*max(sz)/degperpix);
nout_x = 256;
nout_y = 128;

density = 0.04;

MovieRate = 30;
frame_duration = 1/MovieRate;
%total_duration=10;
total_duration = 1200;

nframes = total_duration*MovieRate;
for s= 1:length(sz);
    f=fspecial('disk',round(sz(s)/degperpix));
    f(f>0)=1;
    filt{s}=f;
    area(s) = sum(sum(f));
end

area_theory = pi*(sz/degperpix).^2;
total_area = sum(area);
each_area = (1/length(sz))*density*nx*ny;
npersize = each_area./area;
probpersize = npersize/sum(npersize);
cumprob = cumsum(probpersize);

probperspeed = sqrt(speeds)/sum(sqrt(speeds));
cumprobspeed = cumsum(probperspeed);

%n_perframe = n*frame_duration/frame_rate;
moviedata = zeros(nout_x,nout_y,nframes,'single');
sz_mov = uint8(moviedata);
sp_mov =uint8(moviedata);
th_mov=moviedata;

spots=zeros(length(npersize),6);
N=0;
for i = 1:sum(npersize);
    N=N+1;
     spots = addspot(sz,speeds,npersize,N,spots,nx,ny,cumprob,cumprobspeed);
end
spots

   
xrange = nx/2-nout_x/2 : nx/2+nout_x/2-1;
yrange = ny/2-nout_y/2 : ny/2+nout_y/2-1;
tic

for f = 1:nframes
   if round(f/100) == f/100
       sprintf('%0.2f%% done, estimate %0.1f secs to go',100*f/nframes, toc*(nframes-f)/f)
   end
  
   frm = zeros(nx,ny,'int8');
    szfrm = zeros(nx,ny,'uint8');
    speedfrm = zeros(nx,ny);
    thetafrm = zeros(nx,ny);
    for s = 1:N;
        drawspot = zeros(nx,ny,'uint8');
        drawspot(round(spots(s,1)),round(spots(s,2)))=1;       
        drawspot = imfilter(drawspot,filt{spots(s,3)});
        pix = drawspot~=0;
        frm(pix) = double(drawspot(pix))*spots(s,6);
        szfrm(pix)=spots(s,3);
        speedfrm(pix)=spots(s,4);
        thetafrm(pix)=spots(s,5);
        spots(s,1) = spots(s,1)+spots(s,4)*cos(spots(s,5))*frame_duration/degperpix;
         spots(s,2) = spots(s,2)+spots(s,4)*sin(spots(s,5))*frame_duration/degperpix;
        if spots(s,1)<1 | spots(s,2)<2 | spots(s,1)>nx | spots(s,2)>ny
            spots = addspot(sz,speeds,npersize,s,spots,nx,ny,cumprob,cumprobspeed);
        end
    end
% figure(movfig);
% imagesc(thetafrm);

    moviedata(:,:,f)=frm(xrange,yrange );
    sz_mov(:,:,f) = szfrm(xrange,yrange);
    sp_mov(:,:,f) = speedfrm(xrange,yrange);
    th_mov(:,:,f)= thetafrm(xrange,yrange);
end

sz_mov= uint8(sz_mov);
sz_mov_ind=sz_mov;
for i =1:length(sz);
    sz_mov(sz_mov_ind==i)=sz(i);
end
sp_mov = uint8(sp_mov);
clear sz_mov_ind;

th_mov = mod(th_mov,2*pi);
th_mov = uint8(255*th_mov/(2*pi));

% for i = 1:nframes
%     imagesc(sz_mov(:,:,i));
%     getframe(gcf);
% end

h = hist(sp_mov(:),[0 10 20 40 80 160]);
sp_hist = h/sum(h)

h = hist(sz_mov(:),[0 1 2 4 8 16]);
sz_hist = h/sum(h)


pname = '.';
fname = 'test.mat'

moviedata = uint8((moviedata+1)*128-1);

[fname pname]=uiputfile('*.mat');
save(fullfile(pname,fname),'moviedata','sz_mov','sp_mov','th_mov','degperpix','MovieRate','MovieMag');
