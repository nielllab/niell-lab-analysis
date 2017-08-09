function [params all_fit all_img] = fit2dgabor_angie(wn)




%parfor w = 1:length(wn)
for w = 1:length(wn)
%close all

%%%Load STA data for all cells
% cells
% channel_no = cells(w,1);
% clust_no = cells(w,2);

STA = wn(w).sta;
if ~isempty(STA)
    
%%%Dtermine time point with maximial response
[m ind] = max(abs(STA(:)-127));
[x y t_lag] = ind2sub(size(STA),ind);

STA1 = STA(:,:,t_lag)-128;
% figure
% colormap(redblue)
% imagesc(STA1); axis equal

%%crop image
width = 61;
img = maxsubregion(STA1,width,'abs');
%img = imresize(img,0.5);

stafig1 = figure;
colormap(redblue);
% imagesc(img)


width=size(img,1);
[xi yi] = meshgrid(1:width,1:width) ;
xi = xi- (width+1)/2;
yi = yi - (width+1)/2;

%img = 9*exp(-0.5*((xi+3).^2 + (yi-3).^2)/3.2.^2).*cos(pi/2+0.2 + 2*pi*(xi+3)/13);

%img = 4 + 10 *exp(-0.5*((xi+4).^2 + (yi-4).^2)/8.^2) +5*rand(size(xi));
% figure
% imagesc(img,[-10 10])

base=min(min(img));
base=0
amp = max(max(abs(img)));

dx =width/10;
dsig = width/15;
dA=amp/5;
dB = 4;
dn = 0.1;
dphase=pi/6;


dtheta = pi/6;

nguess=0.001:dn:0.7
phase_guess=0:dphase:pi - dphase/2;
xguess=-4*dx:dx:4*dx;
yguess=xguess;
xsig_guess =dsig:dsig:dsig*5;
ysig_guess=xsig_guess;
Aguess= [amp*0.5:dA:amp*1.5]; %amp =0?
Aguess = [Aguess -Aguess]
%Bguess =( -dB:dB:dB) +base;
Bguess=0;
angle_guess=0:dtheta:pi-dtheta/2;

test_img = img;
test_img(abs(test_img)<0.25*amp)=0;

 [z_est x0 y0 A0 B0 sigx0 sigy0 angle0 n0 phase0] = search2dgabor(xguess,yguess,xsig_guess,ysig_guess,Aguess,Bguess,angle_guess,nguess,phase_guess,xi,yi,test_img)


figure
subplot(1,2,1)
im_min =min(img(:));
im_max=max(img(:));
lim = max(abs(im_min),abs(im_max));
imagesc(img,[-lim lim]);
subplot(1,2,2);
imagesc(z_est,[-lim lim]);
drawnow

exp_var_initial(w) = 1- (var(z_est(:)-img(:))/var(img(:)))


%%%%  secondround
xyguess = -dx:dx/2:dx;
xguess = xyguess+x0;
yguess=xyguess+y0;
sig_guess = -dsig:dsig/2:dsig;
xsig_guess = sig_guess+sigx0;
ysig_guess =sig_guess+sigy0;
%Bguess = (-dB:dB:dB) +B0;
Bguess=B0;
Aguess = (-dA:dA/3:dA)+ A0;
angle_guess = (-dtheta:dtheta/3:dtheta)+angle0 ;
phase_guess = (-dphase:dphase/3:dphase)+phase0;
nguess = (-dn:dn/3:dn)+n0;

 [z_est x0 y0 A0 B0 sigx0 sigy0 angle0 n0 phase0] = search2dgabor(xguess,yguess,xsig_guess,ysig_guess,Aguess,Bguess,angle_guess,nguess,phase_guess,xi,yi,test_img)

figure
subplot(1,2,1)
im_min =min(img(:));
im_max=max(img(:));
lim = max(abs(im_min),abs(im_max));
imagesc(img,[-lim lim]); axis equal
subplot(1,2,2);
imagesc(z_est,[-lim lim]); axis equal

exp_var_final(w) = 1- (var(z_est(:)-img(:))/var(img(:)))

all_img{w}=img;
all_fit{w}=z_est;
all_test_img{w}=test_img;

params(w).x = x0;
params(w).y=y0;
params(w).A=A0;
params(w).B=B0;
params(w).sigx=sigx0;
params(w).sigy=sigy0;
params(w).tilt=angle0;
params(w).nx = n0;
params(w).ny = n0*sigy0/sigx0;
params(w).phase=phase0;
params(w).exp_var=exp_var_final(w);

end

toc
%delete(gcp('nocreate'))

% save(afile,'params','all_img','all_fit','all_test_img','-append');

for w = 1:length(all_img)
  img = all_img{w};
  z_est = all_fit{w};
 if ~isempty(img)
     figure
subplot(2,2,1)
im_min =min(img(:));
im_max=max(img(:));
%lim = max(abs(im_min),abs(im_max));
lim=30;
imagesc(img,[-lim lim]); axis equal
subplot(2,2,2);
imagesc(all_test_img{w},[-lim lim]); axis equal
subplot(2,2,3);
imagesc(z_est,[-lim lim]); axis equal

title(sprintf('Gabor fit'));
% set(gcf, 'PaperPositionMode', 'auto');
% print('-dpsc',psfilename,'-append');

 end
end

end

close all
