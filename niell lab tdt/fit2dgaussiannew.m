
close all

%%%Load STA data for all cells

STA = wn.sta;

%%%Dtermine time point with maximial response
[m ind] = max(abs(STA(:)-127));
[x y t_lag] = ind2sub(size(STA),ind)

STA1 = STA(:,:,t_lag)-128;
figure
colormap(redblue)
imagesc(STA1)

%%%crop image
width = 41;
img = maxsubregion(STA1,width,'neg');
img = imresize(img,0.5);

stafig1 = figure;
colormap(redblue);
imagesc(img)


width=size(img,1);
[xi yi] = meshgrid(1:width,1:width) ;
xi = xi- (width+1)/2;
yi = yi - (width+1)/2;


%img = 4 + 10 *exp(-0.5*((xi+4).^2 + (yi-4).^2)/8.^2) +5*rand(size(xi));
figure
imagesc(img)

dx =width/10;
dsig = width/10;
dA=amp/5;
dB = 4;


base=min(min(img));
base=0
amp = max(max(img-base));

dtheta = pi/6;

xguess=-4*dx:dx:4*dx;
yguess=xguess;
xsig_guess =dsig:dsig:dsig*5;
ysig_guess=xsig_guess;
Aguess= [amp*0.5:dA:amp*1.5];
Aguess = [Aguess -Aguess]
Bguess =( -12:dB:12) +base;
angle_guess=0:dtheta:pi-dtheta/2;


 [z_est x0 y0 A0 B0 sigx0 sigy0 angle0] = search2dgaussian(xguess,yguess,xsig_guess,ysig_guess,Aguess,Bguess,angle_guess,xi,yi,img)

figure
subplot(1,2,1)
imagesc(img);
subplot(1,2,2);
imagesc(z_est);

exp_var = 1- (var(z_est(:)-img(:))/var(img(:)))


%%%%  secondround
xyguess = -dx:dx/3:dx;
xguess = xyguess+x0;
yguess=xyguess+y0;
sig_guess = -dsig:dsig/3:dsig;
xsig_guess = sig_guess+sigx0;
ysig_guess =sig_guess+sigy0;
Bguess = (-dB:dB/3:dB) +B0;
Aguess = (-dA:dA/3:dA)+ A0;
angle_guess = (-dtheta:dtheta/3:dtheta)+angle0 ;

 [z_est x0 y0 A0 B0 sigx0 sigy0 angle0] = search2dgaussian(xguess,yguess,xsig_guess,ysig_guess,Aguess,Bguess,angle_guess,xi,yi,img)

figure
subplot(1,2,1)
imagesc(img);
subplot(1,2,2);
imagesc(z_est);

exp_var = 1- (var(z_est(:)-img(:))/var(img(:)))

