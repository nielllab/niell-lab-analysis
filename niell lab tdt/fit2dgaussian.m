
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
stafig1 = figure;
colormap(redblue);
imagesc(img)

width=41;
[xi yi] = meshgrid(1:width,1:width) ;
xi = xi- (width+1)/2;
yi = yi - (width+1)/2;


%img = 4 + 10 *exp(-0.5*((xi+4).^2 + (yi-4).^2)/8.^2) +5*rand(size(xi));
figure
imagesc(img)

dx =4;
dsig = 4;
dA=5;
dB = 4;

xguess=-4*dx:dx:4*dx;
yguess=xguess;
sig_guess = 4:dsig:20;
base=min(min(img));
base=0
amp = max(max(img-base));

dtheta = pi/6;

Aguess= [amp-20:dA:amp+20];
Aguess = [Aguess -Aguess]
Bguess =( -12:dB:12) +base;
angle_guess=0:dtheta:pi-dtheta/2;

res = zeros(length(xguess),length(yguess),length(sig_guess),length(Aguess),length(Bguess),length(angle_guess));
res = res(:);
[x y sigx sigy A B angle] = ndgrid(xguess,yguess,sig_guess,sig_guess,Aguess,Bguess,angle_guess);
x = x(:); y = y(:); sigx = sigx(:); sigy = sigy(:); A=A(:); B=B(:);angle=angle(:);

length(x)

matlabpool
tic


    
parfor i = 1:length(x);
xir = (xi-x(i))*cos(angle(i)) + (yi-y(i))*sin(angle(i));
yir =-(xi-x(i))*sin(angle(i)) + (yi-y(i))*cos(angle(i));

z = A(i)*exp(-0.5*((xir.^2)/sigx(i)^2 + (yir.^2)/sigy(i).^2)) + B(i);
    
    res(i) = sum(sum((z-img).^2));
end
toc

[val ind] = min(res);
x0 = x(ind)
y0= y(ind)
sigx0 = sigx(ind)
sigy0 = sigy(ind)
B0 = B(ind)
A0 = A(ind)
angle0 = angle(ind)

xir = (xi-x0)*cos(angle0) + (yi-y0)*sin(angle0);
yir =-(xi-x0)*sin(angle0) + (yi-y0)*cos(angle0);
z_est = A0*exp(-0.5*((xir.^2)/sigx0^2 + (yir.^2)/sigy0.^2)) + B0;
figure
imagesc(z_est);
exp_var = 1- (var(z_est(:)-img(:))/var(img(:)))


%%%%  secondround
xyguess = -dx:dx/3:dx;
sig_guess = -dsig:dsig/3:dsig;
Bguess = -dB:dB/3:dB;
Aguess = -dA:dA/3:dA;
angle_guess = -dtheta:dtheta/3:dtheta;

res = zeros(length(xyguess),length(xyguess),length(sig_guess),length(sig_guess),length(Aguess),length(Bguess),length(angle_guess));
res = res(:);
[x y sigx sigy A B angle] = ndgrid(xyguess+x0,xyguess+y0,sig_guess+sigx0,sig_guess+sigy0,Aguess+A0,Bguess+B0,angle_guess+angle0);
x = x(:); y = y(:); sigx = sigx(:); sigy = sigy(:); A=A(:); B=B(:);angle=angle(:);

length(x)
    
parfor i = 1:length(x);
xir = (xi-x(i))*cos(angle(i)) + (yi-y(i))*sin(angle(i));
yir =-(xi-x(i))*sin(angle(i)) + (yi-y(i))*cos(angle(i));

z = A(i)*exp(-0.5*((xir.^2)/sigx(i)^2 + (yir.^2)/sigy(i).^2)) + B(i);
    
    res(i) = sum(sum((z-img).^2));
end
toc
[val ind] = min(res);
x0 = x(ind)
y0= y(ind)
sigx0 = sigx(ind)
sigy0 = sigy(ind)
B0 = B(ind)
A0 = A(ind)
angle0 = angle(ind)

matlabpool close

figure
subplot(1,2,1)
imagesc(img);
subplot(1,2,2)

xir = (xi-x0)*cos(angle0) + (yi-y0)*sin(angle0);
yir =-(xi-x0)*sin(angle0) + (yi-y0)*cos(angle0);
z_est = A0*exp(-0.5*((xir.^2)/sigx0^2 + (yir.^2)/sigy0.^2)) + B0;
imagesc(z_est);

imagesc(z_est);

exp_var = 1- (var(z_est(:)-img(:))/var(img(:)))

