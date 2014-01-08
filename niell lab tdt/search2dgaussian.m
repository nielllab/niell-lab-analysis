function [z_est x0 y0 A0 B0 sigx0 sigy0 angle0] = search2dgaussian(xguess,yguess,xsig_guess,ysig_guess,Aguess,Bguess,angle_guess,xi,yi,img)


res = zeros(length(xguess),length(yguess),length(xsig_guess),length(ysig_guess),length(Aguess),length(Bguess),length(angle_guess));
res = res(:);
[x y sigx sigy A B angle] = ndgrid(xguess,yguess,xsig_guess,ysig_guess,Aguess,Bguess,angle_guess);
x = x(:); y = y(:); sigx = sigx(:); sigy = sigy(:); A=A(:); B=B(:);angle=angle(:);

length(x)

    tic    
for i = 1:length(x);
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