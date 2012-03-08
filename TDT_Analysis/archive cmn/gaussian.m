function y = gaussian(p,x);
    A = p(1);
    x0 = p(2);
    w = p(3);
    b = p(4);
   if x0<1
       x0=1;
   end
   if A<0
       A=0;
   end
    if b<-1.5
        b=-1.5;
    end
    y = b + A*exp(-1*((x-x0).^2)./(2*w^2));
    