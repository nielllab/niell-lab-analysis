function y = 1dgaussian(p,x);
    A = p(1);
    x0 = p(2);
    w = p(3);
    b = p(4);
    if b<0
        b=0;
    end
    y = b + A*exp(-1*((x-x0).^2)./(2*w^2));
    