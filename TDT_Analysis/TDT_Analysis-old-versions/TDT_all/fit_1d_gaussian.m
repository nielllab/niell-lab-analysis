function [A  x0 w B yfit] = fit_1d_gaussian(y);
    [p0(1) p0(2)] = max(y);
    p0(3)=4;
    p0(4) = min(y);
    y_interp = interp1(1:size(y,2),y,1:.1:size(y,2),'spline');
    p = nlinfit(1:.1:size(y,2),y_interp,@gaussian,p0);
    
    A= p(1);
    x0= p(2);
    w = p(3);
    B = p(4);
   
    yfit = gaussian(p,1:size(y,2));
    