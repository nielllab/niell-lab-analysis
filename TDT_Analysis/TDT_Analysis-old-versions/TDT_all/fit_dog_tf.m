function [A1 w1 A2 w2] = fit_dog_tf(y,x);
%%% difference of gaussians for temp freq
    [y0 x0]=max(y);
    p0(1) = 2*y0;
    p0(2) = x(x0);
    if p0(2)==0    %%% low pass
        p0(2)=1;
    end
    p0(3) = p0(1)-y(1);
    p0(4) = p0(2)/2;
    p0
    
    p = nlinfit(x,y,@diff_of_gauss,p0);
    
    p
    
    A1=p(1);
    w1=p(2);
    A2=p(3);
    w2=p(4);
    
    xaxis = x;
    xaxis(1) = x(2)/2;
%     figure
%     plot(log2(xaxis/x(2)),y);
%     hold on
%     plot(log2(xaxis/x(2)),diff_of_gauss(p,x),'g');
    diff_of_gauss(p,x);
    
    
    