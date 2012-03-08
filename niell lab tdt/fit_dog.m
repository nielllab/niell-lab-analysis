function [A1 w1 A2 w2] = fit_dog(y,x);
    [y0 x0]=max(y);
    p0(1) = 2*y0;
    p0(2) = x(x0);
    if p0(2)==0    %%% low pass
        p0(2)=.02;
    end
    p0(3) = p0(1)-y(1);
    p0(4) = p0(2)/4;
    p0
    
    y
    y(y<-1)=0;
    p = nlinfit(x,y,@diff_of_gauss,p0);
    
    p
    
    A1=abs(p(1))
    w1=p(2)
    A2=abs(p(3))
    w2=p(4)
    
    xaxis = x;
    xaxis(1) = x(2)/2
%     figure
%     plot(log2(xaxis/x(2)),y);
%     hold on
%     plot(log2(xaxis/x(2)),y,'bo');
%     plot(log2(xaxis/x(2)),diff_of_gauss(p,x),'g--');
    
    
    