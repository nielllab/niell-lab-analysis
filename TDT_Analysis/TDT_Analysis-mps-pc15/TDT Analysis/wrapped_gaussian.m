function y = wrapped_gaussian(p,x);
A1 = p(1);
A2 = p(2);
w = p(3);
B = p(4);
if B<-1.5
    B=-1.5;
end

if A1<0 & A2<0
    A2=-A2;
end
w=abs(w);
if w>pi/2
    w=pi/2;
end
if w<.1
    w=0.1;
end
    
% if A1<0
%     A1=0;
% end
% if A2<0
%     A2=0;
% end
theta0 = x(2,1);
y = A1*exp(-1*dtheta(x(1,:),theta0).^2 / (2*w^2)) +  ...
    A2*exp(-1*dtheta(x(1,:),theta0+pi).^2 /(2*w^2)) + B;



