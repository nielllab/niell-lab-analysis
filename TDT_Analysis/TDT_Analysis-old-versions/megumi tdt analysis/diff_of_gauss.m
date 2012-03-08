function y= diff_of_gauss(params, x);
A1 = abs(params(1));
A2 =abs(params(3));
w1 = abs(params(2));
w2= abs(params(4));
if A1<0
    A1=0;
end
if A2<0
    A2=0;
end
% if A1>200;
%     A1=200;
% end
% if A2>200
%     A2=200;
% end
% if w1>1;
%     w1 =1;
% end
% if w2>0.5;
%     w2=0.5;
% end
% if w2>w1;
%     w2=w1;
% end
 
y = abs(A1)*exp(-(x.^2/w1.^2)) - abs(A2)*exp(-(x.^2/w2.^2));
y(y<0)=y(y<0)/5;