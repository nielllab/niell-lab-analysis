function R=naka_rushton(p,x)
halfC=p(1);
b = p(2);

R = 1./(1+(halfC./x).^b);