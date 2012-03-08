
close all
clear all
n=0;
for A = 0:0.1:4;
n=n+1;    
x = [0:pi/30:pi*20];
y = A*sin(x)+1;
y= y.*(y>0);

figure
plot(y);
F = abs(fft(y));
figure
plot(F)
OSI(n) = 2*F(11)/F(1);

pref=max(y);
null = min(y);

prefnull(n) = (pref-null)./(pref+null);
end
figure
plot(OSI);
hold on
plot(prefnull,'g');