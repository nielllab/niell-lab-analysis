function sf_values = sf_from_dog(params);

A1=abs(params(1));
A2=abs(params(3));
w1= abs(params(2));
w2=abs(params(4));
if w1>100;
    w1=100;
end
if w2>200
    w2=200;
end
[A1 w1 A2 w2]
x=-100:100;
y = A1*exp(-1*(x.^2)/(w1^2)) - A2*exp(-1*(x.^2)/(w2^2));
figure
plot(y);
sf = abs(fft(y));
figure
plot(sf(1:50));
sf_values = [sf(1) sf(2) sf(3) sf(5) sf(9) sf(17) sf(33)];