function y= dog_peak(params);
A1 = abs(params(1));
A2 = abs(params(3));
w12 = params(2)^2;
w22 = params(4)^2;

y= sqrt( ((w12 * w22)/( w12 - w22)) * log((A2*w12)/(A1*w22)));