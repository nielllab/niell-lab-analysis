function out = spikeFFT(t,f)
out = sum(exp(2*pi*sqrt(-1)*f*t));