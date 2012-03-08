function R = boxcarRate(tsamp, spikeT,w);
dt = zeros(size(spikeT));
for i = 1:length(tsamp)
    dt = spikeT-tsamp(i);   %%% for each time sample, calculate relative time of spikes
    F = abs(dt)<(w/2);   %%% this can be other functions as well
    R(i) = sum(F)/w;
end

