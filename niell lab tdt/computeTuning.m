function f = computeTuning(params,tuning)
f = (params(:,1:5)*tuning(1:5,1)) .* (params(:,6:10)*tuning(6:10,1)) .* (params(:,11:18).*tuning(11:18,1));