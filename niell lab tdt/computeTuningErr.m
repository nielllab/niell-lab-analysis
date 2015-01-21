function f = computeTuningErr(tuning,params,data)
%%% fit = tf*tftuning * sf*sftuning * ori*orituning  * (1+gain)*running + baseline
f = (params(:,1:5)*tuning(1:5,1)) .* (params(:,6:10)*tuning(6:10,1))...
    .* (params(:,11:18)*tuning(11:18,1)).* (1 + params(:,19)*tuning(20,1)) +tuning(19);
%%% squared error relative to measured data
f = sum((f-data').^2);