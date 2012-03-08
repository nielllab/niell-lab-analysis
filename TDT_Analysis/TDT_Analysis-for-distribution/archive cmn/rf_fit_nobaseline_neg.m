function obs = rf_fit_nobaseline_neg( coeff ,x);
%%% rf fit without a baseline (spontaneous rate has been subtracted off)
%%% written by cmn, last modified 02/2008

if coeff(2)<0
    coeff(2)=0;
end

if coeff(2)>max(x)
    coeff(2)=max(x)
end
if coeff(3)>0.5
    coeff(3)=0.5;
end
obs =  coeff(1) * (exp(-0.5*((x-coeff(2))/coeff(3)).^2));