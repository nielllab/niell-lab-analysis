function obs = rf2fit(coeff, x);
    obs = coeff(1) + coeff(2)*(exp(-0.5*(((x(:,1)-coeff(3))/coeff(4)).^2 + ((x(:,2)-coeff(5))/coeff(6)).^2)));