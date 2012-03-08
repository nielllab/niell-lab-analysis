function y = solve_half_dog(w,params);
   
    y = diff_of_gauss(params(1:4),w)-params(5)/2;