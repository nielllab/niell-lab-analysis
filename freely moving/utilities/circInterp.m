function y = circInterp(x0,y0,x);
%%% perform interpolation on circular variables (from -pi to pi)
%%% uses cumulative sum to unwrap values

dy0 = diff(y0);
dy0(dy0>pi)=dy0(dy0>pi)-2*pi;
dy0(dy0<-pi)=dy0(dy0<-pi)+2*pi;
y0unwrap = cumsum([y0(1) dy0]);


yUnwrap = interp1(x0,y0unwrap,x);
y = mod(yUnwrap + pi,2*pi)-pi;

