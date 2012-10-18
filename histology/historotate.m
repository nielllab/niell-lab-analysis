function [XAlign YAlign]=historotate(X, Y, XC, YC, theta);

% % rotate image so that electrode is vertical
XAlign=(X-XC)*cos(pi/2-theta)-(Y-YC)*sin(pi/2-theta);
YAlign=(X-XC)*sin(pi/2-theta)+(Y-YC)*cos(pi/2-theta);
