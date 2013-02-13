function [dX dY] = GetMouseTraj( dXY )

% have to do this stupid trick to keep all bits when converting to uint
% otherwise all negatives are converted to 0
idxNeg= find(dXY<0);
Ysignbit = zeros(size(dXY));
Ysignbit(idxNeg) = 1;

dXY(idxNeg) = dXY(idxNeg) +  hex2dec('80000000'); % reset sign bit
uiXY = uint32(dXY);  % now we CAN convert

uiMask15 = uint32(hex2dec('7FFF')); % keep low 15 bit

uiX = bitand(uiXY, uiMask15);
uiXsignBit = bitand(bitshift(uiXY, -15), 1 ); 
dX = UnwrapCoords(uiX, uiXsignBit);

uiY = bitand(bitshift(uiXY, -16), uiMask15);
uiYsignBit = uint32(Ysignbit); 
dY = UnwrapCoords(uiY, uiYsignBit);

% add the starting point, always 0
dX = [0 dX]; 
dY = [0 dY];

