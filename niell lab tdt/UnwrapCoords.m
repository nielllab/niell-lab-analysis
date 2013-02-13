function dVal = UnwrapCoords( uiVal, uiSignBit )

dMaxDiff = hex2dec('4000');   % 16384
dAdd = hex2dec('8000');      % 32768


Vsign = 1 - 2*double(uiSignBit);
idxNeg = find(Vsign<0);

uiVal(idxNeg) = bitcmp(uiVal(idxNeg), 15); % bit-complement neg values
dMod = double(uiVal) .*Vsign;              % restore sign
dDiff = diff(dMod);

idxOverflow = find(abs(dDiff) > dMaxDiff); % find overflows
% correct overflows
dDiff(idxOverflow) = dDiff(idxOverflow) - dAdd * sign(dDiff(idxOverflow));

dVal = cumsum(dDiff);




