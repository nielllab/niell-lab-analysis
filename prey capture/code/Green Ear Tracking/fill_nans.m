function A = fill_nans(A)
% Replaces the nans in each column with 
% previous non-nan values.
for ii = 1:size(A,2)
    I = A(1,ii);
    for jj = 2:size(A,1)
        if isnan(A(jj,ii))
            A(jj,ii) = I;
        else
            I  = A(jj,ii);
        end
    end
end