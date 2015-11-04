function [len] = nanlength(X,dim)

len = sum(~isnan(X),dim);