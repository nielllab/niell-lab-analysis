function [ med s_med ] = distribution( speed )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

figure
hist(x,1:1:40)

med=nanmedian(x)
s_med=semedian(x)
end

