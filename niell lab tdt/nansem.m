function [sem] = nansem(x,dim)

if nargin<2
    dim = 1;
end
sem = nanstd(x,[],dim)./sqrt(nanlength(x,dim));