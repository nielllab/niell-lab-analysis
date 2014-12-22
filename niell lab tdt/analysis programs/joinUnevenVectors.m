function out=joinUnevenVectors(varargin)
%#Horizontally catenate multiple column vectors by appending zeros 
%#at the ends of the shorter vectors
%#
%#SYNTAX: out = joinUnevenVectors(vec1, vec2, ... , vecN)

    maxLength=max(cellfun(@numel,varargin));
    out=cell2mat(cellfun(@(x)cat(1,x,zeros(maxLength-length(x),1)),varargin,'UniformOutput',false));