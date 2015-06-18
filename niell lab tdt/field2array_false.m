function a = field2array_false(data,field)
a = arrayfun(@(x)(getfield(x,field)),data,'UniformOutput',false);