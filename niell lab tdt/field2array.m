function a = field2array(data,field)
a = arrayfun(@(x)(getfield(x,field)),data,'UniformOutput',false);