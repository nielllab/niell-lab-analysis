function neurons = som_netlab(data, k)

[nrow, dim] = size(data);
ns = k/2;
nsize = ceil(ns);
if nsize > ns
  nsize = [1 k];
else
  nsize = [2 nsize];
end
net = som(dim, nsize);
options = zeros(1,18);

% Ordering phase
options(1) = 0;               % printing iterations
options(2) = 0.0001;
options(3) = 0.0001;
options(4) = 0.000001;
options(14) = 50;           % iterations
options(18) = 0.9;          % Initial learning rate
options(16) = 0.05;        % Final learning rate
options(17) = 0.5*nsize(2);  % Initial neighbourhood size
options(15) = 0.01;        % Final neighbourhood size
net = somtrain(net, options, data);

options(14) = 100+20*k;
options(18) = 0.05;
options(16) = 0.001;
options(17) = 0;
options(15) = 0;
if 1
net = somtrain(net, options, data);      % Convergence phase
else
options(6) = 1;
net = somtrain(net, options, data);      % batch training
end

Re = net.map;
neurons = zeros(dim, k); 
for ni = 1:nsize(2)
    for nj = 1:nsize(1)
      ns = (ni-1)*nsize(1)+nj;
      neurons(:, ns) = Re(:, nj, ni); 
    end
end
neurons = neurons'; 