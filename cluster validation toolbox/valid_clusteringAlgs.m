function [label,NC,Dist,dmax,Cmerge] = valid_clusteringAlgs(data,alg,k,dtype,Dist,dmax,kp)
% Running different clustering algorithms

k2 = k(2);
k = k(1);
[nrow, dim] = size(data);
Cmerge = [];
NC = [];

% Pearson correlation coefficients [-1,1] is normalized to Pearson distance [0,1]
if (size(Dist,1) ~= nrow) && (alg ~= 8)
    if dtype == 1
        [Dist, dmax] = similarity_euclid(data);
    else
        Dist = 1-(1+similarity_pearson(data'))/2;
        dmax = 1;
    end
end

if alg == 1                       % PAM from LIBRA
    dissim = [];
    for i = 2:nrow
        dissim = [dissim Dist(i,1:i-1)];         % dissimilarity vector
    end
    Scluster = pam(dissim', k);  %Scluster = pam(data, k, 4*ones(1,dim));
    label = double(Scluster.ncluv)';

elseif alg == 2               %  k-means from Mathworks
    R = 'sqEuclidean';
    if dtype == 2
        R = 'correlation';
    end%'Start', 'uniform'
    label = kmeans(data, k, 'distance', R, 'replicates', kp,'emptyaction','drop');

elseif alg > 2 && alg < 6    % agglomerative hierarchical from David Corney
    Dist = Dist+diag(ones(1,nrow).*Inf);
    if alg == 3
        [label, Cmerge] = agg_hierarchical(Dist, 'single', [k k2]);
    elseif alg == 4
        [label, Cmerge] = agg_hierarchical(Dist, 'complete', [k k2]);
    else
        [label, Cmerge] = agg_hierarchical(Dist, 'centroid', [k k2]);
    end

elseif alg == 6
    neurons = som_netlab(data, k);   % SOM from Netlab Toolbox
    %neurons = neural_gas(data,k,4*dtype,1.0,k/2,dtype);  from SOM Toolbox
    R = ones(k, nrow);
    for j = 1:k                      % calculating distances from neurons
        if dtype == 1
            D = data-ones(nrow,1)*neurons(j,:);
            D = (D.^2)';
            R(j,:) = sum(D);
        else
            R(j,:) = similarity_pearsonC(data', neurons(j,:)');
        end
    end
    [D, label] = min(R);       % those best match neruons
    label = label';

elseif alg == 7
    valid_apcluster;

elseif alg == 8
    % new clustering algorithm here
end