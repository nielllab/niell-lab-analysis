data = handles.data;
distype = get(handles.popMeasure, 'Value');
R = 'euclid';
if distype == 2
   R = 'correlation';
end
Y= pdist(data, R);
ialg = get(handles.popAlgorithm, 'Value');
R = 'average';
if ialg == 3
   R = 'single';
elseif ialg == 4
   R = 'complete';
elseif ialg == 5
   R = 'centroid';
end
Z= linkage(Y, R);     %T = cluster(Z,'cutoff', 1.8); 
figure('color','white');
[H,T] = dendrogram(Z,'colorthreshold','default');
box on;
