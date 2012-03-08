function Y = condenseData(X,f);
%%%% condenses along 1st axis
n = size(X,1)
n_new = floor(n/f)
Y= zeros(n_new,size(X,2));
for i = 1:n_new;
    Y(i,:) = mean(X((i-1)*f+1 : (i-1)*f+f,:),1);
end
