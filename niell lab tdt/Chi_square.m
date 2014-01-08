function [ chi2 pval df ] = Chi_square(Cat,Group1,Group2)

%Perform Chi-square test when you know: number of categories (Cat), count
%in each category (Group1 and Group2) or can just pass in the list of
%values from two groups (Group1 and Group2) 
 
if length(Group1)<= Cat; %determine if data are counts in each category for each test group
    

M = sum(Group1);%add up all the counts to get total number of observed events
N = sum(Group2);



phat = (Group1+Group2) ./ (M+N); % calculates probabilty of each event
eGroup1 = phat*M; eGroup2 = phat*N;% caculates expected probability based on contingency table principles
chi2 = sum(([Group1 Group2] - [eGroup1 eGroup2]).^2 ./ [eGroup1 eGroup2]);% caclulates Chi-square statistic
df = Cat-1; %calculates degress of freedom (#rows-1)*(#colum-1) row and column in the contingency table, ie. rows=groups and columns = categories
pval = 1 - chi2cdf(chi2,df);

%%to calculate effect size (phi for 2 categories and cV for more than two
%%categories
% phi =  sqrt(chi2/(M+N));
% cV=  sqrt(chi2/((M+N)*(k-1))); %% k will be the number of the smaller number of either rows or col, Cramér's V  is computed by taking the square root of the chi-squared statistic divided by the sample size and the length of the minimum dimension (k is the smaller of the number of rows r or columns c).

%using crosstab method to calc Chi-square from lists of frequency data in
%multiple categories
  
    
    [chi2,pval] = crosstab([Group1 Group2],[ones(size(Group1)) 2*ones(size(Group2))])
    
    %calc by hand got from http://www.mathworks.com/matlabcentral/newsreader/view_thread/290714
        % p = rand(1,k); p = p./sum(p);
    % M = 200; N = 250;
    % x = randsample(1:k,M,true,p); m = histc(x,1:k);
    % y = randsample(1:k,N,true,p); n = histc(y,1:k);
    % 
    % phat = (m+n) ./ (M+N);
    % em = phat*M; en = phat*N;
    % chi2 = sum(([m n] - [em en]).^2 ./ [em en]);
    % df = k-1;
    % pval = 1 - chi2cdf(chi2,df);
    
    
end 



end

 
   
   

