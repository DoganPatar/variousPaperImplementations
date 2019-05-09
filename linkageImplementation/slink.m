function [pies,lambda] = slink(dataSet)

%Implementation of the single linkage algorithm from
%SLINK: An optimally efficient algorithm for the single-link cluster
%method by R.Sibson

[samples,~] = size(dataSet);
pies = zeros([1,samples]);
lambda = Inf([1,samples]);

% incrementally feeding the data to the SLINK
for s=1:samples
    n = s-1;
    %1st Step
    pies(n+1) = n+1;
    lambda(n+1) = Inf;
    indx=1:1:n;
    %2nd Step
    M = pdist2(dataSet(s,:),dataSet(1:n,:));
    %3rd Step
    for i=1:n
        if lambda(i)>=M(i)
            M(pies(i)) = min(M(pies(i)),lambda(i));
            lambda(i) = M(i);
            pies(i) = n+1;
        else
            M(pies(i)) = min(M(pies(i)),M(i));
        end    
    end
    %4th Step
    temp = lambda(indx) >= lambda(pies(indx));
    pies(temp) = n+1;
end