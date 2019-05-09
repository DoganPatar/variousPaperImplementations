function [pies,lambda] = clink(dataSet)

%Importance: Defays CLINK algorithm is input VARIANT.  
%Implementation of the complete linkage algorithm from
%An efficient algorithm for a complete link method by D.Defays
[samples,~] = size(dataSet);
pies=ones([1,samples]);
lambda=Inf([1,samples]);


% incrementally feeding the data to the CLINK
for s=2:samples
    
    n = s-1;
    pies(n+1) = n+1;
    lambda(n+1) = Inf;
    
    %newSamples eucledian distance to others
    M = pdist2(dataSet(s,:),dataSet(1:n,:));
    for i=1:n
        if lambda(i)< M(i)
            M(pies(i)) = max(M(pies(i)),M(i));
            M(i) = Inf;
        end
    end
    a = n;
    for i=1:n
        if (lambda(n-i+1) >= M(pies(n-i+1))) && (M(n-i+1) < M(a))
            a = n-i+1;
        end
        if (lambda(n-i+1) < M(pies(n-i+1)))
            M(n-i+1) = Inf;
        end
    end
    
    b = pies(a);
    c = lambda(a);
    pies(a) = n+1;
    lambda(a) = M(a);
    
    rerun = true;
    while rerun
        rerun = false;
        if (a < n) && (b < n)
            d = pies(b);
            e = lambda(b);
            pies(b) = n+1;
            lambda(b) = c;
            b=d;
            c=e;
            rerun = true;
        end
    end
    
    if(a < n) && (b == n)
        pies(b) = n+1;
        lambda(b) = c;
    end

    for i=1:n
        if (pies(pies(i))==n+1) && (lambda(i) >= lambda(pies(i)))
            pies(i) = n+1;
        end
    end 
end



