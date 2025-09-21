function [resulIntervInfNorm] = normIntervInf(A)
    [n,m] = size(A.inf);
    
    vecSumLines = zeros(n,1);
    for i=1:n
       for j=1:m
           auxVec = [-A.inf(i,j) A.sup(i,j)];
           vecSumLines(i,1) = vecSumLines(i,1) + max(auxVec);
       end
    end
    
    resulIntervInfNorm = max(vecSumLines);
end