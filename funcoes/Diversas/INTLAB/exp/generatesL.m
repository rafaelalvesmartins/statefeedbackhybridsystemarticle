function L = generatesL(A,K)
    
    
    L = abs(log2(normIntervInf(A)/(K+2)));
    L = ceil(ceil((L + 1))*(25/16));
    
%     Article says L>log2(norm(A)/(k+2))
    
end