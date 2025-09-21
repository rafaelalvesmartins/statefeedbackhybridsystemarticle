function [residueIntMat] = calcResidue(A,K)
        
    infNorm = normIntervInf(A);
    rho = (infNorm^(K+1))/(factorial(K+1)*(1-(infNorm/(K+2))));

    n = length(A.inf);
    E = midrad(zeros(n),1);

    residueIntMat = E;
    
    residueIntMat = residueIntMat* rho;
    
end