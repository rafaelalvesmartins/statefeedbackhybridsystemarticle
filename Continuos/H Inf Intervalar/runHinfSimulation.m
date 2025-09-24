function [averageCust,stdDeviatCust,maxCust] = runHinfSimulation(A,B2,B1,C,D2,D1,k,nSimulation)
    
    custVector = zeros(nSimulation,1);
    for i=1:nSimulation
        if i == 1
            AIn     = A.inf;
            B2In    = B2.inf;
            B1In    = B1.inf;
            CIn     = C.inf;
            D2In    = D2.inf;
            D1In    = D1.inf;
        elseif i == 2
            AIn     = A.sup;
            B2In    = B2.sup;
            B1In    = B1.sup;
            CIn     = C.sup;
            D2In    = D2.sup;
            D1In    = D1.sup;
        else
            AIn = genRandMatInTheInterv(A);
            B2In = genRandMatInTheInterv(B2);
            B1In = genRandMatInTheInterv(B1);
            CIn = genRandMatInTheInterv(C);
            D2In = genRandMatInTheInterv(D2);
            D1In = genRandMatInTheInterv(D1);   
         end
        sys = ss((AIn+B2In*k),B1In,(CIn+D1In*k),D2In);
        custVector(i) = norm(sys,inf);
    end

    averageCust = mean(custVector);
    stdDeviatCust = std(custVector);
    maxCust = max(custVector);
