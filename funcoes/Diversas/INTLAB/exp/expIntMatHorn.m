function [resulExpIntMatHorn] =  expIntMatHorn(A,K)
    
    K2 = generatesK(A);
   
    if(nargin == 1 || (K2>K))
        if(nargin == 1)
            K = K2;
        end
        
        if(K2>K)
            K = K2;
        end
        
    end
    
    n = length(A);

    eyeIntMat =  midrad(eye(n),0);
    
    hornerTil = eyeIntMat;
    for i=K:-1:1   
        hornerTil = eyeIntMat + (hornerTil*(A*1/i));
    end

    residueIntMat = calcResidue(A,K);
    hornerIntMat = hornerTil + residueIntMat;

    resulExpIntMatHorn = hornerIntMat + residueIntMat;

end