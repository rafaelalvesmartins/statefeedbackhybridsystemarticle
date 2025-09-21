function [resulExpIntMatScalingSquaring] =  expIntMatScalingSquaring(A,K,L)
    
    K2 = generatesK(A);
    L2 = generatesL(A,K);
    
    if(nargin == 1 || (K2>K) || (L2>L) )
        if(nargin == 1)
            K = K2;
            L = L2;
        end
        
        if(K2>K)
            K = K2;
        end
        
         if(L2>L)
            L = L2;
         end
    end
  
    newA = A*(1/(2^L));
    resulExpIntMatHorn =  expIntMatHorn(newA,K);
    
    resulExpIntMatScalingSquaring =  resulExpIntMatHorn^(2^L);
   
end