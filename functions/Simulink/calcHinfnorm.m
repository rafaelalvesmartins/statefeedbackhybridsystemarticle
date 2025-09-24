function [hInfNorm] = calcHinfnorm(A,B2,B1,C,D2,D1,K,param)

  ACl = A+B2*K;
  BCl = B1;
  CCl = C+D2*K;
  DCl = D1;
   
   
    if(nargin == 8 && isfield(param,'disc'))
         if sum(mag(eig(ACl))>1)
             error('ERROR - Eigen values mag greater than 1 :/');
         end
         sysClMidPoly = ss(ACl,BCl,CCl,DCl,-1);
    else
         if sum(eig(ACl)>0)
             error('ERROR - Eigen value on the RSP :/');
         end
        sysClMidPoly = ss(ACl,BCl,CCl,DCl);
    end
    hInfNorm = hinfnorm(sysClMidPoly);
end