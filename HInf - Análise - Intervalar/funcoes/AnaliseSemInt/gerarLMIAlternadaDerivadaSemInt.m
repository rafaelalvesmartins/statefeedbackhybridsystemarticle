function  TLMI = gerarLMIAlternadaDerivadaSemInt(ACal,ECal,CCal,der,X,mu,nW,i)
  
    U11 = der + ACal'*X{i}+X{i}*ACal;
    U21 = ECal'*X{i};
    U22 = -mu*eye(nW);
    U31 = CCal;
    U32 = zeros(size(U31,1),size(U22,2));
    U33 = -eye(size(U32,1));
    
    TLMI =  [U11 U21' U31';
             U21 U22  U32';
             U31 U32  U33];
end