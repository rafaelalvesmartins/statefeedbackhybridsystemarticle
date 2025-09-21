function TLMI = gerarSinLMIAlternadaDerivadaPoly(A,B,E,C,D,S,Z,L,X,derS,derL,derZ,mu)
                                                     
    U11 = -derS +X.X1_11*A' + A*X.X1_11';      
    U11 = U11 + X.X1_12*B' + B*X.X1_12';
    
    U21 = -derL' + X.X1_21*A' ;
    U21 = U21 + X.X1_22*B';
    
    U22 = -derZ;
    
    U31 = E' + X.X2_11*A' + X.X2_12*B';
    U32 = zeros(size(U31,1),size(U22,2));
    U33 = -mu*eye(size(U31,1));
   
    U41 = X.X3_11*A' + X.X3_12*B' + C*X.X1_11' + D*X.X1_12';
    U42 = C*X.X1_21' + D*X.X1_22';
    U43 = C*X.X2_11' + D*X.X2_12';
    U44 = -eye(size(U41,1)) + C*X.X3_11' + X.X3_11*C' + D*X.X3_12' + X.X3_12*D';
  
    U51 = S + X.X4_11*A' + X.X4_12*B' -X.X1_11';
    U52 = L - X.X1_21'; 
    U53 =  -X.X2_11';
    U54 = X.X4_11*C'+X.X4_12*D'-X.X3_11';
    U55 = -X.X4_11' - X.X4_11;
   
   
   
    
    U61 = L' + X.X4_21*A' + X.X4_22*B' - X.X1_12';
    U62 = Z -X.X1_22'; 
    U63 = -X.X2_12';
    U64 = X.X4_21*C'+X.X4_22*D'-X.X3_12';
    U65 = -X.X4_12'-X.X4_21;
    U66 = -X.X4_22'-X.X4_22;
    
     
    TLMI =  [U11 U21' U31' U41' U51' U61';
             U21 U22  U32' U42' U52' U62';
             U31 U32  U33  U43' U53' U63';
             U41 U42  U43  U44  U54' U64';
             U51 U52  U53  U54  U55  U65';
             U61 U62  U63  U64  U65  U66];
    
end