function TLMI = gerarLMIAlternadaDerivadaInt(ACal0,ECal0,CCal0,der,X,S,T,U,Alphas,Gammas,Epsilons,mu,nXi,nW,nY,i)
                                        
    %   XaVec
    for j=1:nXi
        if(j==1)
            XaVec = X{i};
        else
            XaVec = [XaVec X{i}];
        end
    end
    XaVec = XaVec';
    
     %   XeVec
    for j=1:nW
        if(j==1)
            XeVec = X{i};
        else
            XeVec = [XeVec X{i}];
        end
    end
    XeVec = XeVec';
    
    %   IVec
    for j=1:nY
        if(j==1)
            IVec = eye(nXi);
        else
            IVec = [IVec eye(nXi)];
        end
    end
    IVec = IVec';
    
        
    U11 = der + ACal0'*X{i}+X{i}*ACal0+S;
    U21 = ECal0'*X{i};
    U22 = -mu*eye(nW)+T;
    U31 = CCal0;
    U32 = zeros(size(U31,1),size(U22,2));
    U33 = -eye(size(U32,1))+U;
    U41 = XaVec;
    U42 = zeros(size(U41,1),size(U32,2));
    U43 = zeros(size(U41,1),size(U33,2));
    U44 = -Alphas;
    U51 = XeVec;
    U52 = zeros(size(U51,1),size(U42,2));
    U53 = zeros(size(U51,1),size(U43,2));
    U54 = zeros(size(U51,1),size(U44,2));
    U55 = -Epsilons;
    U61 = IVec;
    U62 = zeros(size(U61,1),size(U52,2));
    U63 = zeros(size(U61,1),size(U53,2));
    U64 = zeros(size(U61,1),size(U54,2));
    U65 = zeros(size(U61,1),size(U55,2));
    U66 = -Gammas;
    TLMI =  [U11 U21' U31' U41'  U51' U61';
             U21 U22  U32' U42'  U52' U62';
             U31 U32  U33  U43'  U53' U63';
             U41 U42  U43  U44   U54' U64';
             U51 U52  U53  U54   U55  U65';
             U61 U62  U63  U64   U65  U66];
end