function  gerarLMIAlternadaDerivadaLMILabInt(X,sX,ACal0,CCal0,ECal0,alphas,gammas,epsilons,deltaACal,deltaCCal,deltaECal,nXi,nY,nW,i,mu,c)
                                             
    %   XaVec
    for j=1:nXi
        if(j==1)
            XaVec = sX{i};
        else
            XaVec = [XaVec sX{i}];
        end
    end
    [XaVec,nVar,sXaVec] = lmivar(3,XaVec');


     %   XeVec
    for j=1:nW
        if(j==1)
            XeVec = sX{i};
        else
            XeVec = [XeVec sX{i}];
        end
    end
    [XeVec,nVar,sXeVec] = lmivar(3,XeVec');


   %   IVec
    for j=1:nY
        if(j==1)
            IVec = eye(nXi);
        else
            IVec = [IVec eye(nXi)];
        end
    end
    IVec = IVec';

    % Entrada (1,1)
    lmiterm([c 1 1 X{i}],1,ACal0,'s');
    for k=1:nXi
       for j=1:nXi
            ej = generatesEi(j,nXi);
            lmiterm([c 1 1 alphas(k,j)],abs(deltaACal(k,j))^2*(ej*ej'),1);
       end
    end
   % Entrada (2,1)
   lmiterm([c 2 1 X{i}], ECal0',1);
   % Entrada (2,2)
   lmiterm([c 2 2 mu], -1, 1);
   for k=1:nXi
       for j=1:nW
            fj = generatesEi(j,nW);
            lmiterm([c 2 2 epsilons(k,j)],abs(deltaECal(k,j))^2*(fj*fj'),1);
       end
   end
   % Entrada (3,1)
   lmiterm([c 3 1 0], CCal0);
   % Entrada (3,3)
   lmiterm([c 3 3 0], -1);
   for k=1:nY
       for j=1:nXi
            gk = generatesEi(k,nY); 
            lmiterm([c 3 3 gammas(k,j)],abs(deltaCCal(k,j))^2*(gk*gk'),1);
       end
    end
    % Entrada (4,1)
    lmiterm([c 4 1 XaVec], 1, 1);
    % Entrada (4,4)
    for k=1:nXi
       for j=1:nXi
           ei = generatesEi(((k-1)*nXi)+j,(nXi*nXi)); 
           lmiterm([c 4 4 alphas(k,j)],-1*(ei*ei'),1);
       end
    end
    % Entrada (5,1)
    lmiterm([c 5 1 XeVec], 1, 1);
    % Entrada (5,5)
    for k=1:nXi
       for j=1:nW
           ei = generatesEi(((k-1)*nW)+j,(nXi*nW)); 
           lmiterm([c 5 5 epsilons(k,j)],-1*(ei*ei'),1);
       end
    end
    % Entrada (6,1)
    lmiterm([c 6 1 0], IVec);
    % Entrada (6,6)
    for k=1:nY
       for j=1:nXi
           ei = generatesEi(((k-1)*nXi)+j,(nY*nXi)); 
           lmiterm([c 6 6 gammas(k,j)],-1*(ei*ei'),1);
       end
    end
end