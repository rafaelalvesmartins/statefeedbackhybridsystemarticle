function saida = genSaveDataEx12NomCont()

%   Optimal H? State Feedback Sampled-data Control
%   Design for Markov Jump Linear Systems
%   G.W. Gabriel, J.C. Geromel & K.M. Grigoriadis
    intvalinit('DisplayInfSup');
    

   m1 = .5;
   m1 = infsup(m1,m1);
   
   m2 = 1;
   m2 = infsup(m2,m2);
   
   b = .2;
   b = infsup(b*0.95,b*1.05);
   
   k1 = 12;
   k1 = infsup(k1,k1);
   
   k2 = 7;
   k2 = infsup(k2*0.98,k2*1.02);
    
  
    A = [0              0       1       0;
         0              0       0       1;
         (-k2-k1)/m1    k2/m1   -b/m1   b/m1;
         k2/m2          -k2/m2  b/m2    -b/m2];
   

    B = [0 0 0 1/m2]';

    
    C = [0 10 0 0;
         0  0 0 1;
         0  0 0 0]; 
    C = infsup(C,C);
    
    D = [0 0 1]';
    D = infsup(D,D);
    
    E = [0 0 1/m1 0]';

    % Sistemas para comparação conforme solicitação do revisor
    % 1) Sistema nominal (valores centrais) para sampled-data
    A_nominal = [0              0       1       0;
                 0              0       0       1;
                 (-k2.mid-k1.mid)/m1.mid    k2.mid/m1.mid   -b.mid/m1.mid   b.mid/m1.mid;
                 k2.mid/m2.mid          -k2.mid/m2.mid  b.mid/m2.mid    -b.mid/m2.mid];

    B_nominal = [0 0 0 1/m2.mid]';
    E_nominal = [0 0 1/m1.mid 0]';
    C_nominal = [0 10 0 0; 0 0 0 1; 0 0 0 0];
    D_nominal = [0 0 1]';

    % 2) Sistema contínuo robusto (para posterior discretização)
    sys_continuous_robust.A = A;
    sys_continuous_robust.B = B;  
    sys_continuous_robust.E = E;
    sys_continuous_robust.C = C;
    sys_continuous_robust.D = D;
    
    % Tamanho das matrizes
    nx = length(A.inf);
    nu = size(B.inf,2);
    nw = size(E.inf,2);
    
    % Matrizes caligrafas
    ACal = [A B ; zeros(nu,nx) zeros(nu,nu)];
    ECal = [E ; zeros(nu,nw)];
    CCal = [C D];
    
    
    % Generate polytope manually
   
   b = [b.inf b.sup];
   
   
   k2 = [k2.inf k2.sup];
    
   
    for i=1:length(b) 
        for j=1:length(k2) 
                              indexJ = j;
                              indexI = (i-1)*length(k2);
                            
                                sysPolyCont{indexJ+indexI}.A = [0              0       1       0;
                                                                 0              0       0       1;
                                                                 (-k2(j)-k1.sup)/m1.sup    k2(j)/m1.sup   -b(i)/m1.sup   b(i)/m1.sup;
                                                                 k2(j)/m2.sup          -k2(j)/m2.sup  b(i)/m2.sup    -b(i)/m2.sup];


                            sysPolyCont{indexJ+indexI}.B2 = [0 0 0 1/m2.sup]';
                            
                            sysPolyCont{indexJ+indexI}.B1 = [0 0 1/m1.sup 0]';

                            sysPolyCont{indexJ+indexI}.C = [0 10 0 0;
                                                             0  0 0 1;
                                                             0  0 0 0]; 
                               

                            sysPolyCont{indexJ+indexI}.D2 = [0 0 1]';
                            sysPolyCont{indexJ+indexI}.D1 = [0 0 0]';
        end
    end


    % Outras variáveis importantes no sistema
    h = 0.02;
    qtdDiv = 5;
    delta = h/qtdDiv;
    
    %% ------------------Parametro para as matrizes Politópicas------------
     numPointsUniSpaced = 5;
     numPointsBy2Points = 2;
     numbPointsUniSpacedSub = 2;
    
    onlyVertice = 1;
       
    
    %% ------------------Aux Vars------------------------------------------
    tol = 1e-5;

    %% ------------------Sava File-----------------------------------------
    % Folder does not exist     
    if(exist(mfilename) ~= 7)
        mkdir(mfilename);
    end
   
     % Sys
    sys.A = A;
    sys.B = B;
    sys.E = E;
    sys.C = C;
    sys.D = D;
     
    sys.ACal = ACal;
    sys.ECal = ECal;
    sys.CCal = CCal;
    sys.h = h;
    sys.delta = delta;
    sys.sysPolyCont = sysPolyCont;
    
    % Sistemas para comparação
    sys.A_nominal = A_nominal;
    sys.B_nominal = B_nominal;
    sys.E_nominal = E_nominal;
    sys.C_nominal = C_nominal;
    sys.D_nominal = D_nominal;
    sys.sys_continuous_robust = sys_continuous_robust;
      
    saida.sys = sys;
    
    % Aux
    saida.aux.tol = tol;
    
    % Simulacao
    sim.numPointsUniSpaced = numPointsUniSpaced;
    sim.numPointsBy2Points = numPointsBy2Points;
    sim.numbPointsUniSpacedSub = numbPointsUniSpacedSub;
    sim.onlyVertice = onlyVertice;
    saida.simSys = sim;

end