function genSaveDataEx12()

   %    LMI relaxations for robust H2 performance analysis of polytopic linear systems
   intvalinit('DisplayInfSup');

   g = infsup(-0.4,0.4);
   
   A = [0    4;
        -1+g -1-g];
    
   B = [0 1]';
   B = infsup(B,B);
   
   C = [1 -2];
   C = infsup(C,C);
    
   D = 0;
   D = infsup(D,D);
    
   E = [1 1]';
   E = infsup(E,E);
       
   K = [-1 -0.75];
%    K = infsup(K,K);
        
    % Matrizes Aumentada
    %     Tamanho das matrizes
    nx = length(A.inf);
    nu = size(B.inf,2);
    nw = size(E.inf,2);
    
    
%     Matrizes caligrafas
    ACal = [A B ; zeros(nu,nx) zeros(nu,nu)];
    ECal = [E ; zeros(nu,nw)];
    CCal = [C D];
    KCal = [eye(nx) zeros(nx,nu);K zeros(nu,nu)];
    
    g = [g.inf g.sup];
    
    for i=1:length(g) 
       

            sysPolyCont{i}.A = [0    4;
                                -1+g(i) -1-g(i)];

            sysPolyCont{i}.B2 = B.inf; 

            sysPolyCont{i}.B1 = E.inf;
clc
            sysPolyCont{i}.C = C.inf;

            sysPolyCont{i}.D2 = D.inf;
            sysPolyCont{i}.D1 = D.inf;
            
            sysCalPoly.APoly{i} = [sysPolyCont{i}.A sysPolyCont{i}.B2 ; zeros(nu,nx) zeros(nu,nu)];
            sysCalPoly.EPoly{i} = ECal.inf;
            sysCalPoly.CPoly{i} = CCal.inf;
                    
        
    end
    
    
    
    % Outras variáveis importantes no sistema
    h = 0.1;
    qtdDiv = 20;
    delta = h/qtdDiv;
       
    
   %% ------------------Parametro para as matrizes Politópicas------------
    numPointsUniSpaced = 500;
    numPointsBy2Points = 200;
    numbPointsUniSpacedSub = 300;
    
%     numPointsUniSpaced = 2;
%     numPointsBy2Points = 2;
%     numbPointsUniSpacedSub = 2;
    
    onlyVertice = 0;
       
    
    %% ------------------Aux Vars------------------------------------------
    tol = 1e-5;

    %% ------------------Sava File-----------------------------------------
    % Folder does not exist     
    if(exist(mfilename) ~= 7)
        mkdir(mfilename);
    end
   
%      % Sys
%     sys.A = A;
%     sys.B = B;
%     sys.E = E;
%     sys.C = C;
%     sys.D = D;
%      
%     sys.ACal = ACal;
%     sys.ECal = ECal;
%     sys.CCal = CCal;
%     sys.KCal = KCal;
%     sys.h = h;
%     sys.delta = delta;
%       
%     saida.sys = sys;
%     
%     % Aux
%     saida.aux.tol = tol;
%     
%     % Simulacao
%     sim.numPointsUniSpaced = numPointsUniSpaced;
%     sim.numPointsBy2Points = numPointsBy2Points;
%     sim.numbPointsUniSpacedSub = numbPointsUniSpacedSub;
%     sim.onlyVertice = onlyVertice;
%     saida.simSys = sim;
    save([mfilename '/data']);

end