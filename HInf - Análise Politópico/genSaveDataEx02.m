function saida = genSaveDataEx02()

    %% ------------------Matrizes do sistema Matriz ACal e KCal estavéis---
   % Configura INTLAB
   intvalinit('DisplayInfSup');
   
    % Matrizes do sistema
   g = infsup(-0.4,0.4);
%    g = infsup(0,0);
   
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

            sysPolyCont{i}.C = C.inf;

            sysPolyCont{i}.D2 = D.inf;
            sysPolyCont{i}.D1 = D.inf;
            
            poly.APoly{i} = [sysPolyCont{i}.A sysPolyCont{i}.B2 ; zeros(nu,nx) zeros(nu,nu)];
            poly.EPoly{i} = ECal.inf;
            poly.CPoly{i} = CCal.inf;
                    
        
    end
    
    % Outras variáveis importantes no sistema
    h = 0.1;
    qtdDiv = 20;
    delta = h/qtdDiv;
    
    %% ------------------Parametro para as matrizes Politópicas------------
%     numPointsUniSpaced = 20000;
%     numPointsBy2Points = 5000;
%     numbPointsUniSpacedSub = 25000;
    
    numPointsUniSpaced = 10;
    numPointsBy2Points = 10;
    numbPointsUniSpacedSub = 10;
    
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
    sys.K = K;
     
    sys.ACal = ACal;
    sys.ECal = ECal;
    sys.CCal = CCal;
    sys.KCal = KCal;
    sys.h = h;
    sys.delta = delta;
      
    saida.sys = sys;
    
    % Aux
    saida.aux.tol = tol;
    
    % Simulacao
    sim.numPointsUniSpaced = numPointsUniSpaced;
    sim.numPointsBy2Points = numPointsBy2Points;
    sim.numbPointsUniSpacedSub = numbPointsUniSpacedSub;
    sim.onlyVertice = onlyVertice;
    
    
    saida.poly = poly;
    saida.simSys = sim;

end