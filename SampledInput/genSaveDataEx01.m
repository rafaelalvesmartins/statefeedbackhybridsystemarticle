function outPut = genSaveDataEx01()

    %% ------------------System's Matrices ------------------------------------
    % Article: Necessary and Sufficient Conditions for the Controllability and 
    % Observability of a Class of Linear, Time-Invariant Systems with Interval 
    % Plants Wang, Kaining Michel, AN - Example 2

    % Set INTLAB
    intvalinit('DisplayInfSup');

    % Matrix A
    A11 = [0.9, 1.1];
    A12 = [-1.1, -0.9];
    A21 = [-0.05, 0.05];
    A22 = [3.8, 4.2];

    A   = [infsup(A11(1),A11(2)) infsup(A12(1),A12(2));
           infsup(A21(1),A21(2)) infsup(A22(1),A22(2))];

    % Matrix B2
    B2 = [infsup(1,1) infsup(-1,-1)]';
    
    %     Martix B1
    B1 = B2;
    
 % Matrix C
    C = [eye(length(A.inf));zeros(size(B2.inf,2)+size(B1.inf,2),length(A.inf))];
    C = infsup(C,C);
    
    %     Matrix D
    D = [zeros(length(A.inf),size(B2.inf,2)+size(B1.inf,2));eye(size(B2.inf,2)+size(B1.inf,2))];
    
    % Matrix D2
    D2 = D(:,size(B2.inf,2));
    D2 = infsup(D2,D2);
    
    % Matrix D1
    D1 = D(:,size(B2.inf,2)+1:end);
    D1 = infsup(D1,D1);
    
    % Initial conditions
    h = 0.1;
    
    %% ------------------Parameters for the Poly Matrices------------------
    numPointsUniSpaced = 20000;
    numPointsBy2Points = 5000;
    numbPointsUniSpacedSub = 25000;
    
%     numPointsUniSpaced = 10;
%     numPointsBy2Points = 10;
%     numbPointsUniSpacedSub = 10;



  
    onlyVertice = 0;
    
      %% ------------------Aux Vars------------------------------------------
    tol = 1e-5;
    maxNormOfK = 1000;
    precision = 0.01;

  %% ------------------Save File-----------------------------------------
     % Folder does not exist     
    if(exist(mfilename) ~= 7)
        mkdir(mfilename);
    end
    % Sys
    sys.A = A;    
    sys.B2 = B2;   
    sys.B1 = B1;
    sys.C = C;   
    sys.D2 = D2; 
    sys.D1 = D1;
    sys.h = h;
    outPut.sys = sys;
    
  %     Poly
%     outPut.poly = sysPoly;
    
    % Simulation
    sim.numPointsUniSpaced = numPointsUniSpaced;
    sim.numPointsBy2Points = numPointsBy2Points;
    sim.numbPointsUniSpacedSub = numbPointsUniSpacedSub;
    sim.onlyVertice = onlyVertice;
    outPut.simSys = sim;
    
     % Aux
    outPut.aux.tol = tol;
    outPut.aux.precision = precision;

end