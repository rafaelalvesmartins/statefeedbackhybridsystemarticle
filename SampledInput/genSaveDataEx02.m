function outPut = genSaveDataEx02()

    %% ------------------System's Matrices ------------------------------------
   % Franklin's Book 6th ed page 450 data from 213 page (page 119 did not work well), DC motor
    
    Jm = 1.13e-2;
    b = 0.028;
    Kt = 0.067;
    La= 1e-1;
    Ra = 0.45;
    Ke = 0.067;
    
    ratioJm = 0.00;   
    minRatioJm = 1-ratioJm;
    maxRatioJm = 1+ratioJm;
    
    ratiob = 0.19;   
    minRatiob = 1-ratiob;
    maxRatiob = 1+ratiob;
    
    ratioKt = 0.00; 
    minRatioKt = 1-ratioKt;
    maxRatioKt = 1+ratioKt;
    
    ratioLa = 0.19;   
    minRatioLa = 1-ratioLa;
    maxRatioLa = 1+ratioLa;
    
    ratioRa = 0.00;   
    minRatioRa = 1-ratioRa;
    maxRatioRa = 1+ratioRa;
    
    ratioKe = 0.19;   
    minRatioKe = 1-ratioKe;
    maxRatioKe = 1+ratioKe;
    
    Jm = infsup(Jm*minRatioJm,Jm*maxRatioJm);
    b = infsup(b*minRatiob,b*maxRatiob);
    Kt = infsup(Kt*minRatioKt,Kt*maxRatioKt);
    La= infsup(La*minRatioLa,La*maxRatioLa);
    Ra = infsup(Ra*minRatioRa,Ra*maxRatioRa);
    Ke = infsup(Ke*minRatioKe,Ke*maxRatioKe);
    
%     Matrix A
    A = [0  1   0;
         0  -(b/Jm)     Kt/Jm;
         0  -(Ke/La)    -(Ra/La)];
     
%      Matrix B2
    B2 = [0 0 1/La]';
    
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
    
    
    % Generate polytope manually
    JmVec = [Jm.inf];
    bVec = [b.inf b.sup];
    KtVec = [Kt.inf];
    LaVec = [La.inf La.sup];
    RaVec = [Ra.inf];
    KeVec = [Ke.inf Ke.sup];
    
    
    
%     for i=1:length(JmVec) %JmVec    
%         for j=1:length(bVec) %bVec
%             for k=1:length(KtVec) %KtVec
%                 for l=1:length(LaVec) %LaVec
%                     for m=1:length(RaVec) %RaVec
%                         for n=1:length(KeVec) %KeVec
%                             indexN = n;
%                             indexM = (m-1)*length(KeVec);
%                             indexL = (l-1)*length(KeVec)*length(RaVec);
%                             indexK = (k-1)*length(KeVec)*length(RaVec)*length(LaVec);
%                             indexJ = (j-1)*length(KeVec)*length(RaVec)*length(LaVec)*length(KtVec);
%                             indexI = (i-1)*length(KeVec)*length(RaVec)*length(LaVec)*length(KtVec)*length(bVec);
%                             
%                             sysPoly{indexN+indexM+indexL+indexK+indexJ+indexI}.A = [0   1 0;
%                                                                                 0   -(bVec(j)/JmVec(i))    KtVec(k)/JmVec(i);
%                                                                                 0   -(KeVec(n)/LaVec(l))      -(RaVec(m)/LaVec(l))];
%                            
%                             sysPoly{indexN+indexM+indexL+indexK+indexJ+indexI}.B2 = [0 0 1/LaVec(l)]'; 
%                             sysPoly{indexN+indexM+indexL+indexK+indexJ+indexI}.B1 = [0 0 1/LaVec(l)]'; 
%                             sysPoly{indexN+indexM+indexL+indexK+indexJ+indexI}.C = C.inf; 
%                             sysPoly{indexN+indexM+indexL+indexK+indexJ+indexI}.D2 = D2.inf; 
%                             sysPoly{indexN+indexM+indexL+indexK+indexJ+indexI}.D1 = D1.inf; 
%                              
%                              
%                         end
%                     end
%                 end
%             end
%         end
%     end

    % Initial conditions
    h = 0.1;
    
    %% ------------------Parameters for the Poly Matrices------------------
    numPointsUniSpaced = 20000;
    numPointsBy2Points = 5000;
    numbPointsUniSpacedSub = 25000;
    onlyVertice = 1;
    
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