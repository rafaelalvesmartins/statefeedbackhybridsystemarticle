function outPut = synHInfKContIntWOD1(A,B2,B1,C,D,param)
% Calculates K minimizing HInf norm, feed back control, from interval matrices: A,B2,B1,C and D
% (continuos)

    if nargin == 6
        if isfield(param,'tol')
            tol = param.tol;
        else
            tol = 1e-7;
        end    
    else
        tol = 1e-7;
    end

%    ---------------------Usefull matrices----------------------------
    A0 = mid(A);
    deltaA = rad(A);
    
    B20 = mid(B2);
    deltaB2 = rad(B2);

    B10 = mid(B1);
    deltaB1 = rad(B1);

    C0 = mid(C);
    deltaC = rad(C);

    D0 = mid(D);
    deltaD = rad(D);

%    ---------------------Dimension of input----------------------------
    nx = length(A.inf);
    nu = size(B2.inf,2);
    nw = size(B1.inf,2);
    ny = size(C.inf,1);
    
%    ---------------------Variable Declariation----------------------------
    outPut.cpusec_m = clock;
    
    %LMI rows counter
    outPut.line = 0;
    
    %New LMI system
    LMIs = set([]);
    
    X = sdpvar(nx,nx,'symmetric');
    LMIs = LMIs + (X >= 0);
    outPut.line = outPut.line + nx;

    L = sdpvar(nu,nx,'full');
    
    %   Set obj
    mu = sdpvar(1,1);
    obj = mu;
       
%   Alphas 
    for i=1:nx
        for j=1:nx
            alphas(i,j) = sdpvar(1,1); % Scalar
        end
    end
    
%    Betas 2
    for i=1:nx
        for j=1:nu
            betas2(i,j) = sdpvar(1,1); % Scalar
        end
    end
    
%    Betas 1
    for i=1:nx
        for j=1:nw
            betas1(i,j) = sdpvar(1,1); % Scalar
        end
    end
    
%    Gammas
    for i=1:ny
        for j=1:nx
            gammas(i,j) = sdpvar(1,1); % Scalar
        end
    end
    
%    Epsilons
    for i=1:ny
        for j=1:nu
            epsilons(i,j) = sdpvar(1,1); % Scalar
        end
    end
   
%    ---------------------Specific Variable Declaration-------------------

%   Calculates S
%   First S' sum
    for i=1:nx
       for j=1:nx       
            ei = generatesEi(i,nx); 
            productOfDeltaAAlphaij = (abs(deltaA(i,j))^2)*alphas(i,j)*ei*ei';
            
%             First time i = j = 1
            if (i==1 && j == 1) 
                firstSum = productOfDeltaAAlphaij;
            else
                firstSum = productOfDeltaAAlphaij + firstSum;
            end
            
       end
    end
    
%   Second S' sum
   for i=1:nx
       for j=1:nu
            ei = generatesEi(i,nx); 
            productOfDeltaBBetaij = (abs(deltaB2(i,j))^2)*betas2(i,j)*ei*ei';
            
%             First time i = j = 1
            if (i==1 && j == 1) 
                secondSum = productOfDeltaBBetaij;
            else
                secondSum = productOfDeltaBBetaij + secondSum;
            end
            
       end
   end
   
   S = firstSum + secondSum;
    
   
   
%   Calculates T
%   First T' sum
    for i=1:ny
       for j=1:nx       
            gi = generatesEi(i,ny); 
            productOfDeltaCGammaij = (abs(deltaC(i,j))^2)*gammas(i,j)*gi*gi';
            
%             First time i = j = 1
            if (i==1 && j == 1) 
                firstSum = productOfDeltaCGammaij;
            else
                firstSum = productOfDeltaCGammaij + firstSum;
            end
            
       end
    end
    
%  Second T' sum
   for i=1:ny
       for j=1:nu
            
            gi = generatesEi(i,ny); 
            productOfDeltaDEpsilonij = (abs(deltaD(i,j))^2)*epsilons(i,j)*gi*gi';
            
%             First time i = j = 1
            if (i==1 && j == 1) 
                secondSum = productOfDeltaDEpsilonij;
            else
                secondSum = productOfDeltaDEpsilonij + secondSum;
            end
            
       end
   end
   
   T = firstSum + secondSum;
   
%   Calculates U
   for i=1:nx
       for j=1:nw
            hj = generatesEi(i,nw); 
            productOfDeltaBBetaij = (abs(deltaB1(i,j))^2)*betas1(i,j)*hj*hj';
            
%             First time i = j = 1
            if (i==1 && j == 1) 
                U = productOfDeltaBBetaij;
            else
                U = productOfDeltaBBetaij + U;
            end
            
       end
   end
   
   
%   Xa, Lb2 and I1
    for i=1:nx
        if(i==1)
            Xa = X;
            Lb2 = L';
        else
            Xa = [Xa X];
            Lb2 = [Lb2 L']; 
        end
    end
    Xa = Xa';
    Lb2 = Lb2';
    
    %   Xc and Ld
    for i=1:ny
        if(i==1)
            Xc = X;
            Ld = L'; 
        else
            Xc = [Xc X];
            Ld = [Ld L'];
        end
    end
    Xc = Xc';
    Ld = Ld';
    
    %  I1
    for i=1:nw
        if(i==1)
            Ib1 = eye(nx);
        else
            Ib1 = [Ib1 eye(nx)];
        end
    end
    Ib1 = Ib1';
    
    Alphas = diag(reshape(alphas',[nx*nx, 1]));
    Gammas = diag(reshape(gammas',[nx*ny, 1]));
    Betas2 = diag(reshape(betas2',[nu*nx, 1]));
    Betas1 = diag(reshape(betas1',[nw*nx, 1]));
    Epsilons = diag(reshape(epsilons',[nu*ny, 1]));
    
    V = A0*X+B20*L;
    U11 = V+V'+S;
    U21 = C0*X+D0*L;
    U31 = B10';
    U41 = Xa;
    U51 = Xc;
    U61 = Lb2;
    U71 = Ld;
    U81 = Ib1;
    
    U22 = -eye(ny)+T;
    U32 = zeros(size(U31,1),size(U22,2));
    U42 = zeros(size(U41,1),size(U22,2));
    U52 = zeros(size(U51,1),size(U22,2));
    U62 = zeros(size(U61,1),size(U22,2));
    U72 = zeros(size(U71,1),size(U22,2));
    U82 = zeros(size(U81,1),size(U22,2));
    
    U33 = -mu*eye(nw)+U;
    U43 = zeros(size(U41,1),size(U33,2));
    U53 = zeros(size(U51,1),size(U33,2));
    U63 = zeros(size(U61,1),size(U33,2));
    U73 = zeros(size(U71,1),size(U33,2));
    U83 = zeros(size(U81,1),size(U33,2));
    
    U44 = -Alphas;
    U54 = zeros(size(U51,1),size(U44,2));
    U64 = zeros(size(U61,1),size(U44,2));
    U74 = zeros(size(U71,1),size(U44,2));
    U84 = zeros(size(U81,1),size(U44,2));
    
    U55 = -Gammas;
    U65 = zeros(size(U61,1),size(U55,2));
    U75 = zeros(size(U71,1),size(U55,2));
    U85 = zeros(size(U81,1),size(U55,2));
    
    U66 = -Betas2;
    U66 = zeros(size(U61,1),size(U66,2));
    U76 = zeros(size(U71,1),size(U66,2));
    U86 = zeros(size(U81,1),size(U66,2));
    
    U77 = -Epsilons;
    U87 = zeros(size(U81,1),size(U77,2));
    
    U88 = -Betas1;
    
    T = [U11 U21' U31' U41' U51' U61' U71' U81';
         U21 U22  U32' U42' U52' U62' U72' U82';
         U31 U32  U33  U43' U53' U63' U73' U83';
         U41 U42  U43  U44  U54' U64' U74' U84';
         U51 U52  U53  U54  U55  U65' U75' U85';
         U61 U62  U63  U64  U65  U66  U76' U86';
         U71 U72  U73  U74  U75  U76  U77  U87';
         U81 U82  U83  U84  U85  U86  U87  U88];
     
    LMIs = LMIs + (T <= 0);
    outPut.line = outPut.line + size(T,1);
    
    outPut.cpusec_m = etime(clock,outPut.cpusec_m);
    outPut.var = size(getvariables(LMIs),2);
    
     sol = optimize(LMIs,obj,sdpsettings('shift','1','verbose',0,'solver','sedumi','sedumi.eps',1e-20,'sedumi.maxiter',2000,'sedumi.numtol',1e-20));
    outPut.cpusec = sol.solvertime;
    
    p = min(checkset(LMIs));
    outPut.delta = p;
    
    outPut.feas = 0;
    if(p > -tol ) % Is feasible
        outPut.feas = 1;
        outPut.X = double(X);
        outPut.L = double(L);
        outPut.K = double(L)/(outPut.X); 
        outPut.mu = double(mu);
    end
end