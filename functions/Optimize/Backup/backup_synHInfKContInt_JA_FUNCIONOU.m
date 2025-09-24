function outPut = synHInfKContInt(A,B2,B1,C,D2,D1,param)
% Calculates K minimizing HInf norm, feed back control, from interval
% matrices: A,B2,B1,C, D2 and D1
% (continuos)

    if nargin == 7
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

    D20 = mid(D2);
    deltaD2 = rad(D2);
    
    D10 = mid(D1);
    deltaD1 = rad(D1);
    
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
    
%    Epsilons 2
    for i=1:ny
        for j=1:nu
            epsilons2(i,j) = sdpvar(1,1); % Scalar
        end
    end
    
%    Epsilons 1
    for i=1:ny
        for j=1:nw
            epsilons1(i,j) = sdpvar(1,1); % Scalar
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
    
%   Second S sum
   for i=1:nx
       for j=1:nu
            ei = generatesEi(i,nx); 
            productOfDeltaB2Beta2ij = (abs(deltaB2(i,j))^2)*betas2(i,j)*ei*ei';
            
%             First time i = j = 1
            if (i==1 && j == 1) 
                secondSum = productOfDeltaB2Beta2ij;
            else
                secondSum = productOfDeltaB2Beta2ij + secondSum;
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
            productOfDeltaD2Epsilon2ij = (abs(deltaD2(i,j))^2)*epsilons2(i,j)*gi*gi';
            
%             First time i = j = 1
            if (i==1 && j == 1) 
                secondSum = productOfDeltaD2Epsilon2ij;
            else
                secondSum = productOfDeltaD2Epsilon2ij + secondSum;
            end
            
       end
   end
   
   T = firstSum + secondSum;
   
%   Calculates First U
   for i=1:nw
       for j=1:nx
            hj = generatesEi(i,nw); 
            productOfDeltaB1Beta1ij = (abs(deltaB1(j,i))^2)*betas1(j,i)*hj*hj';
            
%             First time i = j = 1
            if (i==1 && j == 1) 
                firstSum = productOfDeltaB1Beta1ij;
            else
                firstSum = productOfDeltaB1Beta1ij + firstSum;
            end
            
       end
   end

 %   Calculates Second U
   for i=1:nw
       for j=1:ny
            hj = generatesEi(i,nw); 
            productOfD1Beta1ij = (abs(deltaD1(j,i))^2)*epsilons1(j,i)*hj*hj';
            
%             First time i = j = 1
            if (i==1 && j == 1) 
                secondSum = productOfD1Beta1ij;
            else
                secondSum = productOfD1Beta1ij + secondSum;
            end
            
       end
   end
  
   
   U = firstSum +secondSum;

   
   
%   Xa and Lb2
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
    
    %   Xc, Ld and Id1
    for i=1:ny
        if(i==1)
            Xc = X;
            Ld2 = L'; 
           
        else
            Xc = [Xc X];
            Ld2 = [Ld2 L'];
            
        end
    end
    Xc = Xc';
    Ld2 = Ld2';
    
    
    %  Ib1
    for i=1:nw
        if(i==1)
            Ib1 = eye(nx);
            Id1 = eye(ny);
        else
            Ib1 = [Ib1 eye(nx)];
            Id1 = [Id1 eye(ny)];
        end
    end
    Ib1 = Ib1';
    Id1 = Id1';
    
    Alphas = diag(reshape(alphas',[nx*nx, 1]));
    Gammas = diag(reshape(gammas',[nx*ny, 1]));
    Betas2 = diag(reshape(betas2',[nu*nx, 1]));
    Betas1 = diag(reshape(betas1,[nw*nx, 1]));
    Epsilons2 = diag(reshape(epsilons2',[ny*nu, 1]));
    Epsilons1 = diag(reshape(epsilons1,[ny*nw, 1]));
    
    V = A0*X+B20*L;
    U11 = V+V'+S;
    U21 = C0*X+D20*L;
    U31 = B10';
    U41 = Xa;
    U51 = Xc;
    U61 = Lb2;
    U71 = Ld2;
    U81 = Ib1;
    U91 = zeros(size(Id1,1),size(U81,2));
    
    U22 = -eye(size(T))+T;
    U32 = D10';
    U32 = zeros(size(U31,1),size(U22,2));
    U42 = zeros(size(U41,1),size(U22,2));
    U52 = zeros(size(U51,1),size(U22,2));
    U62 = zeros(size(U61,1),size(U22,2));
    U72 = zeros(size(U71,1),size(U22,2));
    U82 = zeros(size(U81,1),size(U22,2));
    U92 = Id1;
    
    U33 = -mu*eye(size(U))+U;
    U43 = zeros(size(U41,1),size(U33,2));
    U53 = zeros(size(U51,1),size(U33,2));
    U63 = zeros(size(U61,1),size(U33,2));
    U73 = zeros(size(U71,1),size(U33,2));
    U83 = zeros(size(U81,1),size(U33,2));
    U93 = zeros(size(U91,1),size(U33,2));
    
    U44 = -Alphas;
    U54 = zeros(size(U51,1),size(U44,2));
    U64 = zeros(size(U61,1),size(U44,2));
    U74 = zeros(size(U71,1),size(U44,2));
    U84 = zeros(size(U81,1),size(U44,2));
    U94 = zeros(size(U91,1),size(U44,2));
    
    U55 = -Gammas;
    U65 = zeros(size(U61,1),size(U55,2));
    U75 = zeros(size(U71,1),size(U55,2));
    U85 = zeros(size(U81,1),size(U55,2));
    U95 = zeros(size(U91,1),size(U55,2));
    
    U66 = -Betas2;
    U66 = zeros(size(U61,1),size(U66,2));
    U76 = zeros(size(U71,1),size(U66,2));
    U86 = zeros(size(U81,1),size(U66,2));
    U96 = zeros(size(U91,1),size(U66,2));
    
    U77 = -Epsilons2;
    U87 = zeros(size(U81,1),size(U77,2));
    U97 = zeros(size(U91,1),size(U77,2));
    
    U88 = -Betas1;
    U98 = zeros(size(U91,1),size(U88,2));
    
    U99 = -Epsilons1;
    
    T = [U11 U21' U31' U41' U51' U61' U71' U81' U91';
         U21 U22  U32' U42' U52' U62' U72' U82' U92';
         U31 U32  U33  U43' U53' U63' U73' U83' U93';
         U41 U42  U43  U44  U54' U64' U74' U84' U94';
         U51 U52  U53  U54  U55  U65' U75' U85' U95';
         U61 U62  U63  U64  U65  U66  U76' U86' U96';
         U71 U72  U73  U74  U75  U76  U77  U87' U97';
         U81 U82  U83  U84  U85  U86  U87  U88  U98';
         U91 U92  U93  U94  U95  U96  U97  U98  U99];
     

    LMIs = LMIs + (T <= 0);
    outPut.line = outPut.line + size(T,1);
    
    outPut.cpusec_m = etime(clock,outPut.cpusec_m);
    outPut.var = size(getvariables(LMIs),2);
    
    sol = optimize(LMIs,obj,sdpsettings('verbose',0,'solver','sedumi'));
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