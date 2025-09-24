function output = calcKHInfIntSeDuMi(A,B2,B1,C,D2,D1,param)
% Calculates K, feed back control, from interval matrices: A,B1,B2,C and D
% (continuos)

    if nargin == 6
        if isfield(param,'mu')
            mu = param.mu;
        else
            mu = 0;
        end
        if isfield(param,'tol')
            tol = param.tol;
        else
            tol = 1e-7;
        end    
    else
        mu = 0;
        tol = 1e-7;
    end
    
    

%    ---------------------Usefull matrices----------------------------
    A0 = (A.sup+A.inf)/2;
    deltaA = (A.sup-A.inf)/2;
    
    B0_2 = (B2.sup+B2.inf)/2;
    deltaB_2 = (B2.sup-B2.inf)/2;
    
    B0_1 = (B1.sup+B1.inf)/2;
    deltaB_1 = (B1.sup-B1.inf)/2;
    
    C0 = (C.sup+C.inf)/2;
    deltaC = (C.sup-C.inf)/2;

    D0_2 = (D2.sup+D2.inf)/2;
    deltaD_2 = (D2.sup-D2.inf)/2;
    
    D0_1 = (D1.sup+D1.inf)/2;
    deltaD_1 = (D1.sup-D1.inf)/2;

%    ---------------------Dimension of input----------------------------
    nx = length(A.inf);
    nu = size(B2.inf,2);
    ny = size(C.inf,1);
    
%    ---------------------Variable Declariation----------------------------

   %New LMI system
    LMIs = set([]);
    
    X = sdpvar(nx,nx,'symmetric');
    LMIs = LMIs + (X >= 0);
    L = sdpvar(nu,nx,'full');
    
%     Set obj
    if mu == 0
        mu = sdpvar(1,1);
        obj = mu;
    else
%         mu = mu;
        obj = trace(X);
    end
    
    
%   Alphas 
    for i=1:nx
        for j=1:nx
            alphas(i,j) = sdpvar(1,1); % Scalar
            LMIs = LMIs + (alphas(i,j) >= 0);
        end
    end
    
%    Betas 1
    for i=1:nx
        for j=1:nu
            betas_1(i,j) = sdpvar(1,1); % Scalar
            LMIs = LMIs + (betas_1(i,j) >= 0);
        end
    end
    
    %    Betas 2
    for i=1:nx
        for j=1:nu
            betas_2(i,j) = sdpvar(1,1); % Scalar
            LMIs = LMIs + (betas_2(i,j) >= 0);
        end
    end
    
    %    Gammas
    for i=1:ny
        for j=1:nx
            gammas(i,j) = sdpvar(1,1); % Scalar
            LMIs = LMIs + (gammas(i,j) >= 0);
        end
    end
    
    %    Epsilons 2
    for i=1:ny
        for j=1:nu
            epsilons_2(i,j) = sdpvar(1,1); % Scalar
            LMIs = LMIs + (epsilons_2(i,j) >= 0);
        end
    end
    
     %    Epsilons 1
    for i=1:ny
        for j=1:nu
            epsilons_1(i,j) = sdpvar(1,1); % Scalar
            LMIs = LMIs + (epsilons_1(i,j) >= 0);
        end
    end

%    ---------------------Specific Variable Declaration-------------------

% Calculates Q
%     First Q' sum
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
    
%     Second Q' sum
   for i=1:nx
       for j=1:nu
            
            ei = generatesEi(i,nx); 
            productOfDeltaBBeta_2ij = (abs(deltaB_2(i,j))^2)*betas_2(i,j)*ei*ei';
            
%             First time i = j = 1
            if (i==1 && j == 1) 
                secondSum = productOfDeltaBBeta_2ij;
            else
                secondSum = productOfDeltaBBeta_2ij + secondSum;
            end
            
       end
   end
   
   %     Third Q' sum
   for i=1:nx
       for j=1:nu
            
            ei = generatesEi(i,nx); 
            productOfDeltaBBeta_1ij = (abs(deltaB_1(i,j))^2)*betas_1(i,j)*ei*ei';
            
%             First time i = j = 1
            if (i==1 && j == 1) 
                thirdSum = productOfDeltaBBeta_1ij;
            else
                thirdSum = productOfDeltaBBeta_1ij + thirdSum;
            end
            
       end
   end
   
   Q = firstSum + secondSum + thirdSum;
    
   
   
   % Calculates R
%     First R' sum
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
    
%     Second R' sum
   for i=1:ny
       for j=1:nu
            
            gi = generatesEi(i,ny); 
            productOfDeltaD2Epsilonij = (abs(deltaD_2(i,j))^2)*epsilons_2(i,j)*gi*gi';
            
%             First time i = j = 1
            if (i==1 && j == 1) 
                secondSum = productOfDeltaD2Epsilonij;
            else
                secondSum = productOfDeltaD2Epsilonij + secondSum;
            end
            
       end
   end
   
   %     Third R' sum
   for i=1:ny
       for j=1:nu
            
            gi = generatesEi(i,ny); 
            productOfDeltaD1Epsilonij = (abs(deltaD_1(i,j))^2)*epsilons_1(i,j)*gi*gi';
            
%             First time i = j = 1
            if (i==1 && j == 1) 
                thirdSum = productOfDeltaD1Epsilonij;
            else
                thirdSum = productOfDeltaD1Epsilonij + thirdSum;
            end
            
       end
   end
   
   R = firstSum + secondSum + thirdSum;
   
   
%   Xa and Lb
    for i=1:nx
        if(i==1)
            Xa = X;
            Lb = L'; % CONFIRMAR SE NÃO TEM TRANSPOSTO
        else
            Xa = [Xa X];
            Lb = [Lb L'];
        end
    end
    Xa = Xa';
    Lb = Lb';
    
    
    %   Xc and Ld
    for i=1:ny
        if(i==1)
            Xc = X;
            Ld = L'; % CONFIRMAR SE NÃO TEM TRANSPOSTO
            Iq = eye(nu);
        else
            Xc = [Xc X];
            Ld = [Ld L'];
            Iq = [Iq eye(nu)];
        end
    end
    Xc = Xc';
    Ld = Ld';
    Iq = Iq';
    
    Alphas = diag(reshape(alphas',[nx*nx, 1]));
    Gammas = diag(reshape(gammas',[nx*ny, 1]));
    Betas_2 = diag(reshape(betas_2',[nu*nx, 1]));
    Betas_1= diag(reshape(betas_1',[nu*nx, 1]));
    Epsilons_2 = diag(reshape(epsilons_1',[nu*ny, 1]));
    Epsilons_1 = diag(reshape(epsilons_2',[nu*ny, 1]));
    
    V = A0*X+B0_2*L;
    U11 = V+V'+Q;
    U21 = C0*X+D0_2*L;
    U31 = B0_1';
    U41 = Xa;
    U51 = Xc;
    U61 = Lb;
    U71 = Lb;
    U81 = Ld;
    U91 = zeros(ny*ny,size(U81,2));
    
    U22 = -eye(length(R))+R;
    U32 = D0_1';
    U42 = zeros(size(U41,1),size(U22,2));
    U52 = zeros(size(U51,1),size(U22,2));
    U62 = zeros(size(U61,1),size(U22,2));
    U72 = zeros(size(U71,1),size(U22,2));
    U82 = zeros(size(U81,1),size(U22,2));
    U92 = zeros(size(U91,1),size(U22,2));
    
    U33 = -mu*eye(nu);
    U43 = zeros(size(U41,1),size(U33,2));
    U53 = zeros(size(U51,1),size(U33,2));
    U63 = zeros(size(U61,1),size(U33,2));
    U73 = zeros(size(U71,1),size(U33,2));
    U83 = zeros(size(U81,1),size(U33,2));
    U93 = Iq;
    
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
    
    U66 = -Betas_2;
    U76 = zeros(size(U71,1),size(U66,2));
    U86 = zeros(size(U81,1),size(U66,2));
    U96 = zeros(size(U91,1),size(U66,2));
    
    U77 = -Betas_1;
    U87 = zeros(size(U81,1),size(U77,2));
    U97 = zeros(size(U91,1),size(U77,2));
    
    U88 = -Epsilons_2;
    U98 = zeros(size(U91,1),size(U88,2));
   
    U99 = -Epsilons_1;
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
    
    optimize(LMIs,obj,sdpsettings('verbose',0,'solver','sedumi'));

    p = min(checkset(LMIs));
    output.delta = p;
    
    output.feas = 0;
    if(p > -tol ) % Is feasible
        output.L = double(L);
        output.X = double(X);
        output.k = output.L*output.X^-1;
        output.mu = sqrt(double(mu));
        output.feas = 1;
    end
end