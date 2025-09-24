function outPut = synHInfKContIntLMILab(A,B2,B1,C,D2,D1,param)
    % Calculates K minimizing HInf norm, feed back control, from interval
    % matrices: A,B2,B1,C, D2 and D1
    % (continuos)

    if nargin == 7
        if isfield(param,'tol')
            tol = param.tol;
            episilonHinf =0.01;
        else
            tol = 1e-7;
            episilonHinf = 0.01;
        end  
    else
        tol = 1e-7;
        episilonHinf = 0.01;
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
      
    % Set the LMIs
    setlmis([]);
	% Set LMI count 
    c = 1;
    
    [X,nVar,sX] = lmivar(1,[nx 1]);
    lmiterm([-c 1 1 X],1,1); % X > 0
    c = c + 1;
    
    
     [L,nVar,sL] = lmivar(2,[nu nx]);
    
    %   Set obj
    if(nargin == 7 && isfield(param,'mu'))
        mu = param.mu;      
    else
        mu = lmivar(1,[1 0]);
    end
    
       
  %   Alphas 
    for i=1:nx
        for j=1:nx
            alphas(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end
    
 %    Betas 2
    for i=1:nx
        for j=1:nu
            betas2(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end
    
 %    Betas 1
    for i=1:nx
        for j=1:nw
            betas1(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end
    
  %    Gammas
    for i=1:ny
        for j=1:nx
            gammas(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end
    
%    Epsilons 2
    for i=1:ny
        for j=1:nu
            epsilons2(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end
    
%    Epsilons 1
    for i=1:ny
        for j=1:nw
            epsilons1(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end

   
%    ---------------------Specific Variable Declaration-------------------

    
%   Xa, Lb and Ib1
    for i=1:nx
        if(i==1)
            Xa = sX;
			Lb = sL'; 
			Ib1 = eye(nw);
        else
            Xa = [Xa sX];
            Lb = [Lb sL'];
            Ib1 = [Ib1 eye(nw)];
        end
    end
    [Xa,nVar,sXa] = lmivar(3,Xa');
    [Lb2,nVar,sLb] = lmivar(3,Lb');
    Ib1 = Ib1';
   
    %   Xc, Ld and Id1
    for i=1:ny
        if(i==1)
            Xc = sX;
			Ld = sL';
            Id1 = eye(nw);
        else
            Xc = [Xc sX];
			Ld = [Ld sL'];
            Id1 = [Id1 eye(nw)];
        end
    end
    [Xc,nVar,sXc] = lmivar(3,Xc');
    [Ld2,nVar,sLd] = lmivar(3,Ld');
    Id1 = Id1';
    
      
    
    %     First column
    lmiterm([c 1 1 X],A0,1,'s'); % A0*X+sym
	lmiterm([c 1 1 L],B20,1,'s'); % B0*L
	
    for i=1:nx
       for j=1:nx
            ei = generatesEi(i,nx); 
            lmiterm([c,1,1,alphas(i,j)],abs(deltaA(i,j))^2*ei*ei',1);
       end
   end
    for i=1:nx
        for j=1:nu
            ei = generatesEi(i,nx);
            lmiterm([c,1,1,betas2(i,j)],abs(deltaB2(i,j))^2*ei*ei',1);
        end
    end
    for i=1:nx
        for j=1:nw
            ei = generatesEi(i,nx);
            lmiterm([c,1,1,betas1(i,j)],abs(deltaB1(i,j))^2*ei*ei',1);
        end
    end
    lmiterm([c 2 1 X],C0,1); % C0*X
	lmiterm([c 2 1 L],D20,1); % D20*L
    lmiterm([c 3 1 0],B10'); % B10
	lmiterm([c 4 1 Xa],1,1); % Xa
    lmiterm([c 5 1 Xc],1,1); % Xc
    lmiterm([c 6 1 Lb2],1,1); % Lb
    lmiterm([c 7 1 Ld2],1,1); % Ld
    
    
    %     Second Column
   lmiterm([c 2 2 0],-1); % - I
   for i=1:ny
       for j=1:nx
            gi = generatesEi(i,ny); 
            lmiterm([c,2,2,gammas(i,j)],abs(deltaC(i,j))^2*gi*gi',1);
       end
   end
    for i=1:ny
        for j=1:nu
            gi = generatesEi(i,ny);
            lmiterm([c,2,2,epsilons2(i,j)],abs(deltaD2(i,j))^2*gi*gi',1);
        end
    end
     for i=1:ny
        for j=1:nw
            gi = generatesEi(i,ny);
            lmiterm([c,2,2,epsilons1(i,j)],abs(deltaD1(i,j))^2*gi*gi',1);
        end
     end
     lmiterm([c 3 2 0],D10'); % B10
    
  % Third Column 
    if(~(nargin == 7)||(nargin == 7 && ~isfield(param,'mu')))
        lmiterm([c 3 3 mu],-1,1);
    else
        lmiterm([c 3 3 0],-mu);
    end
    lmiterm([c 8 3 0],Ib1); % Ib1
    lmiterm([c 9 3 0],Id1); % Id1
    
    
    %     Fourth Column
   for i=1:nx
       for j=1:nx
           ei = generatesEi(((i-1)*nx)+j,(nx*nx)); 
           lmiterm([c,4,4,alphas(i,j)],-1*ei*ei',1);
       end
   end
   
   %    Fifth line
    for i=1:ny
       for j=1:nx
           ei = generatesEi(((i-1)*nx)+j,(ny*nx)); 
           lmiterm([c,5,5,gammas(i,j)],-1*ei*ei',1);
       end
    end
    
    %    Sixth line
    for i=1:nx
       for j=1:nu
           ei = generatesEi(((i-1)*nu)+j,(nx*nu)); 
           lmiterm([c,6,6,betas2(i,j)],-1*ei*ei',1);
       end
    end
    
      %    Seventh line
    for i=1:ny
       for j=1:nu
           ei = generatesEi(((i-1)*nu)+j,(ny*nu)); 
           lmiterm([c,7,7,epsilons2(i,j)],-1*ei*ei',1);
       end
    end
    
     %    Eigth line
    for i=1:nx
       for j=1:nw
           ei = generatesEi(((i-1)*nw)+j,(nx*nw)); 
           lmiterm([c,8,8,betas1(i,j)],-1*ei*ei',1);
       end
    end
    
       %    Nineth line
    for i=1:ny
       for j=1:nw
           ei = generatesEi(((i-1)*nw)+j,(ny*nw)); 
           lmiterm([c,9,9,epsilons1(i,j)],-1*ei*ei',1);
       end
    end
    
     c = c + 1;
     
    if((nargin == 7 && isfield(param,'mu')))
        [Y,nVar,sY] = lmivar(1,[nx 1]);
        lmiterm([-c 1 1 Y],1,1); % Y > 0
        lmiterm([-c 2 1 0],1); % I > 0
        lmiterm([-c 2 2 X],1,1); % X > 0
    end

    
    
    
    outPut.cpusec_m = etime(clock,outPut.cpusec_m);
    LMIs = getlmis;

    n = decnbr(LMIs);
    c = zeros(n,1);
    
    if(~(nargin == 7)||(nargin == 7 && ~isfield(param,'mu')))
         obj = mu;
        c(diag(decinfo(LMIs,obj))) = 1;
    else
        for j=1:n  
            [Yj] = defcx(LMIs,j,Y);
            c(j) = trace(Yj);
        end
    end
    

    
    options = [tol,2000,0,200,1];
    [copt,xopt] = mincx(LMIs,c,options);
    
    if isempty(xopt)
       outPut.feas = 0;
    else
       outPut.feas = 1;
        
       if(~(nargin == 7)||(nargin == 7 && ~isfield(param,'mu')))
           param.mu = dec2mat(LMIs,xopt,mu)+episilonHinf;
           synInt = synHInfKContIntLMILab(A,B2,B1,C,D2,D1,param);
           outPut.mu = sqrt(synInt.mu);
           outPut.X =  synInt.X;
           outPut.L = synInt.L;
           outPut.K = synInt.K;
           
           
           outPut.withoutExtLMI.mu = sqrt(dec2mat(LMIs,xopt,mu));
           outPut.withoutExtLMI.X =  dec2mat(LMIs,xopt,X);
           outPut.withoutExtLMI.L = dec2mat(LMIs,xopt,L);
           outPut.withoutExtLMI.K =  outPut.withoutExtLMI.L*outPut.withoutExtLMI.X^-1;
           
           
       else
           outPut.mu = mu;
           outPut.X =  dec2mat(LMIs,xopt,X);
           outPut.L = dec2mat(LMIs,xopt,L);
           outPut.K = outPut.L*outPut.X^-1;
       end
       
      
    end
end