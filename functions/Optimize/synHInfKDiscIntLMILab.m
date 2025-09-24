function outPut = synHInfKDiscIntLMILab(A,B2,B1,C,D2,D1,param)
    % Calculates K minimizing HInf norm, feed back control, from interval
    % matrices: A,B2,B1,C, D2 and D1
    % (discrete)

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
      
    % Set the LMIs
    setlmis([]);
	% Set LMI count 
    c = 1;
    
    [X,nVar,sX] = lmivar(1,[nx 1]);
    lmiterm([-c 1 1 X],1,1); % X > 0
    c = c + 1;
    
    
     [L,nVar,sL] = lmivar(2,[nu nx]);
    
    %   Set obj
        mu = lmivar(1,[1 0]);
  
       
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
    lmiterm([-c 1 1 X],1,1); 
   for i=1:nx
       for j=1:nx
            ei = generatesEi(i,nx); 
            lmiterm([-c,1,1,alphas(i,j)],-abs(deltaA(i,j))^(2)*(ei*ei'),1);
       end
   end
    for i=1:nx
        for j=1:nu
            ei = generatesEi(i,nx);
            lmiterm([-c,1,1,betas2(i,j)],-abs(deltaB2(i,j))^(2)*(ei*ei'),1);
        end
    end
    for i=1:nx
        for j=1:nw
            ei = generatesEi(i,nx);
            lmiterm([-c,1,1,betas1(i,j)],-abs(deltaB1(i,j))^(2)*(ei*ei'),1);
        end
    end
    lmiterm([-c 2 1 X],1,A0'); 
	lmiterm([-c 2 1 -L],1,B20');
	lmiterm([-c 4 1 0],B10'); % B10'
   
    
    %     Second Column
   lmiterm([-c 2 2 X],1,1); 
   lmiterm([-c 3 2 X],C0,1); % C0*X
   lmiterm([-c 3 2 L],D20,1); % D20*L
   lmiterm([-c 5 2 Xa],1,1); % Xa
   lmiterm([-c 6 2 Xc],1,1); % Xc
   lmiterm([-c 7 2 Lb2],1,1); % Lb
   lmiterm([-c 8 2 Ld2],1,1); % Ld
    
  % Third Column 
   lmiterm([-c 3 3 0],1); % I
   for i=1:ny
       for j=1:nx
            gi = generatesEi(i,ny); 
            lmiterm([-c,3,3,gammas(i,j)],-abs(deltaC(i,j))^(2)*(gi*gi'),1);
       end
   end
    for i=1:ny
        for j=1:nu
            gi = generatesEi(i,ny);
            lmiterm([-c,3,3,epsilons2(i,j)],-abs(deltaD2(i,j))^(2)*(gi*gi'),1);
        end
    end
     for i=1:ny
        for j=1:nw
            gi = generatesEi(i,ny);
            lmiterm([-c,3,3,epsilons1(i,j)],-abs(deltaD1(i,j))^(2)*(gi*gi'),1);
        end
     end
    lmiterm([-c 4 3 0],D10'); 
        
    %     Fourth Column
   lmiterm([-c 4 4 mu],1,1);
   lmiterm([-c 9 4 0],Ib1); % Ib1
   lmiterm([-c 10 4 0],Id1); % Id1
   
    
   %    Fifth line
   for i=1:nx
       for j=1:nx
           ei = generatesEi(((i-1)*nx)+j,(nx*nx)); 
           lmiterm([-c,5,5,alphas(i,j)],1*ei*ei',1);
       end
   end
   
    
    %    Sixth line
      for i=1:ny
       for j=1:nx
           ei = generatesEi(((i-1)*nx)+j,(ny*nx)); 
           lmiterm([-c,6,6,gammas(i,j)],1*ei*ei',1);
       end
      end
    
    
      %    Seventh line
    for i=1:nx
       for j=1:nu
           ei = generatesEi(((i-1)*nu)+j,(nx*nu)); 
           lmiterm([-c,7,7,betas2(i,j)],1*ei*ei',1);
       end
    end
    
    
     %    Eigth line
     for i=1:ny
       for j=1:nu
           ei = generatesEi(((i-1)*nu)+j,(ny*nu)); 
           lmiterm([-c,8,8,epsilons2(i,j)],1*ei*ei',1);
       end
    end
   
    
       %    Nineth line
     for i=1:nx
       for j=1:nw
           ei = generatesEi(((i-1)*nw)+j,(nx*nw)); 
           lmiterm([-c,9,9,betas1(i,j)],1*ei*ei',1);
       end
     end
    
%    Tenth line
    for i=1:ny
       for j=1:nw
           ei = generatesEi(((i-1)*nw)+j,(ny*nw)); 
           lmiterm([-c,10,10,epsilons1(i,j)],1*ei*ei',1);
       end
    end
    
    
        
    
    
    outPut.cpusec_m = etime(clock,outPut.cpusec_m);
    LMIs = getlmis;

    n = decnbr(LMIs);
    c = zeros(n,1);
    
    obj = mu;
    c(diag(decinfo(LMIs,obj))) = 1;

    
    options = [tol,2000,0,200,1];
    [copt,xopt] = mincx(LMIs,c,options);
    
    if isempty(xopt)
       outPut.feas = 0;
    else
       outPut.feas = 1;
        
       
   outPut.mu = sqrt(dec2mat(LMIs,xopt,mu));
   outPut.X =  dec2mat(LMIs,xopt,X);
   outPut.L = dec2mat(LMIs,xopt,L);
   outPut.K = outPut.L*outPut.X^-1;
           
       
       
      
    end
end