function outPut = synHInfKContPolyLMILab(APoly,B2Poly,B1Poly,CPoly,D2Poly,D1Poly,param)
    % Calculates hinf norm based on lema 31 from Ricardo's slides 

   if nargin == 7
        if isfield(param,'tol')
            tol = param.tol;
        else
            tol = 1e-7;
        end    
    else
        tol = 1e-7;
    end
  
    
    %    ---------------------Dimension of input----------------------------
    nx = length(APoly{1});
    nu = size(B2Poly{1},2);

    N = length(APoly);
    
    outPut.cpusec_m = clock;
    % Set the LMIs
    setlmis([]);
    %     Set LMI count 
    c = 1;
       
    %LMI rows counter
    outPut.line = 0;
    
    [W,nVar,sW] = lmivar(1,[nx 1]);
    lmiterm([-c 1 1 W],1,1); % W > 0
    c = c + 1;
    outPut.line = outPut.line + nx;
    
    [Z,nVar,sZ] = lmivar(2,[nu nx]);
    
     %   Set obj
    mu = lmivar(1,[1 0]);
    obj = mu;

    for i=1:N
        lmiterm([c 1 1 W],APoly{i},1,'s');
        lmiterm([c 1 1 Z],B2Poly{i},1,'s');
        outPut.line = outPut.line + nx;
        
        lmiterm([c 2 1 W],CPoly{i},1);
        outPut.line = outPut.line + nx;
        
        lmiterm([c 2 1 Z],D2Poly{i},1);
        outPut.line = outPut.line + nx;
        
        lmiterm([c 3 1 1],B1Poly{i}',1);
        outPut.line = outPut.line + nx;
        
        lmiterm([c 2 2 0],-1);
        lmiterm([c 3 2 1],D1Poly{i}',1);
        lmiterm([c 3 3 mu],-1,1);
        
        c = c + 1;
    end
    LMIs = getlmis;
    c = zeros(decnbr(LMIs),1);
    c(diag(decinfo(LMIs,obj))) = 1;
    options = [tol,10000,0,200,1];
    [copt,xopt] = mincx(LMIs,c,options);
    
    
    if isempty(xopt)
       outPut.feas = 0;
    else
        outPut.feas = 1;
        outPut.Z = dec2mat(LMIs,xopt,Z);
        outPut.W =  dec2mat(LMIs,xopt,W);
        outPut.K = outPut.Z*outPut.W^-1;
        outPut.mu = dec2mat(LMIs,xopt,mu);
    end
end