function outPut = synHInfKDiscPoly(APoly,B2Poly,B1Poly,CPoly,D2Poly,D1Poly,param)
    % Calculates hinf norm based on lema 32 from Ricardo's slides 

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
    %LMI rows counter
    outPut.line = 0;
    
    %   New LMI system
    LMIs = set([]);
    G = sdpvar(nx,nx,'full');
    Z = sdpvar(nu,nx,'full');
   
    
    outPut.line = outPut.line + nx;
  
 
      %   Set obj
%     if(nargin == 7 && isfield(param,'mu'))
%          Y = sdpvar(nx,nx,'symmetric');
%          U11 = Y;
%          U21 = eye(nx);
%          U22 = W;
%          T = [U11    U21';
%               U21    U22];
%          LMIs = LMIs + (T >= 0);
%          
%         mu = param.mu;   
%         obj = trace(Y);
%     else
        mu = sdpvar(1,1);
        obj = mu;
%     end
%     

    for i=1:N
        W{i} = sdpvar(nx,nx,'symmetric');
         
        U32 = CPoly{i}*G+D2Poly{i}*Z;
        
        U11 = W{i};
        U21 = G'*(APoly{i})'+Z'*(B2Poly{i})';
        U31 = zeros(size(U32,1),size(U21,2));
        U41 = (B1Poly{i})';
        
        U22 = G+G'-W{i};
%       U32  
        U42 = zeros(size(U41,1),size(U32,2));
        
        U33 = eye(size(U32,1));
        U43 = (D1Poly{i})';
        
        U44 = mu*eye(size(U43,1));

        T = [U11    U21'     U31'   U41';
             U21    U22     U32'    U42'; 
             U31    U32     U33     U43';
             U41    U42     U43     U44];
        LMIs = LMIs + (T >= 0);
        
        outPut.line = outPut.line + size(T,1);
       
    end
    
   
    outPut.cpusec_m = etime(clock,outPut.cpusec_m);
    outPut.var = size(getvariables(LMIs),2);
    
    sol = optimize(LMIs,obj,sdpsettings('verbose',0,'solver','sedumi'));
    
    
    outPut.cpusec = sol.solvertime;
    p = min(checkset(LMIs));
    outPut.delta = p;

    if p > -tol
        outPut.feas = 1;
        
%         if(~(nargin == 7)||(nargin == 7 && ~isfield(param,'mu')))
%            param.mu = double(mu)+episilonHinf;
%            synPoly = synHInfKDiscPoly(APoly,B2Poly,B1Poly,CPoly,D2Poly,D1Poly,param);
%            outPut.mu = sqrt(synPoly.mu);
%            outPut.Z =  synPoly.Z;
%            outPut.G = synPoly.G;
%            outPut.K = synPoly.K;
%            
%            outPut.withoutExtLMI.Z = double(Z);
%            outPut.withoutExtLMI.G =  double(G);
%            outPut.withoutExtLMI.K = outPut.withoutExtLMI.Z*outPut.withoutExtLMI.G^-1;
%            outPut.withoutExtLMI.mu = sqrt(double(mu));
%           
%            
%            outPut.withoutExtLMI.K =  outPut.withoutExtLMI.Z*outPut.withoutExtLMI.G^-1;
%            
%            
%         else
           outPut.Z = double(Z);
           outPut.G =  double(G);
           outPut.K = outPut.Z*outPut.G^-1;
           outPut.norm = double(sqrt(mu)); 
%         end
%        
        
        
    else
        outPut.feas = 0;        
    end
end