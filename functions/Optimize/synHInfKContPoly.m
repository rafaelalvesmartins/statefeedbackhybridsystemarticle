function outPut = synHInfKContPoly(APoly,B2Poly,B1Poly,CPoly,D2Poly,D1Poly,param)
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
    %LMI rows counter
    outPut.line = 0;
    
    %   New LMI system
    LMIs = set([]);

    Z = sdpvar(nu,nx,'full');
    W = sdpvar(nx,nx,'symmetric');
    LMIs = LMIs + (W >= 0);
    outPut.line = outPut.line + nx;
  
 
      %   Set obj
    if(nargin == 7 && isfield(param,'mu'))
         Y = sdpvar(nx,nx,'symmetric');
         U11 = Y;
         U21 = eye(nx);
         U22 = W;
         T = [U11    U21';
              U21    U22];
         LMIs = LMIs + (T >= 0);
         
        mu = param.mu;   
        obj = trace(Y);
    else
        mu = sdpvar(1,1);
        obj = mu;
    end
    
   
    

    for i=1:N
        V = APoly{i}*W+B2Poly{i}*Z;
        U11 = V'+V;
        U21 = CPoly{i}*W+D2Poly{i}*Z;
        U31 = (B1Poly{i})';

        U22 = -eye(size(U21,1));
        U32 = (D1Poly{i})';
        U33 = -mu*eye(size(U32,1));

        T = [U11    U21'     U31';
             U21    U22     U32'; 
             U31    U32     U33];
        LMIs = LMIs + (T <= 0);
        
        outPut.line = outPut.line + size(T,1);
       
    end
    
  
    
   
    outPut.cpusec_m = etime(clock,outPut.cpusec_m);
    outPut.var = size(getvariables(LMIs),2);
    
    % sol = optimize(LMIs,obj,sdpsettings('shift','1','verbose',0,'solver','sedumi','sedumi.eps',1e-20,'sedumi.maxiter',2000,'sedumi.numtol',1e-20));
     sol = optimize(LMIs,obj,sdpsettings('verbose',0,'solver','mosek'));
    
    
    outPut.cpusec = sol.solvertime;
    p = min(checkset(LMIs));
    outPut.delta = p;

    if p > -tol
        outPut.feas = 1;
        outPut.Z = double(Z);
        outPut.W =  double(W);
        outPut.K = outPut.Z*outPut.W^-1;
        outPut.norm = sqrt(double(mu)); 
           
        if(~(nargin == 7)||(nargin == 7 && ~isfield(param,'mu')))
           if(~isfield(param,'mu4SmallK')) 
               episilonHinf = 0.001;  
               K = outPut.K;
               maxNormOfK = param.maxNormOfK;
               while(norm(K)>maxNormOfK)
                   param.mu = double(mu)+episilonHinf;
                   synPoly = synHInfKContPoly(APoly,B2Poly,B1Poly,CPoly,D2Poly,D1Poly,param);
                   if(~synPoly.feas)
                       break;
                   end
                   K = synPoly.K;
                   episilonHinf = episilonHinf*1.5;
               end
               
           else
               param.mu = param.mu4SmallK;
               synPoly = synHInfKContPoly(APoly,B2Poly,B1Poly,CPoly,D2Poly,D1Poly,param);
           end
           
         
           
           if(exist('synPoly') && synPoly.feas)
               smallK.feas = 1;
               smallK.norm = synPoly.norm;
               smallK.W =  synPoly.W;
               smallK.Z = synPoly.Z;
               smallK.K = synPoly.K;
           else
               smallK.feas = 0;
           end
           outPut.smallK = smallK;             
        end
    else
        outPut.feas = 0;        
    end
end