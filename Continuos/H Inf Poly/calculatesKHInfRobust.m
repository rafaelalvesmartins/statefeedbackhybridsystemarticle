function output = calculatesKHInfRobust(A,B2,B1,C,D2,D1,param)
% Return de control based on Lema 31 - Robust HInf norm
 if nargin == 7 && isfield(param,'mu')
        obj = [];
        muSquared = param.mu^2;
    else
        muSquared=sdpvar();
        obj = muSquared;
 end
% -------------------Get the dimenstion--------------------------------- 
    
%   Matrices' Length
    lengthA = length(A{1});
    lengthB1 = size(B1{1},2);
    lengthB2 = size(B2{1},2);
    lengthC = size(C{1},1);
   
% -------------------Variable declaration---------------------------------       
    
%     Vars declaration
    Z = sdpvar(lengthB2,lengthA,'full');
    W = sdpvar(lengthA,lengthA,'symmetric');
   
%   Build the matrices and LMIs    
    
    %new LMI system
    LMIs = set([]);
    
%   1St LMI
    LMIs = LMIs + (W>=0);

%   2Nd LMI
    for i=1:length(A)
        for j = 1:length(B2)
            for k = 1:length(B1)
                for l = 1:length(C)
                    for m = 1:length(D2)
                        for n = 1:length(D1)
                            
                            U11 = A{i}*W + W*A{i}' + B2{j}*Z  + Z'*B2{j}';
                            U12 = W*C{l}' + Z'*D2{m}';
                            U13 = B1{k};
                            U22 = -eye(lengthC);
                            U23 = D1{n};
                            U33 = -muSquared*eye(lengthB1);

                            U = [U11    U12     U13;
                                 U12'   U22     U23
                                 U13'   U23'    U33];

                            LMIs = LMIs + (U<=0);
                        end
                    end
                end
            end
        end
    end
 

    % -------------------Solve the LMIs---------------------------------

    optimize(LMIs,obj,sdpsettings('verbose',0,'solver','sedumi'));
    
    p = min(checkset(LMIs));
    tol = 1e-7;
    output.feas = 0;
    if(tol > -p ) % Is feasible
        output.Z = double(Z);
        output.W = double(W);
        output.k = output.Z*output.W^-1;
        output.mu = sqrt(double(muSquared));
        output.feas = 1;
    end
end
