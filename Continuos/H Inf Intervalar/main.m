% Unicamp - FEEC - 25/01/2018
close all;clear;clc;
addpath('../Aux functions');

A.inf = [30 80; 90 200];
A.sup = A.inf + [0.1 0; 0 0.1];

prec = 1e-4;

B2.inf = [1 1]';
B2.sup = B2.inf + [1 1]';

B1.inf = [1 1]';
B1.sup = B1.inf + [1 1]';

C.inf = [1 1];
C.sup = C.inf + [2 1];

D2.inf = 1;
D2.sup = D2.inf + 0.1;

D1.inf = 0;
D1.sup = D1.inf;

% A.sup = A.inf;
% B1.sup = B1.inf;
% B2.sup = B2.inf;
% C.sup = C.inf;
% D.sup = D.inf;



param.tol = 1e-7;
outputSeDuMi = calcKHInfIntSeDuMi(A,B2,B1,C,D2,D1,param);
outputSeDuMi.mu
% 
% outputLMILAB = calcKStateFeedBackContIntABCDLMILAB(A,B,C,D,x0,param);
% outputLMILAB.mu

k = outputSeDuMi.k;

% outputLMILABRafael = calcKStateFeedBackContIntABCDSeDuMiRafael(A,B,C,D,x0,param);
% outputLMILABRafael.mu
% 
% sys2 = ss((A.inf+B.inf*k),x0,(C.inf+D.inf*k),zeros(ny,nu));
% norm(sys2)^2
% 
% Q = C.inf'*C.inf;
% R = D.inf'*D.inf;
% N = C.inf'*D.inf;
% 
% [K,S,E] = lqr(A.inf,B.inf,Q,R,N);
% x0'*S*x0
% 
% 
%     [eig(A.inf+B.inf*K)  eig(A.inf+B.sup*K) eig(A.sup+B.if*K) eig(A.sup+B.sup*K)] >= 0
%     [eig(A.inf+B.inf*K)  eig(A.inf+B.sup*K) eig(A.sup+B.inf*K) eig(A.sup+B.sup*K)] >= 0
%     

    
nSimulation = 1000;
custVector = zeros(nSimulation,1);
for i=1:nSimulation
    AIn = genRandMatInTheInterv(A);
    B2In = genRandMatInTheInterv(B2);
    B1In = genRandMatInTheInterv(B1);
    CIn = genRandMatInTheInterv(C);
    D2In = genRandMatInTheInterv(D2);
    D1In = genRandMatInTheInterv(D1);

    sys = ss((AIn+B2In*k),B1In,(CIn+D1In*k),D1In);
    custVector(i) = norm(sys,inf);
    
end

averageCust = mean(custVector)
stdDeviatCust = std(custVector)
maxCust = max(custVector)




