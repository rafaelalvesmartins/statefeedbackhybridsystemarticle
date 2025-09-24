% Unicamp - FEEC - 25/01/2018
close all;clear;clc;

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

% 
% A.sup = A.inf;
% B2.sup = B2.inf;
% B1.sup = B1.inf;
% C.sup = C.inf;
% D2.sup = D2.inf;
% D1.sup = D1.inf;

verticesAPoly = genVerticesPolyFromUncerMat(A,prec);
verticesB2Poly = genVerticesPolyFromUncerMat(B2,prec);
verticesB1Poly = genVerticesPolyFromUncerMat(B1,prec);
verticesCPoly = genVerticesPolyFromUncerMat(C,prec);
verticesD2Poly = genVerticesPolyFromUncerMat(D2,prec);
verticesD1Poly = genVerticesPolyFromUncerMat(D1,prec);

x0 = ones(2,1);

nx = length(A.inf);
nu = size(B2.inf,2);
ny = size(C.inf,1);

param.tol = 1e-7;
outputSeDuMi = calculatesKHInfRobust(verticesAPoly,verticesB2Poly,verticesB1Poly,verticesCPoly,verticesD2Poly,verticesD1Poly,param);
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



numPointsUniSpaced = 1000;
numPointsBy2Points = 100;
numbPointsUniSpacedSub = 150;

APolytopic = genPoly(verticesAPoly, numPointsUniSpaced, numPointsBy2Points, numbPointsUniSpacedSub);
B2Polytopic = genPoly(verticesB2Poly, numPointsUniSpaced, numPointsBy2Points, numbPointsUniSpacedSub);
B1Polytopic = genPoly(verticesB1Poly, numPointsUniSpaced, numPointsBy2Points, numbPointsUniSpacedSub);
CPolytopic = genPoly(verticesCPoly, numPointsUniSpaced, numPointsBy2Points, numbPointsUniSpacedSub);
D2Polytopic = genPoly(verticesD2Poly, numPointsUniSpaced, numPointsBy2Points, numbPointsUniSpacedSub);
D1Polytopic = genPoly(verticesD1Poly, numPointsUniSpaced, numPointsBy2Points, numbPointsUniSpacedSub);


% A close loop
ACl= cellfun(@(x,y) x+y*k, APolytopic.polytopicMatrices,B2Polytopic.polytopicMatrices, 'UniformOutput',false);
CCl= cellfun(@(x,y) x+y*k, CPolytopic.polytopicMatrices,D2Polytopic.polytopicMatrices, 'UniformOutput',false);


SYSPoly= cellfun(@ss,ACl,B1Polytopic.polytopicMatrices,CCl,D1Polytopic.polytopicMatrices, 'UniformOutput',false);

custs = cell2mat(cellfun(@(x) norm(x,inf), SYSPoly,'UniformOutput',false));

averageCust = mean(custs)
stdDeviatCust = std(custs)
maxCust = max(custs)



