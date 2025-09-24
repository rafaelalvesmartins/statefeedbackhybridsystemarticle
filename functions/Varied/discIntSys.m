function [Ad,Bd,Cd,Dd] = discIntSys(A,B,C,D,h)
    nx = length(A);
    nu = size(B,2);
    
    Acal = [A B; zeros(nu,nx+nu)];
    Qcal = [C'*C C'*D;D'*C D'*D];
%     Rcal = chol(Qcal);
    
    % Scaling and squaring Exponential 
    K = 10;
    L = 10;
    
    % Discretize A e B first
    expDisc =  expIntMatScalingSquaring(Acal*h,K,L);
    Ad = expDisc(1:nx,1:nx);
    Bd = expDisc(1:nx,nx+1:end);
    
    f1 = @(t) Acal*t; 
    f2 = @(t) expIntMatScalingSquaring(f1(t),K,L);
%     f3 = @(t) Rcal*f2(t);
%     f4 = @(t) f3(t)'*f3(t);
  f4 = @(t) f2(t)'*Qcal*f2(t);

    a = 0;
    b = h;
    n = 500;
    
    Qbar = integral_sum(infsup(a,b),f4,n);
    %     Qbar = integralSimpsons(f3,a,b,n);
    L=decCholesky(Qbar);

    Cd = (L(1:nx,1:nx+nu))';
    Dd = (L(nx+1:nx+nu,1:nx+nu))';
end