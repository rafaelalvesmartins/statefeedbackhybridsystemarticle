function riccatiOutput =  calcKRiccatiHinf(A,B2,B1,C,D2,D1,gamma,param)
  

	if(nargin == 8 && isfield(param,'disc'))
        
%         FURUTA, Katsuhisa; PHOOJARUENCHANACHAI, Suthee. An algebraic approach to discrete-time H? control problems. In: American Control Conference, 1990. IEEE, 1990. p. 3067-3072. 
        B = [B1 B2];
        S =[C'*D2 C'*D1];
        R = [D1'*D1-gamma^2*eye(size(D1,2)) D1'*D2;D2'*D1 D2'*D2];
        Q = C'*C;
		[X] = dare(A,B,Q,R,S);
        B1Ch = B1'*X*A+D1'*C;
        B2Ch = B2'*X*A+D2'*C;
        B3Ch = B1'*X*B2+D1'*D2;
        
        theta1 = eye(size(B1,2))*gamma^2-B1'*X*B1-D1'*D1;
        theta2 = D2'*D2+B2'*X*B2+B3Ch'*inv(theta1)*B3Ch;
        theta3 = B2Ch+B3Ch'*inv(theta1)*B1Ch;
        
        K = -inv(theta2)*theta3;
        
    else
         B = [B2 B1];
        Q = C'*C;
        R = [D2'*D2 D2'*D1;D1'*D2 D1'*D1-gamma^2*eye(size(D1,2))];
        S = [C'*D2 C'*D1];
	    [X] = care(A,B,Q,R,S);
% 		K = -inv(D2'*D2)*(D2'*C+B2'*X);
        mat1 = -inv(D2'*inv(eye(size(D1,1))-gamma^(-2)*D1*D1')*D2);
        mat2 = D2'*C+B2'*X+D2'*D1*inv(gamma^2*eye(size(D1',1))-D1'*D1)*(D1'*C+B1'*X);
        K = mat1*mat2;
    end
    
    riccatiOutput.X = X;
    riccatiOutput.K = K;
end