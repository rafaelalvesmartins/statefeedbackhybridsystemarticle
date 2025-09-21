function gerarLMIAlternadaDerivadaSint(A,B,E,C,D,S,Z,L,mu,i,c)
    lmiterm([c 1 1 S{i}],1,A','s');
    lmiterm([c 1 1 L{i}],1,B','s');
    lmiterm([c 2 1 -L{i}],1,A');
    lmiterm([c 2 1 -Z{i}],1,B');
    lmiterm([c 3 1 S{i}],C,1);
    lmiterm([c 3 1 -L{i}],D,1);
    lmiterm([c 3 2 L{i}],C,1);
    lmiterm([c 3 2 Z{i}],D,1);
    lmiterm([c 3 3 0],-1);
    lmiterm([c 4 1 0],E');
    lmiterm([c 4 4 mu],-1,1);
end