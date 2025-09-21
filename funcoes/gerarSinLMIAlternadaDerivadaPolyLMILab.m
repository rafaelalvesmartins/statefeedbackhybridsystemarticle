function gerarSinLMIAlternadaDerivadaPolyLMILab(A,B,E,C,D,S,Z,L,X,mu,c)
                                                     
     
    lmiterm([c 1 1 X.X1_11],1,A','s');
    lmiterm([c 1 1 X.X1_12],1,B','s');
    
    lmiterm([c 2 1 X.X1_21],1,A');
    lmiterm([c 2 1 X.X1_22],1,B');
    
    lmiterm([c 3 1 0],E');
    lmiterm([c 3 1 X.X2_11],1,A');
    lmiterm([c 3 1 X.X2_12],1,B');
    lmiterm([c 3 3 mu],-1,1);
    
    lmiterm([c 4 1 X.X3_11],1,A');
    lmiterm([c 4 1 X.X3_12],1,B');
    lmiterm([c 4 1 -X.X1_11],C,1);
    lmiterm([c 4 1 -X.X1_12],D,1);
	lmiterm([c 4 2 -X.X1_21],C,1);
	lmiterm([c 4 2 -X.X1_22],D,1);
    lmiterm([c 4 3 -X.X2_11],C,1);
    lmiterm([c 4 3 -X.X2_12],D,1);
    lmiterm([c 4 4 0],-1);
    lmiterm([c 4 4 -X.X3_11],C,1,'s');
    lmiterm([c 4 4 -X.X3_12],D,1,'s');
    
    lmiterm([c 5 1 S],1,1);
    lmiterm([c 5 2 L],1,1);
    lmiterm([c 6 1 -L],1,1);
    lmiterm([c 6 2 Z],1,1);
    lmiterm([c 5 1 X.X4_11],1,A');
    lmiterm([c 5 1 X.X4_12],1,B');
    lmiterm([c 6 1 X.X4_21],1,A');
    lmiterm([c 6 1 X.X4_22],1,B');
    lmiterm([c 5 1 -X.X1_11],-1,1);
    lmiterm([c 6 1 -X.X1_12],-1,1);
    lmiterm([c 5 2 -X.X1_21],-1,1);
    lmiterm([c 6 2 -X.X1_22],-1,1);
    lmiterm([c 5 3 -X.X2_11],-1,1);
    lmiterm([c 6 3 -X.X2_12],-1,1);
    lmiterm([c 5 4 -X.X3_11],-1,1);
    lmiterm([c 6 4 -X.X3_12],-1,1);
    lmiterm([c 5 4 X.X4_11],1,C');
    lmiterm([c 5 4 X.X4_12],1,D');
    lmiterm([c 6 4 X.X4_21],1,C');
    lmiterm([c 6 4 X.X4_22],1,D');
    
    
    lmiterm([c 5 5 -X.X4_11],-1,1,'s');
    lmiterm([c 6 5 -X.X4_12],-1,1);
    lmiterm([c 6 5 X.X4_21],-1,1);
    lmiterm([c 6 6 -X.X4_22],-1,1,'s');
end