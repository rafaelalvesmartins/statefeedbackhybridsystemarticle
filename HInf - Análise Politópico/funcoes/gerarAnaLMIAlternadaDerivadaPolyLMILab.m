function gerarAnaLMIAlternadaDerivadaPolyLMILab(ACalPoly,ECalPoly,CPoly,Y,X,mu,c)
     
    lmiterm([c 1 1 X{1}],1,ACalPoly','s');
    lmiterm([c 2 1 X{2}],1,ACalPoly');
    lmiterm([c 2 1 0],ECalPoly');
    lmiterm([c 2 2 mu],-1,1);
    lmiterm([c 3 1 X{3}],1,ACalPoly');
    lmiterm([c 3 1 -X{1}],CPoly,1);
    lmiterm([c 3 2 -X{2}],CPoly,1);
    lmiterm([c 3 3 0],-1);
    lmiterm([c 3 3 X{3}],1,CPoly','s');
    lmiterm([c 4 1 Y],1,1);
    lmiterm([c 4 1 X{4}],1,ACalPoly');
    lmiterm([c 4 1 -X{1}],-1,1);
    lmiterm([c 4 2 -X{2}],-1,1);
    lmiterm([c 4 3 X{4}],1,CPoly');
    lmiterm([c 4 3 -X{3}],-1,1);
    lmiterm([c 4 4 X{4}],-1,1,'s');
    
end