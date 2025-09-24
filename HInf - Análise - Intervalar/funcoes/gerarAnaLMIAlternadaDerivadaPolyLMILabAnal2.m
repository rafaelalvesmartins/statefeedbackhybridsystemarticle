function gerarAnaLMIAlternadaDerivadaPolyLMILabAnal2(ACalPoly,ECalPoly,CPoly,X,Y,mu,c)
                                                 
     
    lmiterm([c 1 1 Y{1}],1,ACalPoly,'s');
    lmiterm([c 2 1 Y{2}],1,ACalPoly);
    lmiterm([c 2 1 -Y{1}],ECalPoly',1);
    lmiterm([c 2 2 mu],-1,1);
    lmiterm([c 2 2 Y{2}],1,ECalPoly,'s');
    lmiterm([c 3 1 0],CPoly);
    lmiterm([c 3 1 Y{3}],1,ACalPoly);
    lmiterm([c 3 2 Y{3}],1,ECalPoly);
    lmiterm([c 3 3 0],-1);
    lmiterm([c 4 1 X],1,1);
    lmiterm([c 4 1 Y{4}],1,ACalPoly);
    lmiterm([c 4 1 -Y{1}],-1,1);
    lmiterm([c 4 2 Y{4}],1,ECalPoly);
    lmiterm([c 4 2 -Y{2}],-1,1);
    lmiterm([c 4 3 -Y{3}],-1,1);
    lmiterm([c 4 4 Y{4}],-1,1,'s');
    
end