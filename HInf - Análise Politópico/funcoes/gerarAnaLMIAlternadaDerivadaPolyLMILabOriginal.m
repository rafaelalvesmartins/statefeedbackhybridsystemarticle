function gerarAnaLMIAlternadaDerivadaPolyLMILabOriginal(ACalPoly,ECalPoly,CPoly,Y,mu,c)
     
    lmiterm([c 1 1 Y],1,ACalPoly','s');
    lmiterm([c 2 1 0],ECalPoly');
    lmiterm([c 2 2 mu],-1,1);
    lmiterm([c 3 1 Y],CPoly,1);
    lmiterm([c 3 3 0],-1);
end