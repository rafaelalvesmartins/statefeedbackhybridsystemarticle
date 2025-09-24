function gerarLMIAlternadaDerivadaLMILabSemInt(ACal,ECal,CCal,X,mu,i,c)

    lmiterm([c 1 1 X{i}],1,ACal,'s');
    lmiterm([c 2 1 X{i}],ECal',1);
	lmiterm([c 2 2 mu],-1,1);
	lmiterm([c 3 1 0],CCal);
	lmiterm([c 3 3 0],-1);
	
end