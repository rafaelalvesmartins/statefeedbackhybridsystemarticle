function TLMI = gerarAnaLMIAlternadaDerivadaPoly(ACalPoly,ECalPoly,CPoly,der,Y,X,mu)

	nW = size(ECalPoly,2);
	nY = size(CPoly,1);
    U11 = der + ACalPoly'*X{1}'+X{1}*ACalPoly;
    U21 = X{2}*ACalPoly+ECalPoly'*X{1}';
	U22 = -mu*eye(nW)+X{2}*ECalPoly+ECalPoly'*X{2}';
	U31 = CPoly+X{3}*ACalPoly;
	U32 = X{3}*ECalPoly;
	U33 = -eye(nY);
	U41 = Y+X{4}*ACalPoly-X{1}';
	U42 = X{4}*ECalPoly-X{2}';
	U43 = -X{3}';
	U44 = -X{4}-X{4}';
	
	TLMI =  [U11 U21' U31' U41';
			 U21 U22  U32' U42';
			 U31 U32  U33  U43';
			 U41 U42  U43  U44];
end