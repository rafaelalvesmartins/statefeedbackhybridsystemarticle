function  TLMI = gerarSinLMIAlternadaDerivada(A0,B0,E0,C0,D0,Alphas,Betas,Gammas,Deltas,Epsilons,S_,T_,U_,nx,nw,ny,S,Z,L,derS,derL,derZ,mu,i)
                                              
    %   S_Vec_nx
    S_Vec_nx = gerarVecRepInt(nx,S{i});
    S_Vec_nx = S_Vec_nx';
    
    %   L_TVec_nx
    L_TVec_nx = gerarVecRepInt(nx,L{i});
    L_TVec_nx = L_TVec_nx';
    
    %   L_Vec_nx
    L_Vec_nx = gerarVecRepInt(nx,L{i}');
    L_Vec_nx = L_Vec_nx';
    
     %   Z_Vec_nx
    Z_Vec_nx = gerarVecRepInt(nx,Z{i});
    Z_Vec_nx = Z_Vec_nx';
    
    %   SVec_ny
    SVec_ny = gerarVecRepInt(ny,S{i});
    SVec_ny = SVec_ny';
    
    %   L_TVec_ny
    L_TVec_ny = gerarVecRepInt(ny,L{i});
    L_TVec_ny = L_TVec_ny';
    
    %   L_Vec_ny
    L_Vec_ny = gerarVecRepInt(ny,L{i}');
    L_Vec_ny = L_Vec_ny';
    
    %   Z_Vec_ny
    Z_Vec_ny = gerarVecRepInt(ny,Z{i});
    Z_Vec_ny = Z_Vec_ny';
    
    %   I_nwVec_nx
    I_nwVec_nx = gerarVecRepInt(nw,eye(nx));
    I_nwVec_nx = I_nwVec_nx';
    
    U11 = -derS + A0*S{i}+B0*L{i}'+S{i}*A0'+L{i}*B0'+S_;
    
    U21 = -derL' + L{i}'*A0'+Z{i}'*B0';
    U22 = -derZ;
    
    U31 = E0';
    U32 = zeros(size(U31,1),size(U22,2));
    U33 = -mu*eye(size(U31,1),size(U31,1))+T_;
    
    U41 = C0*S{i}+D0*L{i}';
    U42 = C0*L{i}+D0*Z{i};
    U43 = zeros(size(U41,1),size(U33,2));
    U44 =  -eye(size(U41,1),size(U41,1))+U_;
    
    U51 = S_Vec_nx;
    U52 = zeros(size(U51,1),size(U42,2));
    U53 = zeros(size(U51,1),size(U43,2));
    U54 = zeros(size(U51,1),size(U44,2));
    U55 = -Alphas;
    
    U61 = L_TVec_nx;
    U62 = zeros(size(U61,1),size(U52,2));
    U63 = zeros(size(U61,1),size(U53,2));
    U64 = zeros(size(U61,1),size(U54,2));
    U65 = zeros(size(U61,1),size(U55,2));
    U66 = -Betas;
    
    U71 = zeros(size(Alphas,1),size(U61,2));
    U72 = L_Vec_nx;
    U73 = zeros(size(U71,1),size(U63,2));
    U74 = zeros(size(U71,1),size(U64,2));
    U75 = zeros(size(U71,1),size(U65,2));
    U76 = zeros(size(U71,1),size(U66,2));
    U77 = -Alphas;
    
    U81 = zeros(size(Betas,1),size(U71,2));
    U82 = Z_Vec_nx;
    U83 = zeros(size(U81,1),size(U73,2));
    U84 = zeros(size(U81,1),size(U74,2));
    U85 = zeros(size(U81,1),size(U75,2));
    U86 = zeros(size(U81,1),size(U76,2));
    U87 = zeros(size(U81,1),size(U77,2));
    U88 = -Betas;
    
    U91 = SVec_ny;
    U92 = zeros(size(U91,1),size(U82,2));
    U93 = zeros(size(U91,1),size(U83,2));
    U94 = zeros(size(U91,1),size(U84,2));
    U95 = zeros(size(U91,1),size(U85,2));
    U96 = zeros(size(U91,1),size(U86,2));
    U97 = zeros(size(U91,1),size(U87,2));
    U98 = zeros(size(U91,1),size(U88,2));
    U99 = -Gammas;
    
    U101 = L_TVec_ny;
    U102 = zeros(size(U101,1),size(U92,2));
    U103 = zeros(size(U101,1),size(U93,2));
    U104 = zeros(size(U101,1),size(U94,2));
    U105 = zeros(size(U101,1),size(U95,2));
    U106 = zeros(size(U101,1),size(U96,2));
    U107 = zeros(size(U101,1),size(U97,2));
    U108 = zeros(size(U101,1),size(U98,2));
    U109 = zeros(size(U101,1),size(U99,2));
    U1010 = -Deltas;
    
    U111 = zeros(size(Gammas,1),size(U101,2));
    U112 = L_Vec_ny;
    U113 = zeros(size(U111,1),size(U103,2));
    U114 = zeros(size(U111,1),size(U104,2));
    U115 = zeros(size(U111,1),size(U105,2));
    U116 = zeros(size(U111,1),size(U106,2));
    U117 = zeros(size(U111,1),size(U107,2));
    U118 = zeros(size(U111,1),size(U108,2));
    U119 = zeros(size(U111,1),size(U109,2));
    U1110 = zeros(size(U111,1),size(U1010,2));
    U1111 = -Gammas;
    
    U121 = zeros(size(Deltas,1),size(U101,2));
    U122 = Z_Vec_ny;
    U123 = zeros(size(U121,1),size(U113,2));
    U124 = zeros(size(U121,1),size(U114,2));
    U125 = zeros(size(U121,1),size(U115,2));
    U126 = zeros(size(U121,1),size(U116,2));
    U127 = zeros(size(U121,1),size(U117,2));
    U128 = zeros(size(U121,1),size(U118,2));
    U129 = zeros(size(U121,1),size(U119,2));
    U1210 = zeros(size(U121,1),size(U1110,2));
    U1211 = zeros(size(U121,1),size(U1111,2));
    U1212 = -Deltas;
    
    U131 = I_nwVec_nx;
    U132 = zeros(size(U131,1),size(U122,2));
    U133 = zeros(size(U131,1),size(U123,2));
    U134 = zeros(size(U131,1),size(U124,2));
    U135 = zeros(size(U131,1),size(U125,2));
    U136 = zeros(size(U131,1),size(U126,2));
    U137 = zeros(size(U131,1),size(U127,2));
    U138 = zeros(size(U131,1),size(U128,2));
    U139 = zeros(size(U131,1),size(U129,2));
    U1310 = zeros(size(U131,1),size(U1210,2));
    U1311 = zeros(size(U131,1),size(U1211,2));
    U1312 = zeros(size(U131,1),size(U1212,2));
    U1313 = -Epsilons;
    
    
    TLMI =  [U11 U21' U31' U41' U51' U61' U71' U81' U91'  U101' U111'  U121'  U131';
             U21 U22  U32' U42' U52' U62' U72' U82' U92'  U102' U112'  U122'  U132';
             U31 U32  U33  U43' U53' U63' U73' U83' U93'  U103' U113'  U123'  U133';
             U41 U42  U43  U44  U54' U64' U74' U84' U94'  U104' U114'  U124'  U134';
             U51 U52  U53  U54  U55  U65' U75' U85' U95'  U105' U115'  U125'  U135';
             U61 U62  U63  U64  U65  U66  U76' U86' U96'  U106' U116'  U126'  U136';
             U71 U72  U73  U74  U75  U76  U77  U87' U97'  U107' U117'  U127'  U137';
             U81 U82  U83  U84  U85  U86  U87  U88  U98'  U108' U118'  U128'  U138';
             U91 U92  U93  U94  U95  U96  U97  U98  U99   U109' U119'  U129'  U139';
             U101 U102 U103 U104 U105 U106 U107 U108 U109 U1010  U1110' U1210' U1310';
             U111 U112 U113 U114 U115 U116 U117 U118 U119 U1110  U1111  U1211' U1311';
             U121 U122 U123 U124 U125 U126 U127 U128 U129 U1210  U1211  U1212  U1312';
             U131 U132 U133 U134 U135 U136 U137 U138 U139 U1310  U1311  U1312  U1313];
         
     
end