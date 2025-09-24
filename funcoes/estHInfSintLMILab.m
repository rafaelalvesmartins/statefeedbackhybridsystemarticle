function saida = estHInfSintLMILab(A,B,E,C,D,h,delta,tol)
    % valEstHInfSint
    %
    % Sintase Ganho K com custo HInf garantido de um sistema híbrido
    % intervalar
    %
    % input:  A,B,E,C,D -> Matrizes do sistema
    %         h -> período de amostragem
    %         delta -> relacionada a aproximação da derivada
    %         tol -> tolerância do solver
    %
    % output: saida.S{i},saida.Z{i},saida.L{i},saida.W -> Matrizes Lyapunov
    %         saida.K -> Ganho de alimentação de estado
    %         saida.gamma -> Custo garantido HInf
    % Criado data: 29 Mar 2020
    % Revisão data: 29 Mar 2020
    % Autor: rafaelalves@dca.fee.unicamp.br;msouza@fee.unicamp.br
    
     % Matrizes Incertas
    A0 = mid(A);
    deltaA = rad(A);
    
    B0 = mid(B);
    deltaB = rad(B);
    
    E0 = mid(E);
    deltaE = rad(E);
    
    C0 = mid(C);
    deltaC = rad(C);
    
    D0 = mid(D);
    deltaD = rad(D);
    
%     Tamanho das matrizes
    nx = length(A.inf);
    nu = size(B.inf,2);
    nw = size(E.inf,2);
    ny = size(C.inf,1);
    
    % Inicialização das variáveis LMIs
    %Inicia o sistema de LMI
    setlmis([]);
	%Configura contador LMI
    c = 1;
    
    % Cria as matrizes S_i,L_i,Z_i     
    quantInt = round(h/delta);
    S = cell(quantInt,1);
    sS = cell(quantInt,1);
    
    L = cell(quantInt,1);
    sL = cell(quantInt,1);
    
    Z = cell(quantInt,1);
    sZ = cell(quantInt,1);
    for i=1:quantInt
		[S{i},nVar,sS{i}] = lmivar(1,[nx 1]);
        [L{i},nVar,sL{i}] = lmivar(2,[nx nu]);
        [Z{i},nVar,sZ{i}] = lmivar(1,[nu 1]);
    end
    
    % Função objetiva
	mu = lmivar(1,[1 0]);
    lmiterm([-c 1 1 mu],1,1); % mu > 0
    c = c + 1;
    obj = mu;
    
    %   alphas
    for i=1:nx
        for j=1:nx
            alphas(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end
    
    %   betas
    for i=1:nx
        for j=1:nu
            betas(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end
    
    %   gammas
    for i=1:ny
        for j=1:nx
            gammas(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end
    
    %   deltas
    for i=1:ny
        for j=1:nu
            deltas(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end
    
    
    %   epsilons
    for i=1:nx
        for j=1:nw
            epsilons(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end
    
    % Montar as restrições LMIs
    % Primeira restrição
    for i=1:(quantInt-1)
        
        % Para i
        % Derivadas
		lmiterm([c 1 1 S{i+1}],-1/delta,1);
        lmiterm([c 1 1 S{i}],1/delta,1);
        lmiterm([c 2 1 -L{i+1}],-1/delta,1);
        lmiterm([c 2 1 -L{i}],1/delta,1);
        lmiterm([c 2 2 Z{i+1}],-1/delta,1);
        lmiterm([c 2 2 Z{i}],1/delta,1);
        
		gerarSinLMIAlternadaDerivadaLMILab(A0,B0,E0,C0,D0,deltaA,deltaB,deltaE,deltaC,deltaD,alphas,betas,gammas,deltas,epsilons,nx,nu,nw,ny,sS,sZ,sL,S,Z,L,mu,i,c);
	    c = c + 1;
        
        % Para i+1
        
		% Primeira restrição
        % Derivadas
		lmiterm([c 1 1 S{i+1}],-1/delta,1);
        lmiterm([c 1 1 S{i}],1/delta,1);
        lmiterm([c 2 1 -L{i+1}],-1/delta,1);
        lmiterm([c 2 1 -L{i}],1/delta,1);
        lmiterm([c 2 2 Z{i+1}],-1/delta,1);
        lmiterm([c 2 2 Z{i}],1/delta,1);
        
		gerarSinLMIAlternadaDerivadaLMILab(A0,B0,E0,C0,D0,deltaA,deltaB,deltaE,deltaC,deltaD,alphas,betas,gammas,deltas,epsilons,nx,nu,nw,ny,sS,sZ,sL,S,Z,L,mu,i+1,c);
	    c = c + 1;
        
    end
    
    % Segunda restrição
    % Restrição para garantir que a Lyapunov sempre decai
    W = lmivar(1,[nx 1]);
    KChap = lmivar(2,[nu nx]);
    lmiterm([-c 1 1 S{quantInt}],1,1);
	lmiterm([-c 2 1 -L{quantInt}],1,1);
	lmiterm([-c 2 2 Z{quantInt}],1,1);
    lmiterm([-c 3 1 S{quantInt}],1,1);
    lmiterm([-c 3 2 L{quantInt}],1,1);
    lmiterm([-c 3 3 W],1,1);
	c = c + 1;
	lmiterm([-c 1 1 W],1,1);
    lmiterm([-c 2 1 W],1,1);
    lmiterm([-c 2 2 S{1}],1,1);
    lmiterm([-c 3 1 KChap],1,1);
    lmiterm([-c 3 2 -L{1}],1,1);
    lmiterm([-c 3 3 Z{1}],1,1);
    c = c + 1;
    
    % Terceira restrição
     for i=1:quantInt
		lmiterm([-c 1 1 S{i}],1,1); % S{i} > 0
        lmiterm([-c 2 1 -L{i}],1,1); 
        lmiterm([-c 2 2 Z{i}],1,1); % Z{i} > 0
		c = c + 1;
    end
        
    LMIs = getlmis;

    c = zeros(decnbr(LMIs),1);
    c(diag(decinfo(LMIs,obj))) = 1;
    options = [tol,2000,0,200,1];
    [copt,xopt] = mincx(LMIs,c,options);

    
    if isempty(xopt)
       saida.feas = 0;
    else
        saida.factivel = 1;
        for i=1:quantInt
            saida.S{i} = dec2mat(LMIs,xopt,S{i});
            saida.Z{i} = dec2mat(LMIs,xopt,Z{i});
            saida.L{i} = dec2mat(LMIs,xopt,L{i});
        end
        saida.mu = dec2mat(LMIs,xopt,mu);
        saida.gamma = sqrt(saida.mu);
        saida.W = dec2mat(LMIs,xopt,W);
        saida.KChap = dec2mat(LMIs,xopt,KChap);
        saida.K = saida.KChap*inv(saida.W);
    end
    

end