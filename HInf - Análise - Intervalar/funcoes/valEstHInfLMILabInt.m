function saida = valEstHInfLMILabInt(ACal,ECal,CCal,KCal,h,delta,tol)
    % validarLMIHInf
    %
    % Avaliar norma HInf de um sistema híbrido
    %
    % input:  ACal,ECal,CCal,KCal -> Matrizes aumentadas
    %         h -> período de amostragem
    %         delta -> relacionada a aproximação da derivada
    %         tol -> tolerância do solver
    %
    % output: saida.X -> Matriz Lyapunov
    %
    % Criado data: 08 Mar 2020
    % Revisão data: 08 Mar 2020
    % Autor: rafaelmartinsalves@gmail.com
    
    % Matrizes Incertas
    ACal0 = mid(ACal);
    deltaACal = rad(ACal);
    
    ECal0 = mid(ECal);
    deltaECal = rad(ECal);
    
    CCal0 = mid(CCal);
    deltaCCal = rad(CCal);
    
    % Dimensões das matrizes de entrada
    nXi = length(ACal.inf);
    nW = size(ECal.inf, 2);
    nY = size(CCal.inf, 1);
    
    % Inicialização das variáveis LMIs
    saida.cpuSeg_m = clock;
    %LMI contador de linhas
    saida.linhas = 0;
    %Inicia o sistema de LMI
    setlmis([]);
	%Configura contador LMI
    c = 1;
    
    % Cria as matrizes X_i     
    quantInt = round(h/delta);
    X = cell(quantInt,1);
    sX = cell(quantInt,1);
    for i=1:quantInt
		[X{i},nVar,sX{i}] = lmivar(1,[nXi 1]);
		lmiterm([-c 1 1 X{i}],1,1); % X{i} > 0
		c = c + 1;
    end
    
    % Função objetiva
	mu = lmivar(1,[1 0]);
    lmiterm([-c 1 1 mu],1,1); % mu > 0
    c = c + 1;
    obj = mu;
    
    %   Alphas 1
    for i=1:nXi
        for j=1:nXi
            alphas(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end
    
    
    %    Gammas
    for i=1:nY
        for j=1:nXi
            gammas(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end
    
    %    Epsilons
    for i=1:nXi
        for j=1:nW
            epsilons(i,j) = lmivar(1,[1 0]); % Scalar
        end
    end
	 
    
    % Montar as restrições LMIs
    for i=1:(quantInt-1)
        
       % PRIMEIRA DERIVADA
        
		% Entrada (1,1)
		lmiterm([c 1 1 X{i+1}],1/delta,1);
        lmiterm([c 1 1 X{i}],-1/delta,1);
        
		gerarLMIAlternadaDerivadaLMILabInt(X,sX,ACal0,CCal0,ECal0,alphas,gammas,epsilons,deltaACal,deltaCCal,deltaECal,nXi,nY,nW,i,mu,c);
	    c = c + 1;
        
        % SEGUNDA DERIVADA
        
		% Entrada (1,1)
		lmiterm([c 1 1 X{i+1}],1/delta,1);
        lmiterm([c 1 1 X{i}],-1/delta,1);
		gerarLMIAlternadaDerivadaLMILabInt(X,sX,ACal0,CCal0,ECal0,alphas,gammas,epsilons,deltaACal,deltaCCal,deltaECal,nXi,nY,nW,i+1,mu,c);
	    c = c + 1;
        
    end
    
        
   
    
    % Restrição para garantir que a Lyapunov sempre decai
	% Entrada (1,1)
	lmiterm([-c 1 1 X{quantInt}], 1, 1);
	% Entrada (2,1)
	lmiterm([-c 2 1 X{1}], 1, KCal);
	% Entrada (2,2)
	lmiterm([-c 2 2 X{1}], 1, 1);
    c = c + 1;

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
            saida.X{i} = dec2mat(LMIs,xopt,X{i});
        end
        saida.mu = dec2mat(LMIs,xopt,mu);
        saida.gamma = sqrt(saida.mu);
    end
    

end